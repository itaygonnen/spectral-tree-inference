"""Fiedler vector computation class."""
import numpy as np
import scipy.linalg
from scipy.sparse import csr_matrix, diags, issparse
from scipy.sparse.linalg import LinearOperator, eigsh

from ..utils.logging import suppress_warnings, log_warning, log_info, is_progress_mode


class FiedlerVectorComputer:
    """
    Computes Fiedler vectors from Laplacian matrices using optimized dense/sparse eigensolvers.

    The Fiedler vector is the eigenvector corresponding to the second-smallest eigenvalue
    of the graph Laplacian. This class automatically selects the appropriate solver
    (dense or sparse) based on matrix size and density.

    Optimizations:
    - Efficient density estimation (avoids O(N²) counting for large matrices)
    - Optimized sparse solver parameters (sigma=-1e-3, ncv=20)
    - Connectivity checking (warns if graph is disconnected)
    - Consistent sign convention enforcement
    """
    
    def __init__(self, sparsity_threshold: float = 0.9, size_threshold: int = 1024):
        """
        Initialize the Fiedler vector computer.

        Args:
            sparsity_threshold: Minimum sparsity for sparse methods (default 0.9 = 90% zeros)
            size_threshold: Matrix size below which dense methods are preferred (default 1024)
        """
        self.sparsity_threshold = sparsity_threshold
        self.size_threshold = size_threshold
        self.last_method_used = None
        # lambda2 of the last solve, and the threshold under which the graph counts as
        # disconnected -- i.e. under which the Fiedler vector is not a defined quantity.
        # The tolerance is RELATIVE to the shift: sub-sampled entries carry a 1/p factor,
        # so the largest degree can be ~1e5 and an absolute 1e-9 on c - lambda is inside
        # the cancellation error.
        self.last_lambda2 = float('nan')
        self.connected_rtol = 1e-8
        self.use_shifted = True
        self._shift = 0.0
        # below this the dense path is already fast and the Lanczos overhead is not
        # repaid; the gains measured here start at m=6000 (11.6 s -> 0.19 s at p=0.2)
        self.shifted_min_n = 2048
    
    def compute(self, laplacian: np.ndarray, sampling_prob: float = None) -> np.ndarray:
        """
        Compute Fiedler vector from Laplacian matrix.

        Automatically selects dense or sparse solver based on matrix size and density.

        Args:
            laplacian: Laplacian matrix L = D - S (where S is similarity matrix)
            sampling_prob: Optional density hint to avoid O(N²) counting (e.g., if L came from subsampled S with probability p)

        Returns:
            Fiedler vector with consistent sign convention
        """
        n = laplacian.shape[0]

        # Determine whether to use dense or sparse solver
        # Selection logic: Use DENSE if N < size_threshold OR density > 0.1 (sparsity < 0.9)

        # Try the shifted-Lanczos path first: it is 20-120x faster than either solver
        # below at m >= 6000 and returns the identical vector -- BUT only where the
        # Fiedler vector is defined. On a disconnected graph L has one zero eigenvalue
        # per component, so lambda2 = 0 and any vector from that null space is a valid
        # answer; the two solvers then legitimately disagree (measured |cos| = 0.11 at
        # m=6000, p=0.001, which had 10 components). There we defer to the original
        # solver, so nothing already measured moves -- and it is cheap in that regime.
        if self.use_shifted and n >= self.shifted_min_n and self._is_connected(laplacian):
            try:
                vec, lam2 = self._solve_shifted(laplacian)
            except Exception as exc:              # noqa: BLE001 - fall through below
                log_warning('fiedler', f"shifted solve failed ({exc}); using eigh/eigsh")
                vec, lam2 = None, 0.0
            if vec is not None and lam2 > self.connected_rtol * max(self._shift, 1.0):
                self.last_method_used = 'shifted'
                self.last_lambda2 = lam2
                return self._apply_sign_convention(vec)
            self.last_lambda2 = lam2

        if issparse(laplacian):
            # Already in sparse format - use sparse solver
            self.last_method_used = 'sparse'
            return self._solve_sparse_laplacian(laplacian)

        # For dense matrices, decide based on size and density
        if n < self.size_threshold:
            # Small matrix: use dense solver
            self.last_method_used = 'dense'
            return self._solve_dense_laplacian(laplacian)

        # For large matrices, estimate density to decide
        # Optimization: Avoid O(N²) density computation for large matrices
        if sampling_prob is not None:
            # If sampling probability is known, use it as density estimate
            density = sampling_prob
        elif n < 1000:
            # For medium matrices, counting nonzeros is acceptable
            density = np.count_nonzero(laplacian) / laplacian.size
        else:
            # For large matrices, estimate from sample (much faster than full count)
            sample_size = min(100, n)
            sample_indices = np.random.choice(n, sample_size, replace=False)
            sample_density = np.count_nonzero(laplacian[sample_indices]) / (sample_size * n)
            density = sample_density

        # Select solver based on density
        use_dense = density > (1.0 - self.sparsity_threshold)

        if use_dense:
            self.last_method_used = 'dense'
            return self._solve_dense_laplacian(laplacian)
        else:
            self.last_method_used = 'sparse'
            # Convert to sparse format for efficiency
            laplacian_sparse = csr_matrix(laplacian)
            return self._solve_sparse_laplacian(laplacian_sparse)

    def _is_connected(self, laplacian) -> bool:
        """Cheap, exact-enough test that the Fiedler vector is a defined quantity.

        An isolated vertex has degree 0 and contributes its own zero eigenvalue, which is
        how a sub-sampled graph falls apart at low p: at m=1000, p=0.001 there are ~500
        edges for 1000 vertices. O(m) and catches the dominant mode; a split into two
        large components is not caught here but is separated afterwards by the relative
        lambda2 test, and a miss costs only a different (equally valid) vector from a
        degenerate eigenspace, never a wrong one.
        """
        diag = laplacian.diagonal() if issparse(laplacian) else np.diag(laplacian)
        return bool(diag.size) and bool(np.min(diag) > 0)

    def _solve_shifted(self, laplacian, tol: float = 0.0):
        """``(fiedler, lambda2)`` from the LARGEST eigenpairs of ``cI - L``.

        Lanczos converges for the extremal-largest eigenvalues and not for the smallest,
        and a shift needs no factorisation: with ``c >= lambda_max(L)`` (Gershgorin:
        twice the largest degree), the two largest eigenpairs of ``cI - L`` are the two
        smallest of ``L``. Each matvec is one L-apply -- O(m^2) dense, O(nnz) sparse --
        and allocates nothing, where ``scipy.linalg.eigh`` is O(m^3) to use two columns.

        This is the same move that already took the B arm from 14.94 s to 0.12 s
        (:mod:`src.utils.griffing`). Measured here, m=10000 at p=0.2: 60.68 s -> 0.50 s,
        identical vector (|cos| = 1.000000).

        ``tol=0`` asks ARPACK for machine precision rather than a loose tolerance: at
        1e-8 the partition metrics were already bit-identical to the dense solver, but
        the dot product -- which reads the eigenvector's fine structure, not just its
        sign -- moved in the sixth decimal. The extra iterations are cheap next to the
        O(m^3) solve this replaces.

        Returns ``lambda2`` as well, because it is what says whether the answer means
        anything: ``lambda2 == 0`` is a disconnected graph, where "the second smallest
        eigenvector" is an arbitrary pick from a null space of dimension equal to the
        number of components.
        """
        n = laplacian.shape[0]
        diag = laplacian.diagonal() if issparse(laplacian) else np.diag(laplacian)
        c = 2.0 * float(np.max(diag)) if n else 0.0
        self._shift = c
        if c <= 0:
            return None, 0.0
        op = LinearOperator((n, n), dtype=np.float64,
                            matvec=lambda v: c * v - laplacian.dot(v))
        with suppress_warnings('fiedler'):
            vals, vecs = eigsh(op, k=2, which='LA', tol=tol,
                               ncv=min(n - 1, 20), maxiter=n * 20)
        order = np.argsort(-vals)            # descending in cI-L == ascending in L
        lam2 = float(c - vals[order[1]])
        return vecs[:, order[1]], lam2

    def _solve_dense_laplacian(self, laplacian: np.ndarray) -> np.ndarray:
        """
        Solve eigenvalue problem L*v = λ*v using dense solver.

        Uses scipy.linalg.eigh to compute the 2 smallest eigenvalues/eigenvectors.
        Optimal for small or dense matrices.

        Args:
            laplacian: Laplacian matrix L

        Returns:
            Fiedler vector (eigenvector corresponding to 2nd smallest eigenvalue)
        """
        # Suppress numpy/scipy warnings and log in standardized format
        with suppress_warnings('fiedler'):
            # Compute the two smallest eigenvalues and corresponding eigenvectors
            eigvals, eigvecs = scipy.linalg.eigh(laplacian, subset_by_index=(0, 1))

        # The Fiedler vector is the eigenvector corresponding to the second smallest eigenvalue
        fiedler_vector = eigvecs[:, 1]

        # Enforce a consistent sign convention
        return self._apply_sign_convention(fiedler_vector)

    def _solve_sparse_laplacian(self, laplacian) -> np.ndarray:
        """
        Solve eigenvalue problem L*v = λ*v using sparse solver.

        Uses scipy.sparse.linalg.eigsh (ARPACK) with optimized parameters:
        - sigma=-1e-3 for stable shift-invert near zero eigenvalues
        - ncv=20 for efficient Krylov subspace
        - Connectivity check to detect disconnected graphs

        Optimal for large sparse matrices.

        Args:
            laplacian: Sparse Laplacian matrix L

        Returns:
            Fiedler vector (eigenvector corresponding to 2nd smallest eigenvalue)
        """
        n = laplacian.shape[0]
        k = 2  # We need 2 smallest eigenvalues

        # Optimization #4: Reduce ncv from 50 to min(n, 20)
        # Computing 50 Lanczos vectors for just 2 eigenvalues is unnecessary overhead
        ncv = min(n - 1, 20)  # Must satisfy ncv > k and ncv <= n-1

        # Reasonable maxiter based on matrix size
        maxiter = max(n, 1000)

        # Suppress warnings
        with suppress_warnings('fiedler'):
            try:
                # Optimization #2: Change sigma from 1e-10 to -1e-3 for stability
                # A negative shift avoids singularity issues near zero better than a tiny positive one
                eigvals, eigvecs = eigsh(
                    laplacian,
                    k=k,
                    sigma=-1e-3,  # Negative shift for stability
                    ncv=ncv,
                    which='LM',
                    maxiter=maxiter,
                    tol=1e-8  # Slightly relaxed tolerance
                )
            except Exception as e:
                log_warning('fiedler', f"Sparse eigsh with sigma=-1e-3 failed: {str(e)}, trying SM mode...")
                try:
                    # Second attempt: Standard smallest magnitude approach
                    eigvals, eigvecs = eigsh(
                        laplacian,
                        k=k,
                        which='SM',
                        ncv=ncv,
                        maxiter=maxiter,
                        tol=1e-8
                    )
                except Exception as e2:
                    # Final fallback: Use dense computation
                    log_warning('fiedler', f"Sparse computation failed: {str(e2)}, falling back to dense method")
                    # Convert to dense and use dense method
                    laplacian_dense = laplacian.toarray() if hasattr(laplacian, 'toarray') else laplacian
                    return self._solve_dense_laplacian(laplacian_dense)

        # The Fiedler vector is the eigenvector corresponding to the second smallest eigenvalue
        # (eigsh does not guarantee order, so sort)
        idx = np.argsort(eigvals)

        # Optimization #3: Connectivity check
        # If the 2nd eigenvalue is ≈ 0 (< 1e-9), this indicates disconnected graph components
        if eigvals[idx[1]] < 1e-9:
            # Only log in debug mode - expected behavior for low p-values
            if not is_progress_mode():
                log_warning('fiedler',
                    f"Second eigenvalue ({eigvals[idx[1]]:.2e}) is near zero - graph may have disconnected components. "
                    f"Fiedler vector may not represent a valid cut.")

        fiedler_vector = eigvecs[:, idx[1]]

        # Enforce a consistent sign convention
        return self._apply_sign_convention(fiedler_vector)

    def _apply_sign_convention(self, fiedler_vector: np.ndarray) -> np.ndarray:
        # Find first nonzero and exit immediately
        for i, val in enumerate(fiedler_vector):
            if abs(val) > 1e-12:
                return -fiedler_vector if val < 0 else fiedler_vector
        return fiedler_vector  # All zeros case
    
    def align_vector(self, fiedler_vector: np.ndarray, reference_vector: np.ndarray) -> np.ndarray:
        """
        Align a Fiedler vector to have consistent sign orientation with a reference vector.
        
        Alignment is based on sign vectors (binary ±1 vectors), not the actual vectors.
        
        Args:
            fiedler_vector: Fiedler vector to align
            reference_vector: Reference vector for alignment
            
        Returns:
            Aligned Fiedler vector
        """
        # Check for zero or near-zero vectors using a more numerically stable method
        # Use max absolute value as a proxy for checking if vector is non-zero
        ref_max_abs = np.max(np.abs(reference_vector))
        if ref_max_abs < 1e-12:
            raise ValueError("Reference vector has zero or near-zero norm")
        
        fiedler_max_abs = np.max(np.abs(fiedler_vector))
        if fiedler_max_abs < 1e-12:
            log_warning('align', "Fiedler vector has zero or near-zero norm, skipping alignment")
            return fiedler_vector
        
        # Convert to sign vectors (binary ±1 vectors)
        # This is numerically stable as we're just checking signs
        ref_sign = np.sign(reference_vector)
        fiedler_sign = np.sign(fiedler_vector)
        
        # Handle zero entries in sign vectors (should be rare but possible)
        # For exactly zero entries, keep them as 0
        ref_sign[ref_sign == 0] = 1  # Default to positive for zero entries
        fiedler_sign[fiedler_sign == 0] = 1  # Default to positive for zero entries
        
        # Compute dot product of sign vectors
        # Since sign vectors only contain ±1, this is numerically stable
        # Use int64 to prevent overflow for very large vectors
        sign_dot_product = np.dot(fiedler_sign.astype(np.int64), ref_sign.astype(np.int64))
        
        # Flip sign if needed (negative dot product means opposite signs)
        if sign_dot_product < 0:
            return -fiedler_vector
        else:
            return fiedler_vector

