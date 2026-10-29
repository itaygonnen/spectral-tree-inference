"""A third data source: IQ-TREE random trees + Sailfish/msasim alignments.

    iqtree -r N -rlen MIN MEAN MAX   ->  <name>/newick/random_tree_0001.nwk
    msasim JC + Gamma, no indels     ->  <name>/fasta/random_tree_0001.fasta

That layout is the whole integration: ``real_datasets.list_datasets`` claims any directory
holding both a ``fasta/`` and a ``newick/`` and pairs them by stem, so what is written here
is read by ``RealLoader`` and swept by the code the downloaded datasets already use.

**Indels are off by default and that is load-bearing.** ``JC_similarity_matrix`` takes its
alphabet size from the data (``spectraltree/similarities.py:39``), so one gap character
makes k=5 and the correction is wrong. At zero rates msasim takes an explicit no-indel
path (``msasim/protocol.py:43``) and every sequence is exactly ``root_seq_size`` long over
{A,C,G,T}. ``scripts/generate_dataset.py`` is the CLI and the gatekeeper for that.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, List, Optional, Tuple

# IQ-TREE 3 dropped ``-r`` from ``-h`` in favour of AliSim's ``-t RANDOM{model,n}``, but
# the old flag still generates a Yule-Harding tree and still honours ``-rlen`` (verified on
# 3.1.4). AliSim would insist on simulating an alignment too, which is msasim's job here.
IQTREE_CANDIDATES = ("iqtree3", "iqtree2", "iqtree")

STEM_FMT = "random_tree_{:04d}"


@dataclass(frozen=True)
class SailfishSpec:
    """Everything that changes the trees or the alignments, and nothing else.

    Every field lands in the dataset name or in ``manifest.json``, because two generations
    that differ in any of them must not share a screen cache -- the same reasoning as
    :func:`src.runners.generated_source.source_name`.
    """

    n_taxa: int
    n_trees: int = 10
    root_seq_size: int = 1000
    rlen: Tuple[float, float, float] = (0.001, 0.1, 0.5)
    gamma_alpha: float = 1.0
    gamma_categories: int = 4
    insertion_rate: float = 0.0
    deletion_rate: float = 0.0
    zipf_alpha: float = 1.7
    zipf_maxlen: int = 50
    tree_seed_base: int = 1000
    msa_seed_base: int = 90000
    iqtree_bin: Optional[str] = None
    # Rejection-sampled to the eta buckets rather than drawn blind. It belongs in the
    # name because such a population is CONDITIONED: sharing a name -- and so a screen
    # cache -- with an unconditioned one would mix two different distributions.
    pooled: bool = False

    def __post_init__(self):
        if self.n_taxa < 4:
            raise ValueError(f"n_taxa must be >= 4, got {self.n_taxa}")
        if self.n_trees < 1:
            raise ValueError(f"n_trees must be >= 1, got {self.n_trees}")
        if self.root_seq_size < 1:
            raise ValueError(f"root_seq_size must be >= 1, got {self.root_seq_size}")
        lo, mean, hi = self.rlen
        if not lo <= mean <= hi:
            raise ValueError(f"rlen must be (min <= mean <= max), got {self.rlen}")
        if min(self.insertion_rate, self.deletion_rate) < 0:
            raise ValueError("indel rates must be non-negative")
        if self.gamma_categories < 1:
            raise ValueError("gamma_categories must be >= 1")

    @property
    def has_indels(self) -> bool:
        return bool(self.insertion_rate or self.deletion_rate)

    @property
    def name(self) -> str:
        """Dataset directory name -- the parameters that make two runs incomparable."""
        parts = [f"sailfish_yule_n{self.n_taxa}",
                 f"L{self.root_seq_size}",
                 f"g{self.gamma_alpha:g}x{self.gamma_categories}"]
        if self.has_indels:
            parts.append(f"indel{self.insertion_rate:g}-{self.deletion_rate:g}")
        if self.pooled:
            parts.append("buckets")
        return "_".join(parts)

    def stems(self) -> List[str]:
        return [STEM_FMT.format(i) for i in range(1, self.n_trees + 1)]

    def seeds_for(self, index: int) -> Tuple[int, int]:
        """The ``(tree, msa)`` seed pair the blind path uses for tree ``index``."""
        return self.tree_seed_base + index, self.msa_seed_base + index


def find_iqtree(spec: SailfishSpec) -> str:
    """The IQ-TREE binary to shell out to, or a message naming how to get one."""
    if spec.iqtree_bin:
        if shutil.which(spec.iqtree_bin) is None:
            raise FileNotFoundError(f"iqtree binary not found: {spec.iqtree_bin!r}")
        return spec.iqtree_bin
    for cand in IQTREE_CANDIDATES:
        if shutil.which(cand):
            return cand
    raise FileNotFoundError(
        "no IQ-TREE binary on PATH (looked for " + ", ".join(IQTREE_CANDIDATES) +
        "). Install one with `brew install iqtree3`, or pass iqtree_bin=.")


def iqtree_version(binary: str) -> str:
    out = subprocess.run([binary, "--version"], capture_output=True, text=True)
    first = (out.stdout or out.stderr).strip().splitlines()
    return first[0] if first else binary


def preflight(spec: Optional["SailfishSpec"] = None) -> List[str]:
    """What is missing before a run can start, as actionable lines (empty = ready).

    Checked up front rather than discovered per attempt: ``attempt`` catches its own
    exceptions so one bad seed cannot end a run, which means a missing ``msasim`` would
    not raise -- it would quietly fail every one of ``--max-attempts`` tries and report
    an empty dataset after burning the whole job.
    """
    problems: List[str] = []
    try:
        import msasim                                      # noqa: F401
    except ImportError:
        problems.append("msasim is not installed   ->  pip install msasim "
                        "(wheels exist for macOS arm64 and manylinux x86_64)")
    try:
        find_iqtree(spec or SailfishSpec(n_taxa=4))
    except FileNotFoundError as exc:
        problems.append(f"{exc}\n      on a cluster: `module load iqtree` or a conda "
                        f"env, then pass --iqtree-bin")
    return problems


def msasim_version() -> str:
    """msasim exposes no ``__version__``; the installed distribution's metadata does."""
    from importlib.metadata import PackageNotFoundError, version
    try:
        return version("msasim")
    except PackageNotFoundError:
        return "unknown"


def make_tree(spec: SailfishSpec, seed: int, stem: str, newick_dir: Path,
              binary: Optional[str] = None) -> Path:
    """One random Yule-Harding tree with ``-rlen`` branch lengths.

    ``seed`` and ``stem`` are separate arguments because rejection sampling needs them to
    be: the seed advances on every *attempt*, while the name is only chosen once a tree is
    *accepted*. The blind path passes ``spec.seeds_for(index)`` and the index's stem, which
    is what it always did.
    """
    binary = binary or find_iqtree(spec)
    newick_dir.mkdir(parents=True, exist_ok=True)
    dest = newick_dir / f"{stem}.nwk"
    lo, mean, hi = spec.rlen
    proc = subprocess.run(
        [binary, "-redo", "-r", str(spec.n_taxa),
         "-rlen", f"{lo:g}", f"{mean:g}", f"{hi:g}",
         "-seed", str(int(seed)), str(dest)],
        capture_output=True, text=True)
    # iqtree leaves a sibling .log next to the tree; it is not a newick and nothing reads it
    log = dest.with_suffix(".nwk.log")
    if log.exists():
        log.unlink()
    if proc.returncode != 0 or not dest.exists() or dest.stat().st_size == 0:
        raise RuntimeError(
            f"iqtree failed for {stem} at seed {seed} (exit {proc.returncode}):\n"
            f"{(proc.stderr or proc.stdout)[-2000:]}")
    return dest


def make_msa(spec: SailfishSpec, tree_path: Path, seed: int, fasta_dir: Path) -> Path:
    """One alignment on ``tree_path``: JC + Gamma, indels per the spec (zero by default).

    In-process rather than a subprocess -- msasim is a Python package, which is what lets
    this drop the environment-activation half of the cluster scripts entirely.
    """
    from msasim import (ALPHABET_CODES, MODEL_CODES, ReplacementModelSpec,
                        SimProtocol, Simulator, SiteRateModelSpec, ZipfDistribution)

    fasta_dir.mkdir(parents=True, exist_ok=True)
    dest = fasta_dir / f"{tree_path.stem}.fasta"

    protocol = SimProtocol(
        tree=str(tree_path),
        root_seq_size=spec.root_seq_size,
        insertion_rate=spec.insertion_rate,
        deletion_rate=spec.deletion_rate,
        insertion_dist=ZipfDistribution(spec.zipf_alpha, spec.zipf_maxlen),
        deletion_dist=ZipfDistribution(spec.zipf_alpha, spec.zipf_maxlen),
        seed=int(seed),
    )
    replacement = ReplacementModelSpec(
        model=MODEL_CODES.NUCJC,
        alphabet=ALPHABET_CODES.DNA,
        site_rate_model=SiteRateModelSpec(gamma_alpha=spec.gamma_alpha,
                                          gamma_categories=spec.gamma_categories),
    )
    Simulator(protocol, replacement_model=replacement)().write_msa(str(dest))

    if not dest.exists() or dest.stat().st_size == 0:
        raise RuntimeError(f"msasim wrote no alignment for {tree_path.name}")
    return dest


def _done(newick_dir: Path, fasta_dir: Path, stem: str) -> bool:
    nwk, fa = newick_dir / f"{stem}.nwk", fasta_dir / f"{stem}.fasta"
    return (nwk.exists() and nwk.stat().st_size > 0
            and fa.exists() and fa.stat().st_size > 0)


def generate(spec: SailfishSpec, out_root: Path, name: Optional[str] = None,
             log: Callable[[str], None] = print) -> Path:
    """Write ``<out_root>/<name>/{fasta,newick}`` and its manifest; return the directory.

    Resumable: an index whose tree and alignment are both already on disk is skipped, so an
    interrupted generation is finished by re-running the same command.
    """
    binary = find_iqtree(spec)
    dataset = Path(out_root) / (name or spec.name)
    newick_dir, fasta_dir = dataset / "newick", dataset / "fasta"
    started = time.time()

    log(f"iqtree : {iqtree_version(binary)}")
    log(f"msasim : {msasim_version()}")
    log(f"into   : {dataset}")

    built = skipped = 0
    for i in range(1, spec.n_trees + 1):
        stem = STEM_FMT.format(i)
        if _done(newick_dir, fasta_dir, stem):
            skipped += 1
            log(f"  {stem}  (already on disk)")
            continue
        t0 = time.time()
        tree_seed, msa_seed = spec.seeds_for(i)
        tree = make_tree(spec, tree_seed, stem, newick_dir, binary)
        make_msa(spec, tree, msa_seed, fasta_dir)
        built += 1
        log(f"  {stem}  {time.time() - t0:5.1f}s")

    write_manifest(dataset, spec, binary, elapsed=time.time() - started)
    log(f"done   : {built} built, {skipped} already there, "
        f"{time.time() - started:.1f}s")
    return dataset


def write_manifest(dataset: Path, spec: SailfishSpec, binary: str,
                   elapsed: float = 0.0, fill: Optional[dict] = None) -> Path:
    """Record every parameter beside the data, so a dataset explains itself later.

    ``fill`` carries the rejection-sampling terms when there were any. Its ``last_seed``
    is what a resumed run continues from, and its presence is how a reader -- or a later
    run -- knows this population is conditioned on eta rather than drawn blind.
    """
    dest = dataset / "manifest.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    body = {
        "generator": "sailfish_sim",
        "spec": asdict(spec),
        "iqtree": iqtree_version(binary),
        "iqtree_bin": binary,
        "msasim": msasim_version(),
        "written": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_sec": round(elapsed, 1),
    }
    if fill:
        body["fill"] = fill
    dest.write_text(json.dumps(body, indent=2, default=str) + "\n")
    return dest
