"""Unified cache + results I/O for sub_sampled_fielder_vec.

One key format, one ``.complete`` sentinel, one set of scopes for every
intermediate cache and every final result location in the repo.

Disk layout::

    sub_sampled_fielder_vec/
      cache/
        full_matrix/<key>/...
        experiment_data/<key>/...
        pool_sample/<param_key>/eta<TT>/sample_NNNN/...
        sweep_trial/<full_key>/p<p>_seed<seed>/...
        bootstrap_sweep/<sweep_id>/...
      results/
        runs/<timestamp>-<run_name>/...
        notebooks/<subdir>/<notebook>/...

Atomicity: every ``save`` writes payload files first and then touches a
``.complete`` sentinel last. Loaders check the sentinel before reading, so
an interrupted write behaves like a cache miss and is recomputed cleanly
on the next call.

Key format: floats are ``f"{value:.4f}".replace(".", "p")``; ints are
zero-padded to 4 digits. Kwargs are sorted before joining so call-site
argument order doesn't matter. Long parameter sets (>12 keys) get a
``sha1[:12]`` suffix to bound filename length.
"""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path
from typing import Any, Callable, Optional

import numpy as np


PROJECT_ROOT = Path(__file__).parent.parent  # src/cache_io.py -> sub_sampled_fielder_vec/
CACHE_ROOT   = PROJECT_ROOT / "cache"
RESULTS_ROOT = PROJECT_ROOT / "results"
SENTINEL     = ".complete"

_MAX_KEY_KWARGS = 12


def _fmt_value(value: Any) -> str:
    """Filesystem-safe stringification: dots become 'p', floats get 4 decimals."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.4f}".replace(".", "p")
    if isinstance(value, int):
        return f"{value:04d}"
    if isinstance(value, (list, tuple)):
        return "_".join(_fmt_value(v) for v in value)
    return str(value).replace(".", "p").replace("/", "_")


def make_key(prefix: Optional[str] = None, **kwargs: Any) -> str:
    """Deterministic, filesystem-safe key from a prefix and sorted kwargs.

    >>> make_key("balanced_binary", alpha=0.9, n=512)
    'balanced_binary_alpha0p9000_n0512'
    >>> make_key(alpha=0.9, n=512) == make_key(n=512, alpha=0.9)
    True
    """
    parts: list[str] = []
    if prefix:
        parts.append(str(prefix))
    sorted_keys = sorted(kwargs)
    for k in sorted_keys:
        parts.append(f"{k}{_fmt_value(kwargs[k])}")
    base = "_".join(parts)
    if len(sorted_keys) > _MAX_KEY_KWARGS:
        digest = hashlib.sha1(
            json.dumps(kwargs, sort_keys=True, default=str).encode("utf-8")
        ).hexdigest()[:12]
        base = f"{base}__{digest}"
    return base


class CacheScope:
    """A namespace under ``CACHE_ROOT``.

    Each scope owns one subdirectory tree. Entries are addressed by a
    primary ``key`` plus optional ``*subkeys`` for nested layouts (e.g.
    per-trial under a per-experiment dir).

    Artifacts inside an entry are written with one of:
      - ``np.ndarray``        → ``<name>.npz`` keyed by ``<name>``
      - ``dict``/``list``/scalar → ``<name>.json``
      - ``str``               → ``<name>.txt``
    """

    def __init__(self, name: str, root: Optional[Path] = None) -> None:
        self.name = name
        self._root = (root or CACHE_ROOT) / name

    @property
    def root(self) -> Path:
        return self._root

    def path(self, key: str, *subkeys: str) -> Path:
        d = self._root / key
        for sk in subkeys:
            d = d / sk
        return d

    def is_complete(self, key: str, *subkeys: str) -> bool:
        return (self.path(key, *subkeys) / SENTINEL).exists()

    def save(self, key: str, *subkeys: str, **artifacts: Any) -> Path:
        """Write artifacts atomically. ``.complete`` sentinel touched last."""
        d = self.path(key, *subkeys)
        d.mkdir(parents=True, exist_ok=True)
        # Drop any pre-existing sentinel so a partial write between save() and
        # touch can never be observed as complete by a parallel reader.
        sentinel = d / SENTINEL
        if sentinel.exists():
            sentinel.unlink()
        for fname, payload in artifacts.items():
            self._write_one(d, fname, payload)
        sentinel.touch()
        return d

    @staticmethod
    def _write_one(d: Path, fname: str, payload: Any) -> None:
        if isinstance(payload, np.ndarray):
            np.savez_compressed(d / f"{fname}.npz", **{fname: payload})
            return
        if isinstance(payload, str):
            (d / f"{fname}.txt").write_text(payload)
            return
        # JSON for dict/list/scalar/None; default=str catches numpy scalars.
        (d / f"{fname}.json").write_text(
            json.dumps(payload, indent=2, default=str)
        )

    def load(self, key: str, *subkeys: str) -> dict[str, Any]:
        """Return every artifact in the entry. Raises if sentinel missing."""
        d = self.path(key, *subkeys)
        if not (d / SENTINEL).exists():
            raise FileNotFoundError(f"cache miss (no sentinel): {d}")
        out: dict[str, Any] = {}
        for f in sorted(d.iterdir()):
            if f.name == SENTINEL:
                continue
            stem = f.stem
            if f.suffix == ".npz":
                with np.load(f, allow_pickle=True) as nz:
                    nz_keys = list(nz.keys())
                    if len(nz_keys) == 1 and nz_keys[0] == stem:
                        out[stem] = nz[stem]
                    else:
                        out[stem] = {k: nz[k] for k in nz_keys}
            elif f.suffix == ".json":
                out[stem] = json.loads(f.read_text())
            elif f.suffix == ".txt":
                out[stem] = f.read_text()
            # Anything else: silently skip
        return out

    def try_load(self, key: str, *subkeys: str) -> Optional[dict[str, Any]]:
        """``None`` on cache miss, dict on hit. Never raises."""
        if not self.is_complete(key, *subkeys):
            return None
        try:
            return self.load(key, *subkeys)
        except Exception:
            return None

    def get_or_compute(
        self,
        key: str,
        *subkeys: str,
        build_fn: Callable[[], dict[str, Any]],
        force: bool = False,
    ) -> tuple[dict[str, Any], bool]:
        """Return ``(data, hit)``. ``hit=True`` on cache hit; ``False`` after computing."""
        if not force:
            cached = self.try_load(key, *subkeys)
            if cached is not None:
                return cached, True
        artifacts = build_fn()
        self.save(key, *subkeys, **artifacts)
        return artifacts, False

    def extend(self, key: str, *subkeys: str, **artifacts: Any) -> Path:
        """Merge new artifacts into an existing entry; re-stamp sentinel.

        Used by the lazy-metric path: load metrics, compute the missing
        ones, write the merged dict back, re-stamp ``.complete``.
        """
        d = self.path(key, *subkeys)
        d.mkdir(parents=True, exist_ok=True)
        sentinel = d / SENTINEL
        if sentinel.exists():
            sentinel.unlink()
        for fname, payload in artifacts.items():
            self._write_one(d, fname, payload)
        sentinel.touch()
        return d

    def list_keys(self) -> list[str]:
        """Direct-child key names with at least one entry (no completeness check)."""
        if not self._root.exists():
            return []
        return sorted(p.name for p in self._root.iterdir() if p.is_dir())

    def clear(self, key: Optional[str] = None) -> None:
        """Remove one entry or the whole scope."""
        if key is None:
            if self._root.exists():
                shutil.rmtree(self._root)
        else:
            d = self._root / key
            if d.exists():
                shutil.rmtree(d)


# Pre-built scopes. These are the ONLY scopes that should exist on disk -- see
# docs/CACHE_AND_RESULTS.md for each one's size, key scheme, writer and readers.
full_matrix     = CacheScope("full_matrix")
experiment_data = CacheScope("experiment_data")
pool_sample     = CacheScope("pool_sample")
sweep_trial     = CacheScope("sweep_trial")
bootstrap_sweep = CacheScope("bootstrap_sweep")
bpart_sweep     = CacheScope("bpart_sweep")
# Declared here rather than constructed ad hoc in src/runners/nj_sweep.py, which is how
# a seventh cache directory appeared on disk outside the declared set.
distance_matrix = CacheScope("distance_matrix")


def run_dir(
    timestamp: str,
    run_name: str = "",
    *,
    tree_model: Optional[str] = None,
    sampling_method: Optional[str] = None,
) -> Path:
    """Path under ``results/runs/`` for one production experiment run.

    Flat layout (no kwargs): ``results/runs/<ts>-<run_name>/``.

    Nested layout when ``tree_model`` and ``sampling_method`` are both
    given: ``results/runs/<tree_model>/<sampling_method>/<ts>-<run_name>/``
    — the same hierarchy ``ExperimentRunner._make_run_dir`` used to
    create at the top of ``results/`` before unification.
    """
    base = RESULTS_ROOT / "runs"
    if tree_model and sampling_method:
        base = base / tree_model / sampling_method
    suffix = f"-{run_name}" if run_name else ""
    return base / f"{timestamp}{suffix}"


def notebook_dir(notebook_relpath: str) -> Path:
    """Path under ``results/notebooks/`` for a notebook's outputs.

    ``notebook_relpath`` is e.g. ``"04_tree_reconstruction/stdr_partition_recovery"``.
    """
    return RESULTS_ROOT / "notebooks" / notebook_relpath


__all__ = [
    "PROJECT_ROOT", "CACHE_ROOT", "RESULTS_ROOT", "SENTINEL",
    "make_key", "CacheScope",
    "full_matrix", "experiment_data", "pool_sample",
    "sweep_trial", "bootstrap_sweep", "bpart_sweep",
    "run_dir", "notebook_dir",
]


if __name__ == "__main__":
    # Smoke checks. Run with: python -m src.cache_io
    import tempfile

    # 1. Key formatting is order-independent and deterministic.
    a = make_key("balanced_binary", alpha=0.9, n=512)
    b = make_key("balanced_binary", n=512, alpha=0.9)
    assert a == b == "balanced_binary_alpha0p9000_n0512", a

    # 2. Sentinel discipline + round-trip via a temp scope.
    with tempfile.TemporaryDirectory() as tmp:
        scope = CacheScope("smoke", root=Path(tmp))
        key = make_key("test", alpha=0.5, n=4)

        arr = np.arange(4, dtype=np.float32)
        scope.save(key, fiedler=arr, metrics={"sign_agreement": 0.95})
        assert scope.is_complete(key)

        data = scope.load(key)
        assert np.allclose(data["fiedler"], arr)
        assert data["metrics"]["sign_agreement"] == 0.95

        # get_or_compute hits cache the second time.
        def build():
            raise AssertionError("should not be called on cache hit")

        result, hit = scope.get_or_compute(key, build_fn=build)
        assert hit is True

        # extend merges metrics without invalidating fiedler.
        scope.extend(key, metrics={"sign_agreement": 0.95, "ari": 0.88})
        data2 = scope.load(key)
        assert data2["metrics"] == {"sign_agreement": 0.95, "ari": 0.88}
        assert np.allclose(data2["fiedler"], arr)

    print("cache_io smoke checks passed.")
