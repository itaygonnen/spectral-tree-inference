"""Rejection-sample Sailfish trees until the eta buckets are full.

Blind generation pays the O(m^3) screen for every tree and then lets ``select_ids`` throw
most of them away -- 49 of 69 valid trees were unbucketed at m=6000. This inverts that:
generate one tree, score it, keep it only if it lands in a bucket that still has room and
its Fiedler partition is a real tree edge, and otherwise delete it and advance the seed.

**Acceptance uses the screen's eta, not the pool's, and that is the whole correctness
argument.** ``scripts/build_eta_pool.py`` accepts on ``estimate_features_from_M`` -- the
sign rule alone -- while the screen scores both cut rules and keeps the more balanced one
(``operators.choose_cut``, ``rule_policy="best"``). A tree admitted on sign-eta can be
re-scored under k-means into a different bucket, or out of every bucket, and the sampling
would be silently pointless. So this calls :func:`~src.runners.operator_screen.screen_one`
and :func:`~src.runners.eta_buckets.bucket_of` -- the same functions ``select_ids`` uses --
and the bucket a tree is accepted into is by construction the bucket it is selected from.

That also makes the row ``screen_one`` returns *the* screen cache row, so an accepted tree
is never screened twice. At m=10000 that is the dominant cost, saved.

The worker/staging/consume split mirrors ``scripts/build_eta_pool_parallel.py``: workers
stage into their own pid-keyed directory, and the main process alone decides capacity and
writes, so there is no index race.
"""
from __future__ import annotations

import os
import shutil
import signal
import sys
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..models.sailfish_sim import SailfishSpec, find_iqtree, make_msa, make_tree
from ..utils.run_paths import screen_cache_path
from .eta_buckets import ETA_TARGETS, ETA_TOL, PER_BUCKET, bucket_of
from .experiment_run import DEFAULT_OPERATORS, gate_operators

_ERROR_ABORT = 10          # consecutive failures that mean the environment, not the seed

_BLAS_VARS = ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")

# set by _init in each worker process
_CTX: Dict = {}


class _Terminated(Exception):
    """SIGTERM, raised in the main loop so it unwinds through the same path as Ctrl+C."""


def _log(msg: str) -> None:
    """Print and flush.

    Python block-buffers stdout when it is not a TTY, which is exactly the cluster case
    (``nohup ... > run.log``). Measured: a run killed after accepting 18 trees had *zero*
    progress lines in its log -- every one was still in the buffer. Flushing costs
    nothing next to an O(m^3) eigensolve.
    """
    print(msg, flush=True)


@contextmanager
def _terminable(log: Callable[[str], None]):
    """Turn SIGTERM into an exception for the duration of the loop.

    SLURM sends SIGTERM at the time limit. The default disposition kills the process
    mid-write, so the manifest and the staging cleanup never run; catching it lets the
    loop finish its checkpoint and tidy up.
    """
    def handler(signum, frame):
        raise _Terminated()

    try:
        previous = signal.signal(signal.SIGTERM, handler)
    except ValueError:                 # not the main thread -- leave the default
        yield
        return
    try:
        yield
    finally:
        signal.signal(signal.SIGTERM, previous)


class Progress:
    """Live bucket occupancy: one row per target, a cell per tree, attempts below.

    Two modes, chosen from the stream, because the per-attempt line this replaces was
    unreadable on a terminal (hundreds of `no bucket` lines to find five keeps) and a
    repainting bar is worse than useless in a cluster log:

    * a TTY  -> redraw the block in place.
    * a file -> one line whenever a tree is kept, which is the only event worth a line,
                plus a heartbeat so a long dry spell still shows the run is alive.
    """

    HEARTBEAT = 250            # attempts between log lines when nothing is kept

    def __init__(self, targets, caps, counts, tty: Optional[bool] = None):
        self.targets = [int(t) for t in targets]
        self.caps = dict(caps)
        self.counts = dict(counts)
        self.attempts = 0
        self.kept = 0
        self.t0 = time.time()
        # getattr, not sys.stdout.isatty(): a run tees stdout into run.log through a
        # wrapper with no isatty, and treating that as "not a terminal" is also the right
        # answer -- repainting would write escape codes into the log file.
        self.tty = (getattr(sys.stdout, "isatty", lambda: False)()
                    if tty is None else tty)
        self._drawn = 0
        self._last_beat = 0

    def _bar(self, target: int) -> str:
        have, cap = self.counts.get(target, 0), self.caps[target]
        return "\u2588" * min(have, cap) + "\u00b7" * max(cap - have, 0)

    def _lines(self) -> List[str]:
        out = []
        for t in self.targets:
            have, cap = self.counts.get(t, 0), self.caps[t]
            mark = " " if have < cap else "\u2713"
            out.append(f"  eta{t:<3} {self._bar(t):<14} {have:>3}/{cap:<3} {mark}")
        rate = self.kept / self.attempts if self.attempts else 0.0
        out.append(f"  {self.attempts} attempts, {self.kept} kept "
                   f"({100 * rate:.0f}%), {time.time() - self.t0:.0f}s")
        return out

    def _redraw(self) -> None:
        if self._drawn:
            sys.stdout.write(f"\033[{self._drawn}A")
        for line in self._lines():
            sys.stdout.write("\033[2K" + line + "\n")
        sys.stdout.flush()
        self._drawn = len(self._lines())

    def attempt(self) -> None:
        self.attempts += 1
        if self.tty:
            self._redraw()
        elif self.attempts - self._last_beat >= self.HEARTBEAT:
            self._last_beat = self.attempts
            print(f"  {self.attempts} attempts, {self.kept} kept, "
                  f"have {[self.counts.get(t, 0) for t in self.targets]}", flush=True)

    def keep(self, target: int, stem: str) -> None:
        self.counts[target] = self.counts.get(target, 0) + 1
        self.kept += 1
        if self.tty:
            self._redraw()
        else:
            print(f"  kept {stem} -- "
                  f"{[self.counts.get(t, 0) for t in self.targets]} of "
                  f"{[self.caps[t] for t in self.targets]} "
                  f"after {self.attempts} attempts", flush=True)

    def start(self) -> None:
        if self.tty:
            self._redraw()

    def note(self, msg: str) -> None:
        """A line that must survive the repaint (an error, the closing summary)."""
        if self.tty and self._drawn:
            sys.stdout.write(f"\033[{self._drawn}A\033[J")
            self._drawn = 0
        print(msg, flush=True)


def stem_for(target: int, index: int) -> str:
    """``eta05_003`` -- the bucket is in the name, so ``ls`` answers which is which."""
    return f"eta{int(target):02d}_{int(index):03d}"


def needs_distance(operators: Sequence[str]) -> bool:
    """Whether any selected operator reads ``D``. Only the griffing arm (``B``) does."""
    from .operators import OPERATORS, resolve
    return any(OPERATORS[k].kind == "griffing" for k in resolve(operators))


# ---------------------------------------------------------------------------
# one attempt
# ---------------------------------------------------------------------------

def attempt(spec: SailfishSpec, seed: int, staging_root: Path,
            operators: Sequence[str] = DEFAULT_OPERATORS,
            binary: Optional[str] = None,
            min_split: int = 5, num_gaps: int = 10,
            rule_policy: str = "best") -> Dict:
    """Generate one tree + alignment and screen it. Never writes to the dataset.

    Returns the screen row with ``_seed`` and ``_staging`` attached, or a row carrying
    ``error`` when generation failed -- one bad seed must not end the run.
    """
    from .benchmark_loaders import RealLoader
    from .operator_screen import screen_one

    stage = Path(staging_root) / f"pid{os.getpid()}_seed{int(seed)}"
    stem = "attempt"
    try:
        tree = make_tree(spec, seed, stem, stage / "newick", binary)
        make_msa(spec, tree, seed, stage / "fasta")
        loader = RealLoader(stage / "fasta", stage / "newick",
                            need_distance=needs_distance(operators))
        row = screen_one(stem, loader, operators=operators, min_split=min_split,
                         rule_policy=rule_policy, num_gaps=num_gaps)
    except Exception as exc:                       # noqa: BLE001 -- the seed is skipped
        shutil.rmtree(stage, ignore_errors=True)
        return {"_seed": int(seed), "error": f"{type(exc).__name__}: {exc}"}

    if row is None:
        shutil.rmtree(stage, ignore_errors=True)
        return {"_seed": int(seed), "error": "loader returned nothing"}
    row["_seed"] = int(seed)
    row["_staging"] = str(stage)
    return row


def _init(spec, staging_root, operators, binary, min_split, num_gaps, policy) -> None:
    _CTX.update(spec=spec, staging_root=staging_root, operators=operators,
                binary=binary, min_split=min_split, num_gaps=num_gaps, policy=policy)


def _worker(seed: int) -> Dict:
    return attempt(_CTX["spec"], seed, _CTX["staging_root"], _CTX["operators"],
                   _CTX["binary"], _CTX["min_split"], _CTX["num_gaps"], _CTX["policy"])


# ---------------------------------------------------------------------------
# the loop
# ---------------------------------------------------------------------------

def verdict(row: Dict, counts: Dict[int, int], caps: Dict[int, int],
            cap_ops: Sequence[str], targets: Sequence[int],
            tol: float) -> Tuple[Optional[int], str]:
    """``(target, reason)`` -- the accept/reject decision for one screened tree.

    Bucket first, validity second, as asked. The order costs nothing either way: the
    O(m) bipartition lookup is already paid inside ``screen_one`` next to the O(m^3)
    eigensolve, so it only decides which word the log line carries.
    """
    if "error" in row:
        return None, "error"
    target = bucket_of(row, cap_ops, targets, tol)
    if target is None:
        return None, "no bucket"
    if counts.get(target, 0) >= caps.get(target, 0):
        return None, f"eta{target:02d} full"
    arm = cap_ops[0] if cap_ops else "S"
    if not row.get(f"valid_{arm}", False):
        return None, f"eta{target:02d} invalid split"
    return target, f"eta{target:02d}"


def existing_state(dataset: Path, source: str, cap_ops: Sequence[str],
                   targets: Sequence[int], tol: float) -> Tuple[Dict[int, int], List[Dict]]:
    """``(bucket counts, cached rows)`` for a dataset that is already part-filled."""
    path = screen_cache_path(source)
    if not path.exists():
        return {int(t): 0 for t in targets}, []
    rows = [dict(r) for r in np.load(path, allow_pickle=True)["rows"]]
    counts = {int(t): 0 for t in targets}
    on_disk = {p.stem for p in (dataset / "fasta").glob("*.fasta")}
    for r in rows:
        if r.get("tree") not in on_disk:
            continue
        t = bucket_of(r, cap_ops, targets, tol)
        if t is not None:
            counts[int(t)] = counts.get(int(t), 0) + 1
    return counts, rows


def fill_buckets(spec: SailfishSpec, out_root: Path, *,
                 name: Optional[str] = None,
                 targets: Sequence[int] = ETA_TARGETS,
                 per_bucket: Sequence[int] = PER_BUCKET,
                 tol: float = ETA_TOL,
                 max_attempts: int = 2000,
                 workers: int = 1,
                 operators: Sequence[str] = DEFAULT_OPERATORS,
                 gate: str = "valid_L",
                 seed_start: Optional[int] = None,
                 log: Callable[[str], None] = _log) -> Dict:
    """Generate until every bucket is at its cap. Resumable; returns run statistics.

    Written to survive a cluster: every accept flushes the screen cache *and* the
    manifest, and SIGTERM -- what SLURM sends at the time limit -- is caught and unwound
    like Ctrl+C. The manifest is the only record of which seeds were already tried, so
    leaving it to the end means a killed job resumes by re-drawing every tree it had
    already rejected.
    """
    targets = [int(t) for t in targets]
    caps = {int(t): int(c) for t, c in zip(targets, per_bucket)}
    cap_ops = gate_operators(gate)

    dataset = Path(out_root) / (name or spec.name)
    source = dataset.name
    staging_root = dataset / "_staging"
    binary = find_iqtree(spec)

    # A SIGKILLed predecessor leaves staging behind, and `dataset_scan.scan_roots`
    # discovers directories *named* fasta/newick -- a stale `_staging/pid*/fasta` would
    # be offered as a dataset of its own. Clear it before anything looks.
    shutil.rmtree(staging_root, ignore_errors=True)

    counts, rows = existing_state(dataset, source, cap_ops, targets, tol)
    prior = _prior_fill(dataset)
    seed = int(seed_start if seed_start is not None else _next_seed(dataset, spec))
    bar = Progress(targets, caps, counts)

    stats = {"attempts": 0, "accepted": 0, "no_bucket": 0, "invalid": 0,
             "full": 0, "errors": 0, "etas": [], "stopped": "exhausted"}
    t0 = time.time()

    def full() -> bool:
        return all(counts.get(t, 0) >= caps[t] for t in targets)

    if full():
        bar.note("every bucket is already at its cap -- nothing to do")
        return _finish(stats, counts, caps, targets, seed - 1, t0, staging_root, prior)

    def checkpoint() -> None:
        """Make the run's state durable. Cheap, so it runs on every accept."""
        _save_rows(source, rows)
        _write_fill(dataset, spec, binary,
                    _finish(dict(stats), counts, caps, targets, seed, t0,
                            staging_root, prior, sweep=False),
                    targets, [caps[t] for t in targets], tol)

    bar.start()
    seeds = iter(range(seed, seed + max_attempts))
    with _terminable(log):
        try:
            for row in _attempts(seeds, spec, staging_root, operators, binary, workers,
                                 spec_policy="best"):
                stats["attempts"] += 1
                bar.attempt()
                seed = max(seed, row.get("_seed", seed))
                target, reason = verdict(row, counts, caps, cap_ops, targets, tol)
                eta = row.get(f"eta_{cap_ops[0] if cap_ops else 'S'}")
                if eta is not None:
                    stats["etas"].append(float(eta))

                if target is None:
                    _tally(stats, reason)
                    shutil.rmtree(row.get("_staging", ""), ignore_errors=True)
                    if reason == "error":
                        bar.note(f"  seed={row['_seed']} FAILED: {row['error']}")
                        # Every attempt erroring is a broken environment, not bad luck:
                        # `attempt` swallows its own exceptions so one bad seed cannot
                        # end a run, which without this would spend all --max-attempts
                        # failing identically and report an empty dataset.
                        if stats["errors"] >= _ERROR_ABORT and not stats["accepted"]:
                            stats["stopped"] = "errors"
                            bar.note(f"abandoning: the first {stats['errors']} "
                                     f"attempts all failed the same way. Fix the error "
                                     f"above -- nothing here is retryable.")
                            break
                    continue

                counts[target] = counts.get(target, 0) + 1
                stem = stem_for(target, counts[target])
                _accept(row, Path(row["_staging"]), dataset, stem)
                rows = [r for r in rows if r.get("tree") != stem] + [row]
                stats["accepted"] += 1
                checkpoint()
                bar.keep(target, stem)
                if full():
                    stats["stopped"] = "complete"
                    break
        except (KeyboardInterrupt, _Terminated) as exc:
            how = "SIGTERM" if isinstance(exc, _Terminated) else "interrupt"
            stats["stopped"] = "signal"
            bar.note(f"\n{how} after {stats['attempts']} attempts -- {stats['accepted']} "
                f"accepted this run are on disk, in the screen cache and in the "
                f"manifest. Re-run the same command to continue from seed {seed + 1}.")

    out = _finish(stats, counts, caps, targets, seed, t0, staging_root, prior)
    _write_fill(dataset, spec, binary, out, targets,
                [caps[t] for t in targets], tol)
    return out


def _attempts(seeds, spec, staging_root, operators, binary, workers, spec_policy):
    """Yield screened rows, serially or across a process pool."""
    if workers <= 1:
        for s in seeds:
            yield attempt(spec, s, staging_root, operators, binary,
                          rule_policy=spec_policy)
        return

    import multiprocessing as mp

    # macOS spawns, so children re-import numpy and pick these up; setting them here
    # rather than at module import leaves the serial path at full BLAS width.
    for var in _BLAS_VARS:
        os.environ.setdefault(var, "1")
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=int(workers), initializer=_init,
                  initargs=(spec, staging_root, operators, binary, 5, 10,
                            spec_policy)) as pool:
        for row in pool.imap_unordered(_worker, seeds, chunksize=1):
            yield row


def _accept(row: Dict, stage: Path, dataset: Path, stem: str) -> None:
    """Move one staged attempt into the dataset under its final name."""
    for sub, ext in (("fasta", ".fasta"), ("newick", ".nwk")):
        (dataset / sub).mkdir(parents=True, exist_ok=True)
        src = next(iter((stage / sub).glob(f"*{ext}")))
        shutil.move(str(src), str(dataset / sub / f"{stem}{ext}"))
    shutil.rmtree(stage, ignore_errors=True)
    row["tree"] = stem
    row.pop("_staging", None)


def _save_rows(source: str, rows: List[Dict]) -> None:
    """Rewrite the screen cache so an accepted tree is never screened again."""
    path = screen_cache_path(source)
    path.parent.mkdir(parents=True, exist_ok=True)
    clean = [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows]
    np.savez(path, rows=np.array(clean, dtype=object))


def _next_seed(dataset: Path, spec: SailfishSpec) -> int:
    """Where to resume seeding.

    From the manifest, not from the screen cache: ``_save_rows`` strips the private
    fields so a cached row stays a plain screen row, which means the cache cannot say
    which seeds were already tried. Reading it from there would restart at the base seed
    and re-draw every tree the previous run had already rejected.
    """
    import json
    path = dataset / "manifest.json"
    if path.exists():
        try:
            last = json.loads(path.read_text()).get("fill", {}).get("last_seed")
        except (OSError, ValueError):
            last = None
        if last is not None:
            return int(last) + 1
    return int(spec.tree_seed_base)


def _tally(stats: Dict, reason: str) -> None:
    if reason == "error":
        stats["errors"] += 1
    elif reason == "no bucket":
        stats["no_bucket"] += 1
    elif reason.endswith("full"):
        stats["full"] += 1
    else:
        stats["invalid"] += 1


def _fmt(eta) -> str:
    return "n/a" if eta is None else f"{float(eta):.2f}"


def _prior_fill(dataset: Path) -> Dict:
    """The ``fill`` block of an earlier run, so a resume adds to it rather than
    replacing it. A no-op resume otherwise rewrites the manifest with zero attempts
    and a 0% acceptance rate, erasing the only record of what the fill cost."""
    import json
    path = dataset / "manifest.json"
    if not path.exists():
        return {}
    try:
        return dict(json.loads(path.read_text()).get("fill", {}))
    except (OSError, ValueError):
        return {}


def _write_fill(dataset: Path, spec, binary: str, stats: Dict,
                targets: Sequence[int], caps: Sequence[int], tol: float) -> None:
    """Persist the ``fill`` block. Called on every accept, not only at the end."""
    from ..models.sailfish_sim import write_manifest
    write_manifest(
        dataset, spec, binary, elapsed=stats.get("elapsed_sec", 0.0),
        fill={"mode": "fill_buckets", "eta_targets": list(targets),
              "per_bucket": list(caps), "eta_tol": float(tol),
              "attempts": stats["attempts_total"],
              "accepted": stats["accepted_total"],
              "acceptance": stats["acceptance_total"],
              "last_seed": stats["last_seed"],
              "counts": stats["counts"], "complete": stats["complete"]})


def _finish(stats, counts, caps, targets, seed, t0, staging_root,
            prior: Optional[Dict] = None, sweep: bool = True) -> Dict:
    prior = prior or {}
    if sweep:                          # a mid-run checkpoint must not delete live staging
        shutil.rmtree(staging_root, ignore_errors=True)
    stats["elapsed_sec"] = round(time.time() - t0, 1)
    stats["last_seed"] = int(seed)
    stats["counts"] = {int(t): int(counts.get(t, 0)) for t in targets}
    stats["caps"] = {int(t): int(caps[t]) for t in targets}
    stats["complete"] = all(counts.get(t, 0) >= caps[t] for t in targets)
    a = stats["attempts"]
    stats["acceptance"] = round(stats["accepted"] / a, 4) if a else 0.0
    # cumulative across every run that built this dataset
    stats["attempts_total"] = int(prior.get("attempts", 0)) + a
    stats["accepted_total"] = int(prior.get("accepted", 0)) + stats["accepted"]
    tot = stats["attempts_total"]
    stats["acceptance_total"] = (round(stats["accepted_total"] / tot, 4)
                                 if tot else 0.0)
    if not a and prior.get("last_seed") is not None:
        stats["last_seed"] = int(prior["last_seed"])   # a no-op must not rewind
    return stats
