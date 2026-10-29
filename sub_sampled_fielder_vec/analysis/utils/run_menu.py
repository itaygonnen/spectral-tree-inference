"""Interactive real-data menu for ``scripts/interactive_run.py``.

The launcher's other branch configures a *simulated* sweep (tree model, mutation rate,
sequence length). Real datasets have none of those knobs -- the alignment and the true
tree are given -- so this asks a different, shorter set of questions:

    which datasets  any <name>/{fasta,newick} under data/tree_sets/, several at a time
    which stage    screening (eta + validity per operator) or recovery sweep (vs p)
    parameters     asked ONCE and applied to every chosen dataset, so sizes stay
                   comparable -- the whole point of running m=1000 beside m=6000

This module is only the questions and the status panel. Selecting trees, planning,
running and exporting all live in :mod:`src.runners.experiment_run`, which
``scripts/run_sweep.py`` calls too -- so a menu session on a login node and a nohup'd
batch run produce the same run directory from the same cache. That runner knows nothing
about FASTA: this module's real job is turning the answers into a ``RunSpec`` and a list
of ``Source`` objects.

Lives beside the dataset helpers in ``analysis/utils/`` -- ``src/`` must never import
from ``analysis/``, so the dependency runs launcher -> analysis -> src, never back.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from dataclasses import replace

_ROOT = Path(__file__).resolve().parents[2]      # sub_sampled_fielder_vec
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.interactive_ui import (                       # noqa: E402
    confirm, get_input, get_menu_choice, get_multi_choice, print_divider, print_error,
    print_header, print_option, print_success, print_warning)

from .real_datasets import (DATASETS_ROOT, describe_search,          # noqa: E402
                            list_datasets, sweep_cache_dir)
from src.runners.eta_buckets import (ETA_TARGETS, ETA_TOL, PER_BUCKET,  # noqa: E402
                                     UNBUCKETED, label,
                                     panel_lines)
from src.runners.experiment_run import (DEFAULT_GATE, EARLY_STOP,   # noqa: E402
                                        GATES, MAX_ETA,
                                        P_MIN, P_POINTS, REPS, RunSpec, as_command,
                                        default_gate, execute, gate_label,
                                        gate_operators, plan, select_buckets, verdicts,
                                        already_swept)
from src.runners.operators import ALL_OPERATORS, ARM_OF, OPERATORS  # noqa: E402


def _cmd_width() -> int:
    """Width to wrap a printed command at, leaving room for the two-space indent."""
    import shutil
    return max(40, min(shutil.get_terminal_size(fallback=(80, 24)).columns, 100) - 8)


def _say(text: str, indent: str = "") -> None:
    """Print prose wrapped to the terminal.

    Long lines are not a cosmetic problem here: the line before a prompt used to be 160
    characters, so on an 80-column terminal it wrapped twice and the answer landed
    mid-sentence -- which reads as the program having crashed.
    """
    import shutil
    import textwrap

    width = max(40, min(shutil.get_terminal_size(fallback=(80, 24)).columns, 100))
    for line in textwrap.wrap(text, width=width, initial_indent=indent,
                              subsequent_indent=indent + "  ") or [indent]:
        print(line)


_SOURCES_BY_NAME: Dict[str, object] = {}     # set by _run so the panel can count
                                             # trees that are not screened yet


def _ask_config(names, is_screen: bool, verdicts_by: Dict[str, dict],
                pooled_as: Optional[tuple] = None) -> RunSpec:
    """One parameter set for every source: show the defaults, edit them only on request.

    ``pooled_as`` is ``(eta_targets, per_bucket)`` when the trees were drawn from the eta
    pool to those targets. The sweep then selects on the same terms they were built on,
    instead of overriding them with the paper's caps -- a run pooled at 10 per cell was
    otherwise cut to 10/10/3/3, silently discarding 14 of the 40 trees it had just
    simulated. Bucketing itself is the same code either way (``select_ids``); what the
    pool changes is only what the buckets are expected to hold.
    """
    spec = RunSpec(stage="screen" if is_screen else "sweep",
                   max_trees=0,                   # 0 = every tree in each dataset
                   workers=max(1, min(8, (os.cpu_count() or 4) // 2)),
                   gate=default_gate(names, verdicts_by),
                   max_eta=MAX_ETA, p_min=P_MIN, p_points=P_POINTS, reps=REPS,
                   early_stop=EARLY_STOP,
                   eta_targets=(pooled_as[0] if pooled_as else ETA_TARGETS),
                   per_bucket=(pooled_as[1] if pooled_as else PER_BUCKET),
                   extra_metrics=True, prefix="", display_mode="progress")

    print_divider()
    print("configuration (every source):")
    if is_screen:
        print(f"  trees      all, {spec.workers} workers")
        print(f"  operators  {', '.join(spec.operators)}")
    else:
        # Only what varies and changes the run. The bucket settings used to be echoed
        # here as well, and they are the first line of the panel immediately below.
        _say(f"gate       [{gate_label(spec.gate)}], eta <= {spec.max_eta:g}",
             indent="  ")
        stop = (f"early stop {spec.early_stop}" if spec.early_stop else "no early stop")
        print(f"  p-grid     {spec.p_points} points, {spec.p_min:g} .. 1.0, "
              f"{spec.reps} reps, {stop}")
        print(f"  operators  {', '.join(spec.operators)}")

    print()
    if not confirm("Edit this configuration?", default=False):
        return spec

    cap = get_input("Max trees per source (blank = all)", default="")
    spec.max_trees = int(cap) if cap and cap.strip() else 0
    if is_screen:
        # screening has no free parameters, so it has one canonical output per dataset
        spec.workers = int(get_input("Workers", default=str(spec.workers)))
        return spec

    # Spell the counts out. A bare "2, 91" reads as a range, and even "2/66" hides that
    # the denominator a gate can act on is the SCREENED trees, not the dataset.
    # No counts on the options: the Screening status table above already gives coverage
    # and verdicts, and the plan printed after this shows what the choice actually selects.
    keys = list(GATES)
    labels = [gate_label(k) for k in keys]
    chosen = get_menu_choice("Reference partition must be a real tree edge under:",
                             labels, default_index=keys.index(spec.gate))
    spec.gate = keys[labels.index(chosen)]
    spec.max_eta = float(get_input(
        "Drop trees whose reference split is more lopsided than eta (0 = keep all)",
        default=str(spec.max_eta)))
    _ask_buckets(spec, names, verdicts_by)
    spec.p_min = float(get_input("Smallest sub-sampling rate p (largest is always 1.0)",
                                 default=str(spec.p_min)))
    spec.p_points = int(get_input(
        f"p-grid points (log-spaced, {spec.p_min:g}..1)", default=str(spec.p_points)))
    spec.reps = int(get_input("Bootstrap reps per p", default=str(spec.reps)))
    spec.early_stop = int(get_input(
        "Settled points before an end of the grid is extrapolated (0 = measure every p)",
        default=str(spec.early_stop)))
    ops = get_input("Operators (comma-separated: " + ",".join(ALL_OPERATORS) + ")",
                    default=",".join(spec.operators))
    spec.operators = tuple(o.strip() for o in ops.split(",") if o.strip())
    spec.extra_metrics = confirm(
        "Also record sigma2, ||sub-full||_2 and the numerical ranks at every p?",
        default=True)
    name = get_input("Name for this run (blank = timestamp only)", default="")
    spec.prefix = (name or "").strip()
    spec.display_mode = get_menu_choice(
        "Output:", ["progress - bars, and every line in the run's experiment.log",
                    "debug - every line on the terminal too, no bars"],
        default_index=0).split(" ")[0]
    return spec


def _ask_buckets(spec: RunSpec, names, verdicts_by: Dict[str, dict]) -> None:
    """Show what each eta bucket holds, then take the caps. Mutates ``spec`` in place.

    This is where a run stops paying for trees it does not need, so it is a dialogue
    rather than a flag: on real data the paper's windows at +-1 keep a minority of the
    valid trees and some buckets come up empty, and the only way to find out is to look.
    The panel is redrawn after the tolerance so a choice can be made against the numbers
    it produces rather than remembered from before it.

    Sets ``spec.bucket``, ``spec.eta_tol``, ``spec.eta_targets`` and ``spec.per_bucket``.
    """
    spec.bucket = confirm(
        "Restrict the sweep to the paper's eta buckets?  (n sweeps every gated tree)",
        default=spec.bucket)
    if not spec.bucket:
        return
    raw = get_input("Eta buckets - comma-separated targets",
                    default=",".join(str(t) for t in spec.eta_targets))
    spec.eta_targets = tuple(int(x) for x in raw.split(",") if x.strip())
    if len(spec.per_bucket) != len(spec.eta_targets):
        # keep the parallel lists parallel: reuse the last cap for any new target
        last = spec.per_bucket[-1] if spec.per_bucket else 3
        spec.per_bucket = tuple(
            (spec.per_bucket[i] if i < len(spec.per_bucket) else last)
            for i in range(len(spec.eta_targets)))

    while True:
        tol = get_input("Bucket half-width (a tree joins the nearest target within it)",
                        default=f"{spec.eta_tol:g}")
        new_tol = float(tol)
        if new_tol == spec.eta_tol:
            break
        spec.eta_tol = new_tol          # redraw once at the new tolerance, then move on
        _bucket_panels(spec, names, verdicts_by, legend=False)
        if confirm("Keep this half-width?", default=True):
            break

    raw = get_input("Trees per bucket - one count per target, in order",
                    default=",".join(str(c) for c in spec.per_bucket))
    caps = [int(x) for x in raw.split(",") if x.strip()]
    if len(caps) == len(spec.eta_targets):
        spec.per_bucket = tuple(caps)
    else:
        print_warning(f"need {len(spec.eta_targets)} counts, got {len(caps)} -- keeping "
                      + ",".join(str(c) for c in spec.per_bucket))
    _bucket_panels(spec, names, verdicts_by, legend=False)


def _bucket_panels(spec: RunSpec, names, verdicts_by: Dict[str, dict],
                   legend: bool = True, is_screen: bool = False) -> None:
    """One occupancy panel per source, under the spec as it stands.

    Shown on every run. A tree's bucket is its MEASURED eta, so on a screening run the
    trees not yet screened have no bucket yet -- they are counted on their own line
    rather than guessed at, and the table then reads as progress towards the buckets
    instead of a selection.

    ``legend=False`` for a redraw: the dialogue shows each panel up to three times while
    the tolerance and the caps are settled, and three copies of the same footer pushes
    the numbers off the screen.
    """
    if len(names) > 1:
        _bucket_matrix(spec, names, verdicts_by, is_screen=is_screen)
        return

    for i, name in enumerate(names):
        rows = verdicts_by.get(name) or {}
        n_all = len(_ids_of(name, names, verdicts_by))
        kept, buckets, _short = (select_buckets(sorted(rows), rows, spec)
                                 if rows else ([], {}, {}))
        # the screen row carries the taxon count, which the sweep-cache meta is keyed on
        m = next((r.get("m", 0) for r in rows.values() if r.get("m")), 0)
        done = already_swept(name, m, kept, spec) if kept else set()
        todo = len(kept) - len(done)
        head = (f"{name}  ->  {len(rows)} of {n_all} screened"
                if is_screen else
                f"{name}  ->  {todo} tree(s) to sweep"
                + (f"  ({len(done)} already done)" if done else ""))
        body = panel_lines(buckets, spec.eta_targets, spec.per_bucket, spec.eta_tol,
                           legend=(legend and i == 0), done=done)
        if is_screen and n_all > len(rows):
            # inserted after the data rows, before TOTAL/legend, so the column that
            # cannot be known yet is visible rather than silently missing
            at = next(j for j, ln in enumerate(body) if ln.startswith("TOTAL"))
            body.insert(at, f"{'unscreened':<11}{n_all - len(rows):>7}{'-':>7}{'-':>7}"
                            f"  bucket unknown until screened")
        _framed([head] + body, rule_after=0)



def _bucket_matrix(spec: RunSpec, names, verdicts_by: Dict[str, dict],
                   is_screen: bool = False) -> None:
    """One ROW per dataset, buckets as columns -- the size-series form of the panel.

    The per-dataset bar panel is right for one dataset and unreadable for several: four
    sizes drew four frames of thirteen lines, each repeating the same window/cap header
    and the same legend, for ~52 lines of mostly duplicate text. The targets and caps are
    constant down every column, so they are stated once in the header, exactly as the
    single panel already does with its window and cap.
    """
    targets = list(spec.eta_targets)
    caps = list(spec.per_bucket)
    shown = _short_names(list(names))
    width = max(max((len(v) for v in shown.values()), default=7), 7)
    head = (f"  {'dataset':<{width}}{'m':>7}"
            + "".join(f"{label(t):>9}" for t in targets)
            + f"{'none':>7}{'to run':>9}")
    print(f"  nearest target within +-{spec.eta_tol:g}, at most "
          + "/".join(str(c) for c in caps) + " per bucket")
    print(head)
    total = 0
    for name in names:
        rows = verdicts_by.get(name) or {}
        n_all = len(_ids_of(name, names, verdicts_by))
        kept, buckets, _short = (select_buckets(sorted(rows), rows, spec)
                                 if rows else ([], {}, {}))
        m = next((r.get("m", 0) for r in rows.values() if r.get("m")), 0)
        done = already_swept(name, m, kept, spec) if kept else set()
        todo = len(rows) if is_screen else len(kept) - len(done)
        total += todo
        cells = "".join(f"{len(buckets.get(t, [])):>4}/{c:<4}"
                        for t, c in zip(targets, caps))
        print(f"  {shown[name]:<{width}}{m:>7}{cells}"
              f"{len(buckets.get(UNBUCKETED, [])):>7}{todo:>9}")
        if is_screen and n_all > len(rows):
            print(f"  {'':<{width}}{'':>7}" + " " * (9 * len(targets))
                  + f"{'':>7}{n_all - len(rows):>9}  unscreened")
    print(f"  {'TOTAL':<{width}}{'':>7}" + " " * (9 * len(targets))
          + f"{'':>7}{total:>9}  <- this run")


def _short_names(names: List[str]) -> Dict[str, str]:
    """Names with the part every selection shares stripped off both ends.

    A size series differs in one field, so printing the full name four times spends the
    column on text that is identical down it -- and then truncates the only part that is
    not. What is common is stated once in the heading above instead.
    """
    if len(names) < 2:
        return {n: n for n in names}
    # Snap to "_" boundaries: a raw common prefix cuts mid-token, turning n200/n250 into
    # "20"/"25" -- the digits that differ, with the field name they belong to removed.
    pre = os.path.commonprefix(names)
    pre = pre[:pre.rfind("_") + 1] if "_" in pre else ""
    rev = os.path.commonprefix([n[::-1] for n in names])[::-1]
    suf = rev[rev.find("_"):] if "_" in rev else ""
    out = {}
    for n in names:
        mid = n[len(pre):len(n) - len(suf)] if len(pre) + len(suf) < len(n) else n
        out[n] = mid or n
    return out

def _ids_of(name, names, verdicts_by) -> list:
    """Every tree id the source offers, screened or not. Falls back to the screened ones
    when the caller did not carry the source objects (the panel is also drawn from
    ``--dry-run``, which has them, and from the menu, which has the names)."""
    src = _SOURCES_BY_NAME.get(name)
    return list(src.ids) if src is not None else list(verdicts_by.get(name) or {})


def _framed(lines: List[str], rule_after: int = -1) -> None:
    """Draw ``lines`` inside a box, so a status panel is not mistaken for a menu.

    A box wider than the terminal is worse than no box: every line wraps and the borders
    land mid-sentence, which is what an 80-column terminal did to this panel once it
    carried a column pair per operator. When it will not fit, rule lines are used
    instead -- the panel still reads as a panel, and nothing is mangled.
    """
    import shutil

    width = max(len(x) for x in lines)
    cols = shutil.get_terminal_size(fallback=(80, 24)).columns
    if width + 4 > cols:
        rule = "─" * min(cols, width)
        print(rule)
        for i, line in enumerate(lines):
            print(line)
            if i == rule_after:
                print(rule)
        print(rule)
        print()
        return
    print("┌" + "─" * (width + 2) + "┐")
    for i, line in enumerate(lines):
        print(f"│ {line.ljust(width)} │")
        if i == rule_after:
            print("├" + "─" * (width + 2) + "┤")
    print("└" + "─" * (width + 2) + "┘")
    print()


def _print_status(sources, verdicts_by: Dict[str, dict]) -> None:
    """One line per source: how much is screened, how much is usable, how much is swept.

    This used to be a framed table of per-operator columns plus eight lines of legend --
    19 lines for two sources, most of it restating what the bucket panel below it says
    better. What survives is what the table alone could tell you: coverage, whether the
    trees are usable at all, and how lopsided they are.
    """
    print_header("Screening status")
    for c in sources:
        v, ids = verdicts_by[c.name], c.ids
        done = [t for t in ids if t in v]
        swept = (len(list(sweep_cache_dir(c.name).glob("*.npz")))
                 if sweep_cache_dir(c.name).is_dir() else 0)
        if not done:
            _say(f"{c.name}: {len(ids)} trees, none screened yet", indent="  ")
            continue
        key = f"eta_{gate_operators(DEFAULT_GATE)[0]}"
        etas = [v[t][key] for t in done if key in v[t]]
        n_ok = sum(1 for t in done if GATES[DEFAULT_GATE][1](v[t]))
        med = f", median eta {float(np.median(etas)):.1f}" if etas else ""
        _say(f"{c.name}: {len(done)}/{len(ids)} screened, {n_ok} usable"
             f"{med}, {swept} swept", indent="  ")
    print()


def _run(sources, *, batch_hint: str, pooled_as: Optional[tuple] = None,
         reconfigure=None) -> None:
    """Status, stage, one configuration, the plan, then the run. Source-agnostic.

    ``reconfigure`` is offered as a third choice on the stage menu and returns a fresh
    ``(sources, pooled_as)`` or None to quit. It is a menu entry rather than a prompt so
    that NOTHING is asked before the status is on screen: changing the settings is a
    thing you do after seeing that the defaults do not suit, not before.
    """
    while True:
        verdicts_by = {c.name: verdicts(c.name) for c in sources}
        _SOURCES_BY_NAME.clear()
        _SOURCES_BY_NAME.update({c.name: c for c in sources})
        _print_status(sources, verdicts_by)

        # The bucket table FIRST -- straight after the source is picked, before any
        # question at all. It is the status this launcher exists to show, and burying it
        # behind prompts meant answering several to find out whether there was anything
        # worth running. Occupancy does not depend on the stage, so there is nothing to
        # wait for.
        _preview = RunSpec(stage="sweep", gate=default_gate(
            [c.name for c in sources], verdicts_by),
            **(dict(eta_targets=pooled_as[0], per_bucket=pooled_as[1])
               if pooled_as else {}))
        _bucket_panels(_preview, [c.name for c in sources], verdicts_by)

        options = ["screening - split the full matrix, check eta and validity",
                   "recovery sweep - sweep trees over p for NMI curve "
                   "(needs screening)"]
        if reconfigure is not None:
            options.append("change the data settings - model, size, pooling")
        stage = get_menu_choice("Run:", options, default_index=0)
        if reconfigure is not None and stage.startswith("change"):
            got = reconfigure()
            if not got or not got[0]:
                return
            sources, pooled_as = got
            continue
        break
    is_screen = stage.startswith("screening")

    spec = _ask_config([c.name for c in sources], is_screen, verdicts_by, pooled_as)
    if spec.max_trees:
        sources = [replace(c, ids=c.ids[:spec.max_trees]) for c in sources]
    rows = plan(spec, sources, verdicts_by)

    print_divider()
    if not is_screen:
        pv = spec.p_values()
        print(f"  p-grid     {len(pv)} log-spaced points, {pv[0]:.4g} .. {pv[-1]:.4g}")
    width = max([12] + [len(r.name) + 1 for r in rows])
    total = sum(r.hours for r in rows)
    for r in rows:
        print(f"  {r.name:<{width}} m={r.m:<5} {len(r.selected):>4} trees"
              f"  ~{r.hours:.1f} h")
        if not is_screen:
            # where the trees went: the two filters, in the order they are applied
            _say(f"{r.n_screened} screened -> {r.n_gate} "
                 f"[{gate_label(spec.gate)}] -> {len(r.selected)} with eta <= "
                 f"{spec.max_eta:g}", indent="      ")
        if r.warning:
            print_warning(f"{r.name}: {r.warning}")
    print()
    _say(f"total ~{total:.1f} h"
         + ("" if is_screen else " (the sweep runs one tree at a time; workers "
                                "parallelise the screen only)")
         + ". Every tree is cached on its own, so this is safe to interrupt and "
           "resume.")
    if not confirm("Run it here now?  (n prints the command instead)",
                   default=total <= 1.0):
        # "n" is not a cancellation, it is a request: the command carrying the answers
        # just given is the whole point of having answered them. Printing a warning here
        # made choosing "run it elsewhere" look like a failure.
        print()
        print_success("Not running here. This is the run you configured:")
        print()
        print("  " + as_command(spec, sources, width=_cmd_width(), indent="      "))
        print()
        print("  # detached, with a log:")
        print("  nohup " + as_command(spec, sources, width=_cmd_width(),
                                      indent="      ")
              + f" \\\n      > logs/{spec.stage}.log 2>&1 &")
        print()
        # A file beats a long command line for an sbatch script: it survives editing,
        # it is diffable, and it is the same shape every run already writes.
        where = get_input("Save this as a config file? (path, or Enter to skip)",
                          default="")
        if where:
            _save_run_config(where, spec, [c.name for c in sources])
        print()
        _say("Add --dry-run to either one to see the selection without running it. "
             "Nothing has been written.")
        return

    def _announce(k, n, row):
        print_divider()
        print_header(f"[{k}/{n}] {row.name}")

    execute(spec, sources, rows, on_source=_announce)
    print_success("done")


def run_real_data_menu() -> None:
    """Datasets on disk: pick the FASTA/newick directories, then run.

    Downloaded and simulated sets share one list because they are the same thing to
    everything downstream -- a directory of alignments beside their trees. They are
    *labelled* apart because they are not the same evidence, and an eta-filled set is
    labelled again because its eta distribution is conditioned by rejection sampling.

    A size series is the normal case once you can simulate one, so sizes of the same
    family are grouped under one letter and sorted by m. Each size stays its OWN source:
    the cost model and every per-row output are keyed on the taxon count, so merging two
    sizes into one source would make both meaningless.
    """
    datasets = list_datasets()
    if not datasets:
        print_error("No datasets found.")
        print(describe_search())
        return

    groups = _group_datasets(datasets)
    ordered = [c for _, members in groups for c in members]

    print_header("Datasets on disk")
    _say("Select with commas, a range ('5-8'), a family letter, or 'all'", indent="  ")
    print()
    letters, n = {}, 0
    for family, members in groups:
        if len(members) > 1:
            letter = chr(ord("a") + len(letters))
            letters[letter] = list(range(n + 1, n + 1 + len(members)))
            sizes = ", ".join(str(c.shape()[0]) for c in members)
            print(f"  [{letter}] {family}  —  {len(members)} sizes: m = {sizes}")
        for c in members:
            n += 1
            m, seq_len = c.shape()
            print_option(str(n),
                         f"{c.name}  ({len(c.ids())} trees, m={m}, L={seq_len})",
                         highlight=False)
            fill = c.bucket_fill()
            print(f"      {c.provenance()}{'   buckets ' + fill if fill else ''}")
        if len(members) > 1:
            print()

    print()
    picks = _pick_datasets(len(ordered), letters)
    if not picks:
        return
    chosen = [ordered[i - 1] for i in picks]
    _say(f"running {len(chosen)}: " + ", ".join(c.name for c in chosen), indent="  ")
    _run([c.source() for c in chosen], batch_hint="python scripts/run_sweep.py")


def _family_of(dataset) -> str:
    """The name with its taxon count removed -- what two sizes of one series share."""
    import re
    return re.sub(r"_n\d+_", "_", dataset.name) or dataset.name


def _group_datasets(datasets):
    """``[(family, [dataset sorted by m])]``, families in first-seen order.

    Sorting matters: plain alphabetical puts n1000 before n200, which reads as a bug the
    first time you generate a size series.
    """
    fams: Dict[str, list] = {}
    for c in datasets:
        fams.setdefault(_family_of(c), []).append(c)
    return [(fam, sorted(members, key=lambda d: d.shape()[0]))
            for fam, members in fams.items()]


def _pick_datasets(count: int, letters: Dict[str, list]):
    """Numbers, ranges, family letters and 'all' -> a list of 1-based indices."""
    while True:
        raw = (get_input("Dataset(s)") or "").strip().lower()
        if not raw:
            continue
        if raw in ("q", "quit", "back"):
            return []
        if raw == "all":
            return list(range(1, count + 1))
        picked, bad = [], []
        for tok in (t.strip() for t in raw.split(",") if t.strip()):
            if tok in letters:
                picked += letters[tok]
            elif "-" in tok and all(part.strip().isdigit()
                                    for part in tok.split("-", 1)):
                lo, hi = (int(part) for part in tok.split("-", 1))
                picked += list(range(min(lo, hi), max(lo, hi) + 1))
            elif tok.isdigit():
                picked.append(int(tok))
            else:
                bad.append(tok)
        picked = [i for i in dict.fromkeys(picked) if 1 <= i <= count]
        if bad:
            print_warning(f"not a number, range, family letter or 'all': "
                          f"{', '.join(bad)}")
        elif picked:
            return picked
        else:
            print_warning(f"nothing in range 1-{count}")


def run_generated_menu() -> None:
    """Simulated trees: ask the model, sizes and sequence length, then run.

    The same experiment as the real branch -- a population of trees, screened, gated and
    swept, with the curve a median over trees. The simulated branch used to run one tree
    per size and take its spread from bootstrap replicates, which is a different
    quantity and was not comparable with a real-data curve.
    """
    from src.runners.generated_source import sources_from_plan
    from src.utils.generated_prompts import (baseline_plan, describe_plan,
                                             prompt_generated_plan, resolve_pool_ids)

    def _build(plan_g):
        """A plan -> (sources, pooled_as), or (None, None) when it yields no trees."""
        _say(describe_plan(plan_g), indent="  ")
        if plan_g.pooled:
            # the eta pool is the only place a rejection-sampled tree exists; top it
            # up before anything asks for ids, because what comes back is what is on
            # disk
            plan_g.tree_ids = resolve_pool_ids(plan_g)
            if not plan_g.tree_ids:
                print_error("the eta pool is empty for these parameters and nothing "
                            "could be built -- nothing to run")
                return None, None
        srcs = sources_from_plan(plan_g)
        if not srcs:
            print_error("no trees to run")
            return None, None
        # selected on the terms they were pooled on, not the paper's caps
        pooled = ((tuple(plan_g.etas), (plan_g.per_cell,) * len(plan_g.etas))
                  if plan_g.pooled and plan_g.etas else None)
        return srcs, pooled

    print_header("Generated data")
    # Start from the model's baseline -- no questions -- so the status is the first
    # thing on screen. The six-step dialogue is the "change the data settings" entry on
    # the menu below.
    sources, pooled_as = _build(baseline_plan())
    if not sources:
        return
    _run(sources, batch_hint="", pooled_as=pooled_as,
         reconfigure=lambda: _build(prompt_generated_plan()))


# ---------------------------------------------------------------------------
# Simulate a dataset (IQ-TREE + Sailfish), filled to the eta buckets
# ---------------------------------------------------------------------------

def run_simulate_menu() -> bool:
    """Ask for a size, measure the hit rate, then fill the buckets. True if it wrote one.

    Calibration is a prompt rather than a flag because the acceptance rate is the number
    that decides the run and it is not knowable in advance: at m=1000 the measured eta
    median is 19.8, above every target, and only ~7% of trees are kept. Showing that
    before the fill is the difference between a 2-minute job and an abandoned overnight
    one. The CLI twin is ``scripts/generate_dataset.py --fill-buckets``.
    """
    from src.models.sailfish_sim import SailfishSpec, preflight
    from src.runners.sailfish_pool import fill_buckets

    print_header("Simulate a dataset")
    missing = preflight()
    if missing:
        print_error("this branch needs two things that are not installed:")
        for line in missing:
            print(f"    - {line}")
        return "back"

    # Start from the defaults and show them, rather than asking seven questions nobody
    # wants to answer twice -- the same shape as the generated branch, where the status
    # is the first thing on screen and editing is a menu entry.
    cfg = dict(sizes=[1000], seq_len=1000, targets=list(ETA_TARGETS),
               caps=list(PER_BUCKET), tol=ETA_TOL, workers=4, max_attempts=2000)
    while True:
        _show_sim_config(cfg)
        pick = get_menu_choice("Run:", ["start", "change these settings", "back"],
                               default_index=0)
        if pick.startswith("back"):
            return "back"
        if pick.startswith("start"):
            break
        cfg = _edit_sim_config(cfg)

    specs = [SailfishSpec(n_taxa=n, root_seq_size=cfg["seq_len"], pooled=True)
             for n in cfg["sizes"]]
    wrote = False
    for spec in specs:
        if len(specs) > 1:
            print_divider()
            _say(f"{spec.n_taxa} taxa", indent="  ")
        stats = fill_buckets(spec, DATASETS_ROOT, targets=cfg["targets"],
                             per_bucket=cfg["caps"], tol=cfg["tol"],
                             workers=cfg["workers"],
                             max_attempts=cfg["max_attempts"])
        got = [stats["counts"][int(t)] for t in cfg["targets"]]
        wrote = True
        print()
        if stats["complete"]:
            print_success(f"{spec.name}: every bucket full {got} "
                          f"({stats['accepted']} kept of {stats['attempts']} tried)")
        else:
            print_warning(f"{spec.name}: incomplete {got} of {cfg['caps']} -- re-run "
                          f"to continue from seed {stats['last_seed'] + 1}, or widen "
                          f"the half-width")
    return wrote


def _show_sim_config(cfg: Dict) -> None:
    """The whole configuration in four lines, so starting needs no reading."""
    sizes = ", ".join(str(n) for n in cfg["sizes"])
    buckets = ", ".join(f"{c} at eta{t}"
                        for t, c in zip(cfg["targets"], cfg["caps"]))
    print()
    _say(f"taxa      {sizes}    alignment L = {cfg['seq_len']}", indent="  ")
    _say(f"collect   {buckets}   (+-{cfg['tol']:g})", indent="  ")
    _say(f"          {sum(cfg['caps'])} trees per size, "
         f"{cfg['workers']} workers, up to {cfg['max_attempts']} attempts",
         indent="  ")
    for n in cfg["sizes"]:
        _say(f"-> sailfish_yule_n{n}_L{cfg['seq_len']}_g1x4_buckets", indent="  ")
    print()


def _edit_sim_config(cfg: Dict) -> Dict:
    """The full dialogue, reached only by asking for it."""
    cfg = dict(cfg)
    cfg["sizes"] = _ask_ints("Taxa per tree - comma-separated for several sizes",
                             ",".join(map(str, cfg["sizes"])))
    cfg["seq_len"] = _ask_ints("Alignment length L (root sequence size)",
                               str(cfg["seq_len"]))[0]
    cfg["targets"] = _ask_ints("Eta targets - the imbalance levels to collect at",
                               ",".join(map(str, cfg["targets"])))
    named = ", ".join(f"eta{t}" for t in cfg["targets"])
    default_caps = (",".join(map(str, cfg["caps"]))
                    if len(cfg["caps"]) == len(cfg["targets"]) else "10")
    while True:
        caps = _ask_ints(f"How many trees at each of {named}"
                         f"\n  (one number per target, or a single number for all)",
                         default_caps)
        if len(caps) == 1:
            caps = caps * len(cfg["targets"])
        if len(caps) == len(cfg["targets"]):
            cfg["caps"] = caps
            break
        print_warning(f"give {len(cfg['targets'])} numbers, one per target ({named}), "
                      f"or a single number for all - got {len(caps)}")
    cfg["tol"] = _ask_float("Bucket half-width - a tree counts as eta5 within +-",
                            str(cfg["tol"]))
    cfg["workers"] = _ask_ints("Workers", str(cfg["workers"]))[0]
    cfg["max_attempts"] = _ask_ints("Max attempts per size",
                                    str(cfg["max_attempts"]))[0]
    return cfg


def _ask_ints(prompt: str, default: str) -> List[int]:
    """One or more integers. Re-prompts instead of raising.

    A menu that tracebacks on a typo is a menu you cannot use: the prompt invites a
    comma-separated list, so ``int(get_input(...))`` crashed on exactly what it asked for.
    """
    while True:
        raw = get_input(prompt, default=default)
        try:
            vals = [int(tok) for tok in str(raw).split(",") if tok.strip()]
        except ValueError:
            print_warning(f"not a list of whole numbers: {raw!r}")
            continue
        if vals:
            return vals


def _ask_float(prompt: str, default: str) -> float:
    while True:
        raw = get_input(prompt, default=default)
        try:
            return float(str(raw).strip())
        except ValueError:
            print_warning(f"not a number: {raw!r}")


def _save_run_config(path: str, spec, names) -> None:
    """Write a config.json the non-interactive runner can replay verbatim."""
    import json
    from dataclasses import asdict
    from src.utils.run_export import _commit
    dest = Path(path).expanduser()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(dict(asdict(spec), datasets=list(names),
                                    stage=spec.stage, commit=_commit()),
                               indent=2, default=str) + "\n")
    print_success(f"wrote {dest}")
    print(f"    python scripts/run_sweep.py --config {dest}")
