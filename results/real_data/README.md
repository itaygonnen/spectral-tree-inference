# results/real_data/

Everything the real-cohort experiments produce, one directory per cohort. This is the
whole output: nothing is written beside the notebooks, and nothing else has to be
collected from a cluster run.

```
results/real_data/<cohort>/          e.g. 6000_taxa/
    screen.npz          screening, machine-readable: one row per tree
    screen.csv          screening, readable: tree, m, eta_L, valid_L, eta_B, valid_B
    screen.log          what the screening run printed
    sweep/<tree>.npz    recovery sweep, per tree: every metric, both arms, per p
    sweep_summary.csv   recovery sweep, readable: per p, median and std across trees
    sweep.log           what the sweep run printed
```

Naming a run (the launcher's "Name for this run", or `--prefix`) puts its sweep in
`sweep_<name>/` with `sweep_summary_<name>.csv` and `sweep_<name>.log`, so two grids or
two validity gates sit side by side instead of one invalidating the other. Screening
takes no name: it has no free parameters, so there is one per cohort.

Small by design — a screen is ~10 KB and a swept tree ~3 KB, so the whole thing is a few
hundred KB and **is tracked by git**. A run on another machine therefore comes back with
`git add results/real_data && git commit && git push`; no rsync needed.

`$STR_RESULTS_DIR` moves the root elsewhere (a scratch filesystem, say).

## Reading it

`screen.csv` — one row per tree. `eta` is the partition imbalance (larger clan / smaller
clan) of that operator's split of the full matrix; `valid` is 1 when that split is a real
single-edge bipartition of the true tree. `_L` is the Fiedler vector of `L(S)` cut by
k-means, `_B` is the leading-|λ| eigenvector of `B = HDH` cut by sign.

`sweep_summary.csv` — one row per sub-sampling rate `p`, columns
`<metric>_<arm>_median` and `_std` across the trees included. Metrics: `nmi`, `ari`,
`agreement` (% of taxa on the same side), `dot` (with the reference eigenvector), plus
`sign_L`. NMI is what the paper figure plots; the rest are there so a follow-up question
does not need a re-run.

## Producing it

```bash
cd sub_sampled_fielder_vec
python scripts/interactive_run.py            # real data -> cohorts -> step 1, then step 2
# or, unattended:
nohup python scripts/run_real_sweep.py --cohort "6000 taxa" --stage screen --workers 16 &
nohup python scripts/run_real_sweep.py --cohort "6000 taxa" --stage sweep --cohort-rule valid_S &
```

Both stages are cached per tree and resumable: interrupt either and re-run the same
command. See `sub_sampled_fielder_vec/scripts/cluster/README.md` for the ssh workflow.
