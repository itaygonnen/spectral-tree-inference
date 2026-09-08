# results/real_data/

Everything the real-cohort experiments produce, one directory per cohort. This is the
whole output: nothing is written beside the notebooks, and nothing else has to be
collected from a cluster run.

```
results/real_data/<cohort>/          e.g. 6000_taxa/
    screen.npz          step 1, machine-readable: one row per tree
    screen.csv          step 1, readable: tree, m, eta_L, valid_L, eta_B, valid_B
    screen.log          what the step-1 run printed
    sweep/<tree>.npz    step 2, per tree: every metric, both arms, one value per p
    sweep_summary.csv    step 2, readable: per p, median and std across trees
    sweep.log           what the step-2 run printed
```

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
