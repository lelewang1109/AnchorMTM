# Planar benchmark protocol

`python -m experiments.planar.run --dataset Ring` runs one sequence. The same entry point supports Wildfires and Storms. Run `python -m experiments.planar.pm25 --budget 0.05 --baselines` for the hourly PM₂.₅ comparison.

- **Input.** The planar persistence fractions are 0.001 for Ring, 0.15 for Wildfires, 0.105 for Storms and 0.05 for hourly PM₂.₅. Each fraction scales the scalar range of its frame. Input acquisition and preprocessing are listed in [Data](../../data/README.md).
- **Exact.** The requested merge budget is zero.
- **Balanced.** Use the first chronological quarter to evaluate budgets from 1% to 10% in steps of 1%. Among budgets that pass the rendered-merge check, minimize the equally weighted sum of min–max normalized NPE and SNS. Within 0.01 of the minimum, prefer the smaller budget, then smaller attained merge error. The hourly PM₂.₅ case uses a fixed 5%.
- **Placement.** The benchmark defaults to Clarabel for the fixed-order QP. `--solver slsqp` selects the archived execution path. Positive-budget benchmark runs use an additional numerical allowance of 10⁻⁷ in the search, and the final map is checked against its nominal budget.
- **Baselines.** The same episode, features and metrics are used for comparison. The planar runner outputs baseline anchors and NPE. The shared formulas for SNS and TW are applied by `scripts/analyses/baseline_sns_tw.py`.
- **Timing.** The timed stages read the source and extract fresh trees. Reported historical timings come from single-threaded runs and are stored separately from new timings.

The benchmark writes to `output/planar/`. The working directory caches features and solved layouts. Choose a fresh `ANCHORMTM_OUTPUT` when changing inputs, parameters or solvers.
