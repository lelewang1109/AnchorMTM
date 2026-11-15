# Reproduction

Run the commands from the repository root in the environment described in the [README](../README.md).

Install the solvers and volume tools before running the full test suite or the dataset experiments.

```sh
python -m pip install -e '.[solvers,volume3d]'
python -m unittest discover -s tests -q
```

## Figures and tables

```sh
python scripts/run_all.py
```

The workflow uses the editable sources and publication exports in `figures_src/`. Figure 2 uses its publication PNG. To export the editable slide again, use `python scripts/figures/fig02_workflow.py SLIDE.pdf` with a PDF export of the slide, or omit the argument when LibreOffice is available. Figures 1 and 7 add their spatial-axis arrows in LaTeX.

Table 1 preserves the recorded baseline TW values for Cylinder TMTM (0.8017) and Tangaroa ST-MTM (0.8137). Direct volume rescoring can differ in the fourth decimal. For example, the Tangaroa motion-axis Exact check measured TW 0.7262 from the prepared inputs, compared with 0.7260 in Table S4. Table S1 exports the planar settings and measured runtimes. The volume timing and parameter note is preserved in [the paper table](../results/paper_tables/table2_parameters_runtime.tex). Table 2 defines metrics and contains no experimental measurements.

Figures S4, S6 and S7 use Matplotlib 3.11.2 in the publication environment. Install this version in a separate directory to retain the default Matplotlib 3.9.4 for the other figures:

```sh
python -m pip install --no-deps --target tmp/matplotlib-supplement matplotlib==3.11.2
python scripts/run_all.py --supplement-matplotlib tmp/matplotlib-supplement
python scripts/check_submission.py /path/to/paper-source --figures
```

The last command compares the regenerated PNGs against the paper source directory. On the validation machine, 13 figures matched pixel for pixel. Figure S7 had a mean absolute RGB difference of 0.218 on a scale of 255.

## Planar experiments

Use the [dataset inputs](../data/README.md) and build the TTK adapter. CMake, Ninja and a C++ compiler with OpenMP are required.

```sh
python -m experiments.volume3d.build_native_ttk
python -m experiments.planar.run --dataset Ring
python -m experiments.planar.run --dataset Wildfires
python -m experiments.planar.run --dataset Storms
python -m experiments.planar.pm25 --budget 0.05 --baselines
```

The planar benchmark uses Clarabel for AnchorMTM's capped placement QP. `--solver slsqp` reproduces the solver path of the archived planar measurements. Use a separate output directory to change the solver. Hourly PM₂.₅ uses SLSQP with a same-QP Clarabel retry on numerical failure. Table S1 retains the timings measured during the recorded runs.

`ANCHORMTM_DATA` points to the dataset root for the planar loaders. `ANCHORMTM_TTK_BIN` selects an existing executable. `ANCHORMTM_OUTPUT` specifies the destination for new planar runs, analyses and figures. The [planar protocol](../experiments/planar/PROTOCOL.md) describes budget selection, map validation and baseline scoring.

## Volume layouts and supplementary analyses

The prepared feature data, tracks, merge matrices and baseline anchors in `results/volume_inputs/` support layout recomputation without the original volumes.

```sh
python scripts/analyses/volumes.py run Cylinder x
python scripts/analyses/volumes.py run Tangaroa x
python scripts/analyses/volumes.py run Tangaroa motion
python scripts/analyses/family_sweep.py Cylinder
python scripts/analyses/tangaroa_x_extras.py
python scripts/analyses/supplement_analyses.py strict
python scripts/analyses/supplement_analyses.py ties
python scripts/analyses/supplement_analyses.py motion
for name in Ring Wildfires Storms PM25 Cylinder Tangaroa; do
    python scripts/analyses/supplement_analyses.py ablation "$name"
    python scripts/analyses/supplement_analyses.py depth "$name"
done
```

For a shorter check, `python scripts/analyses/volumes.py validate` runs Exact and 10% layouts along the Cylinder streamwise axis and the Tangaroa motion axis. The analytical reruns write to `output/analyses/`. The figure and table scripts use the bundled records in `results/`. The source-data preprocessing and baseline pipelines for [Cylinder](../experiments/paper1_suite/README.md) and [Tangaroa](../experiments/volume3d/README.md) require the original flow volumes and additional disk space. The Cylinder hierarchy-guided layout can be rerun from the prepared trees and TMTM profiles with `python -m experiments.volume3d.central_tree PREPARED_DIR`.

The earlier planar temporal audit in Table S3 uses separate OLR, daily PM₂.₅ and pressure inputs. The paper tables use its bundled metric records. The full raw-volume acquisition-to-rendering workflows and this earlier raw-data audit have separate input requirements from the layout reruns described here.
