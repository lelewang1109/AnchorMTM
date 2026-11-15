# Datasets

| Dataset | Input and preprocessing | Frames |
| --- | --- | ---: |
| Ring | Spreading-ring generator in `experiments/ring/source/`, evaluated on a 14 × 14 grid | 40 |
| Wildfires | Fire radiative power from [Zenodo](https://zenodo.org/records/11234747), accumulated over equal-duration bins from November–December 2019 on a 64 × 64 grid | 32 |
| Storms | ERA5 hourly mean sea-level pressure for December 1999. Subtract the 192-hour moving mean, then smooth with a Gaussian of 2 grid cells | 744 |
| PM₂.₅ | January 2017 hourly CAQRA, northern China. Interpolate onto a 160 × 144 grid over 110–120°E and 32–41°N | 744 |
| Cylinder | Square-cylinder velocity archive. Compute Q on the 192 × 64 × 48 grid and crop to 135 × 64 × 48 | 102 |
| Tangaroa | FAU velocity archive. Use the 300 × 180 × 120 grid, treating near-zero velocity as inflow and clipping speeds to [0, 2] | 201 |

Use the following paths under the data directory:

- Wildfires: `real/wildfire/wildfire.json`.
- Storms: `real/ERA5_MSLP/ERA5_MSLP_19991117_20000114.nc`. Request hourly mean sea-level pressure (30–75°N, 30°W–40°E, 0.25°) for 17 November 1999 through 14 January 2000 from the Copernicus Climate Data Store. The extra days supply the moving-mean margins. The netCDF variables are `msl`, `valid_time`, `latitude` and `longitude`.
- PM₂.₅: `raw/pm25_hourly_source.npz`. Prepare the ChinaVis 2021 CAQRA hourly CSVs with `python -m experiments.planar.prepare_pm25_source CSV_DIRECTORY data/raw/pm25_hourly_source.npz`. The regional subset includes a 0.5° interpolation margin.
- Cylinder: `real/SquareCylinder/SquareCylinder.7z`. See the [preparation commands](../experiments/grid/README.md) for the archive and extraction.
- Tangaroa: `real/Tangaroa-FAU/`. See the [preparation commands](../experiments/volume3d/README.md) for archive acquisition and the pinned member inventory.

`python scripts/hash_inputs.py` checks the three external planar source files against the recorded SHA-256 values. Set `ANCHORMTM_DATA` when these files are stored elsewhere.

The earlier temporal audit uses a 12-hourly ERA5 series from ARCO-ERA5, NOAA daily OLR and daily CN-Reanalysis PM₂.₅ (DOI 10.57760/sciencedb.35623). Their preparation scripts are `experiments/analysis/era5_arco.py`, `experiments/analysis/era5_arco_expanded.py`, `scripts/fetch_public_scalar.py` and `scripts/prepare_cn_pm25.py`.

## Prepared records

`results/inputs/layout_records/` stores the labeled pairwise merge heights, spatial anchors, leaf widths and per-budget layouts used for the bounds and ablation analyses. The `V`, `q` and `w` arrays in each `per_frame_inputs.npz` are split by the frame offsets in the archive. `metadata.json` gives the canvas and merge-height scale. The Tangaroa basin-depth audit uses the motion axis. The main volume comparison uses the streamwise x-axis.

`results/inputs/volume_inputs/` stores extracted leaf geometry, merge heights, tracks and baseline positions. `results/figures/storms_event/` and `results/figures/pm25_baseline_profiles/` contain the baseline profiles used in the case figures.
