# AnchorMTM

Spatially Anchored Merge Tree Maps with Topology–Position Tradeoffs

AnchorMTM lays out merge trees from a time-varying scalar field along a fixed spatial axis. A merge budget controls how far the displayed merge heights may depart from the input trees. The evaluation covers Ring, Wildfires, Storms, hourly PM₂.₅, Cylinder and Tangaroa.

## Quick start

Use Python 3.11 or 3.12. Figure export requires Poppler (`pdftocairo`) and either Tectonic or pdfLaTeX with TikZ and `standalone`.

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[solvers,volume3d]'
python -m unittest discover -s tests -q
python scripts/run_all.py
```

The last command generates the paper's figures and tables in `output/` and checks the data tables against the paper. It uses the included experiment records and figure sources. See the [reproduction guide](docs/reproduction.md) for running experiments, checking figure exports, and separating measured runtimes from recorded runtimes.

## Repository layout

| Path | Contents |
| --- | --- |
| `src/rgmtm/`, `analysis/` | Merge trees, spatial reference, order certificates, placement and scalar reconstruction |
| `experiments/` | Controlled trees, planar and volumetric experiments, and baselines |
| `results/` | Feature and layout records, intermediate results and paper tables |
| `scripts/` | Analysis, figure generation, table generation and verification |
| `figures_src/` | Editable figures, publication exports and LaTeX overlays |
| `data/` | [Dataset sources, expected paths and preprocessing](data/README.md) |

Author information is in [CITATION.cff](CITATION.cff). The source and experiment records include material from [ydvislab/AnchorMTM](https://github.com/ydvislab/AnchorMTM).
