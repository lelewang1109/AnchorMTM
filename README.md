# AnchorMTM

Spatially Anchored Merge Tree Maps with Topology–Position Tradeoffs

AnchorMTM lays out merge trees from a time-varying scalar field along a fixed spatial axis. A merge budget controls how far the displayed merge heights may depart from the input trees. The evaluation covers Ring, Wildfires, Storms, hourly PM₂.₅, Cylinder and Tangaroa.

[Paper (PDF)](https://github.com/lelewang1109/AnchorMTM/releases/download/paper-artifacts/AnchorMTM-paper.pdf) · [Supplementary material (PDF)](https://github.com/lelewang1109/AnchorMTM/releases/download/paper-artifacts/AnchorMTM-supplement.pdf) · [Supplementary video (MP4)](https://github.com/lelewang1109/AnchorMTM/releases/download/paper-artifacts/AnchorMTM-supplementary-video.mp4)

## Run a small example

Use Python 3.11 or 3.12. Run these commands from the repository root.

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -e .
python -m experiments.controlled.run
```

The example enumerates the leaf orders of a four-leaf tree and checks the displacement bounds against independent linear programs. It writes CSV files and a summary to `output/controlled/`. It needs no datasets or native TTK build.

## Build the paper figures and tables

Install Poppler (`pdftocairo`) and either Tectonic or pdfLaTeX with TikZ and `standalone`, then run

```sh
python scripts/run_all.py
```

This uses the bundled results and figure sources. It writes figures to `output/figures/` and tables to `output/tables/`, and checks the tables against the paper. For raw-data experiments, optional dependencies and tests, see the [reproduction guide](REPRODUCTION.md).

## Read the implementation

Start with [`experiments/analysis/anchormtm.py`](experiments/analysis/anchormtm.py). `prepare_fields` extracts frame data, and `anchor_mtm` runs the layout and reconstruction. Follow one frame through these steps:

1. [Order certificates and candidates](src/rgmtm/order_certificate.py) compute bounds and propose leaf orders. [Interval order](src/rgmtm/interval_order.py) accounts for leaf widths.
2. [Budgeted placement](src/rgmtm/unified.py) passes the chosen order to the [placement solver](src/rgmtm/error_budget.py), which applies the displacement cap.
3. [Map reconstruction](src/rgmtm/reference_guided.py) turns the layout into scalar columns.

To trace a paper result, start with [Table 1's script](scripts/tables/table_01.py), which lists the records it reads from `results/`. The [planar benchmark](experiments/planar/run.py) and [volume analysis](scripts/analyses/volumes.py) rerun the layouts. See [Data](data/README.md) for input sources and preparation.

| Directory | Contents |
| --- | --- |
| `src/rgmtm/` | Core algorithms and baselines |
| `experiments/` | Controlled, planar and volume experiments, with shared analysis and support code |
| `scripts/` | Reproduction commands, plotting and table generation |
| `data/` | Input instructions and editable figure sources |
| `results/` | Prepared inputs and measurements used by the paper |
| `tests/` | Code tests |

Author information is in [CITATION.cff](CITATION.cff). The source and experiment records include material from [ydvislab/AnchorMTM](https://github.com/ydvislab/AnchorMTM).
