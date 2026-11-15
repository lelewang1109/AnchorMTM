# Tangaroa preprocessing

The full-resolution analysis uses the 201 Tangaroa velocity fields from the [FAU dataset catalogue](https://www.vc.tf.fau.de/publications/datasets/). The native adapter uses TTK 0.9.9.

```sh
python -m experiments.volume3d.setup_inputs
python -m experiments.volume3d.build_native_ttk
python -m experiments.volume3d.prepare_tangaroa_paper --workers 6
python -m experiments.volume3d.paper_trees --threads 6
python -m experiments.volume3d.run_paper
python -m experiments.support.stmtm_completion.run tangaroa
```

`setup_inputs` copies the pinned archive-member inventory. The downloader checks the member sizes and CRCs. The scalar preparation uses the original 300 × 180 × 120 grid and absolute persistence 0.2. The archive, scalar fields and augmented trees require substantial disk space. They are stored beneath `data/real/Tangaroa-FAU/`.

Use `python scripts/analyses/volumes.py run Tangaroa x` for the streamwise comparison on the supplied features or `python scripts/analyses/volumes.py run Tangaroa motion` for the motion-axis variant.
