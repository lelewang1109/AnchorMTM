# Cylinder preprocessing

The square-cylinder dataset has 102 velocity fields. The preparation script computes Q on the original grid, keeps streamwise indices 57–191, and extracts the merge trees.

```sh
python -m experiments.volume3d.build_native_ttk
python -m experiments.grid.acquire_cylinder
python -m experiments.grid.prepare_3d
python -m experiments.grid.run square-cylinder3d
python -m experiments.support.stmtm_completion.run square-cylinder3d
```

The acquisition script reads the [KTH archive](https://tinoweinkauf.net/notes/squarecylinder.html) and extracts it beneath `data/real/SquareCylinder/`. Prepared fields and trees are under `data/real/paper1_suite/square-cylinder3d/`. The baseline runs write to `results/paper1_suite/` and `results/stmtm_completion20261026/`.

Run `python scripts/analyses/volumes.py run Cylinder x` for the streamwise comparison on the supplied features. Run `python -m experiments.volume3d.central_tree data/real/paper1_suite/square-cylinder3d` for the hierarchy-guided case after generating the TMTM profiles.
