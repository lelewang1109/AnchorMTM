"""Build the paper figures and tables from the bundled results."""

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEPS = [
    ("Figure 3", ["scripts/figures/fig03_selective_flattening.py"]),
    ("Figure 4", ["scripts/figures/fig04_forced_error.py"]),
    ("Figure 5", ["scripts/figures/fig05_volume_positions.py"]),
    ("Figure 6", ["scripts/figures/fig06_budget_response.py"]),
    ("Figure 7 map panels (check)", ["scripts/figures/fig07_storms_budget_maps.py"]),
    ("Figures 1 and 7 (East arrow)", ["scripts/figures/east_arrows.py"]),
    ("Figure S4", ["scripts/figures/fig_s04_budget_sweep.py"]),
    ("Figure S5 panels (check)", ["scripts/figures/fig_s05_cylinder_variant.py"]),
    ("Figure S6", ["scripts/figures/fig_s06_cylinder_compare.py"]),
    ("Figure S7", ["scripts/figures/fig_s07_pm25_compare.py"]),
    ("Table 1", ["scripts/tables/table_01.py"]),
    ("Table S1", ["scripts/tables/table_s01.py"]),
    ("Tables S2-S4", ["scripts/tables/table_s02_s04.py"]),
    ("Tables S5–S8", ["scripts/tables/table_s05_s08.py"]),
    ("Paper tables", ["scripts/check_submission.py"]),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--supplement-matplotlib",
        type=Path,
        help="Directory containing Matplotlib 3.11.2 for Figures S4, S6 and S7",
    )
    args = parser.parse_args()
    steps = list(STEPS)
    steps.append(("PNG export", ["scripts/figures/export_png.py"]))
    failed = []
    for name, cmd in steps:
        env = os.environ.copy()
        if args.supplement_matplotlib and name in {
            "Figure S4",
            "Figure S6",
            "Figure S7",
        }:
            paths = [str(args.supplement_matplotlib.resolve())]
            if env.get("PYTHONPATH"):
                paths.append(env["PYTHONPATH"])
            env["PYTHONPATH"] = os.pathsep.join(paths)
        r = subprocess.run(
            [sys.executable, *cmd], cwd=ROOT, env=env, capture_output=True, text=True
        )
        print(f"{'ok  ' if r.returncode == 0 else 'FAIL'} {name}")
        if r.returncode:
            failed.append(name)
            print(r.stderr[-1500:])
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
