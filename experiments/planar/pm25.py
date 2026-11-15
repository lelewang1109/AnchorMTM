"""Build and run the January 2017 hourly PM2.5 sequence."""

import argparse
from pathlib import Path

import numpy as np
from scipy.spatial import Delaunay
from threadpoolctl import threadpool_limits

import anchormtm_paths as AP

from . import qp_fallback, run

NAME = "PM25"
BUDGETS = (0.0, 0.01, 0.02, 0.05, 0.10, 0.20)


def prepare(source):
    with np.load(source, allow_pickle=False) as data:
        values, locations = data["pm25"], data["xy"]
    x, y = np.meshgrid(np.linspace(110, 120, 160), np.linspace(32, 41, 144))
    grid = np.column_stack((x.ravel(), y.ravel()))
    triangulation = Delaunay(locations)
    simplex = triangulation.find_simplex(grid)
    if np.any(simplex < 0):
        raise ValueError("The source samples do not cover the 160 × 144 grid")
    transform = triangulation.transform[simplex]
    weights = np.einsum("nij,nj->ni", transform[:, :2], grid - transform[:, 2])
    weights = np.column_stack((weights, 1 - weights.sum(axis=1)))
    fields = (
        np.einsum("tnk,nk->tn", values[:, triangulation.simplices[simplex]], weights)
        .reshape(744, 144, 160)
        .astype(np.float32)
    )
    coordinates = grid * [111.32 * np.cos(np.deg2rad(36.5)), 110.57]
    coordinates -= coordinates.min(axis=0)
    return fields, coordinates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path, default=AP.DATA / "raw/pm25_hourly_source.npz"
    )
    parser.add_argument("--output", type=Path, default=AP.OUTPUT / "pm25_hourly")
    parser.add_argument("--budget", type=float, choices=BUDGETS, default=0.05)
    parser.add_argument("--baselines", action="store_true")
    args = parser.parse_args()
    run.PARAMS[NAME] = dict(nx=160, ny=144, T=744, p=0.05, L=4096, k=1, Start=0)
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    run.O = root
    with threadpool_limits(limits=1):
        fields, coordinates = prepare(args.source)
        trees = root / "trees"
        run.extract(NAME, fields, trees, cache=False)
        area = np.prod(coordinates.max(axis=0) / [159, 143])
        sequence = run.load(NAME, fields, coordinates, area, trees)
        dest = root / f"epsilon_{args.budget:g}"
        with qp_fallback.retry_failed_frames() as certificates:
            result = run.solve(sequence, args.budget)
        summary = run.save_result(dest, NAME, sequence, result, args.budget)
        run.dump(dest / "solver_fallback_certificates.json", certificates)
        if not summary["success"] or summary["rendered_budget_violations"]:
            raise SystemExit(f"PM2.5 failed for budget {args.budget:g}; see {dest}")
        if args.baselines:
            baseline = run.baseline_npe(NAME, sequence)
            if any(not row["success"] for row in baseline):
                raise SystemExit(
                    f"PM2.5 baseline failed; see {root / 'results/PM25/baseline_npe'}"
                )
        print(summary)


if __name__ == "__main__":
    main()
