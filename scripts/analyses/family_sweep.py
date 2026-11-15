"""Moving within the common model on the two volumes (data of Table S2).

TMTM and ST-MTM are other settings of the common model. The AnchorMTM default loses to ST-MTM on TW (and on SNS on
Tangaroa) because the default is a different point of the family, and because a map that follows a 1-D reference
inherits the neighbourhood quality of that reference. This script measures both, on the trees, references, widths and
tracks of the volume rows of Table 1, with the same interval QP (rgmtm.unified.solve_unified_frame):

  reference    x = q: the reference positions themselves (the limit of a map that draws every feature where it is)
  default      (lambda_R, lambda_G, lambda_M, lambda_C) = (4, 1, .5, .1), tau-selected order, cap tau + beta [Table 1, Exact]
  geometry     (0, 1, .5, .1), same order and cap
  uncapped     (0, 1, .5, .1), same order, no position cap
  ST order     (0, 1, .5, .1), no position cap, leaf order of the ST-MTM run, temporal target x_{t-1} as in ST-MTM
  ST-MTM       the stored ST-MTM anchors

All at zero merge error. Weights are the same on both volumes; nothing is tuned. Scores: metrics_n1c.score (SNS, TW k=1, TD, ME).
Run from the repository root:  python scripts/analyses/family_sweep.py [Tangaroa Cylinder]
Output: output/analyses/family_sweep/<name>.json
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
import anchormtm_paths as AP
from experiments.metrics_n1c.run import score
from experiments.order_layout.run import TIE
from rgmtm.framework import interval_costs
from rgmtm.interval_order import width_order
from rgmtm.order_certificate import budget_order
from rgmtm.unified import Controls, solve_unified_frame
from scripts.analyses import volumes

OUT = AP.OUTPUT / "analyses" / "family_sweep"


def sequence(d, orders, R, p, weights, capped, motion_fraction=1.0):
    lr, lg, lm, lc = weights
    rows = []
    for t, o in enumerate(orders):
        q = np.asarray(d["references"][t], float)
        w = p.width_scale * np.asarray(d["measures"][t], float)
        tau = float(interval_costs(np.array([o]), w, q, p)[0])
        prev = prevq = mask = None
        if t:
            old = {k: i for i, k in enumerate(d["tracks"][t - 1])}
            mask = np.array([k in old for k in d["tracks"][t]])
            prev = np.array(
                [rows[-1]["x"][old[k]] if k in old else 0.0 for k in d["tracks"][t]]
            )
            prevq = np.array(
                [
                    d["references"][t - 1][old[k]] if k in old else 0.0
                    for k in d["tracks"][t]
                ]
            )
        c = Controls(
            topology_budget=TIE,
            position_budget=(tau + p.extra_budget) / p.canvas if capped else np.inf,
            reference_weight=lr,
            geometry_weight=lg,
            temporal_weight=lm,
            eccentricity_weight=lc,
            motion_fraction=motion_fraction,
        )
        rows.append(
            solve_unified_frame(
                d["centers"][t],
                d["measures"][t],
                d["hierarchies"][t],
                q,
                d["merge"][t],
                R,
                c,
                p,
                candidate_orders=[o],
                previous=prev,
                previous_q=prevq,
                matched=mask,
            )
        )
    return [r["x"] for r in rows]


def run(name, axis=None, *, refresh_st=False):
    axis = axis or ("motion" if name == "Tangaroa" else "x")
    d, p, feats = volumes.load(name, axis)
    R, L = d["scalar_range"], p.canvas
    folder = volumes.FOLDER[name]
    occupied = name == "Tangaroa"
    if occupied:
        diag = json.loads((folder / "reference.json").read_text())["domain_diagonal"]
    else:
        with np.load(folder / "coordinates.npz", allow_pickle=False) as c:
            diag = float(np.linalg.norm([np.ptp(c[k]) for k in "xyz"]))
    Q = [np.asarray(q, float) for q in d["references"]]

    def scored(xs, reflect=False):
        s = score(d, xs, p, feats, float(diag), occupied)
        npe = lambda V: float(
            np.mean([np.max(np.abs(np.asarray(x) - q)) / L for x, q in zip(V, Q)])
        )
        # as in Table 1, an axis with no direction of its own (ST-MTM order) is reflected once in its favour
        s["NPE"] = (
            min(npe(xs), npe([2 * p.canvas_origin + L - np.asarray(x) for x in xs]))
            if reflect
            else npe(xs)
        )
        return s

    orders = []
    for t, (V, q) in enumerate(
        zip(d["merge"], d["references"]) if not refresh_st else ()
    ):
        w = p.width_scale * np.asarray(d["measures"][t], float)
        orders.append(
            width_order(
                V, q, w, p, 0.0, TIE * R, point=budget_order(V, q, L, 0.0, TIE * R)
            )["order"]
        )
    nat = json.loads((folder / "baseline_landmarks.json").read_text())["positions"][
        "ST-MTM"
    ]
    stored = [p.canvas_origin + L * np.asarray(u, float) for u in nat]
    st_orders = [
        [int(i) for i in np.argsort(np.asarray(u, float), kind="stable")] for u in nat
    ]
    out = dict(
        dataset=name,
        frames=len(Q),
        canvas=L,
        weights="(lambda_R, lambda_G, lambda_M, lambda_C)",
        instances={},
    )
    if refresh_st:
        out = json.loads(
            (AP.FIGDATA / "family_sweep" / f"{name.lower()}.json").read_text()
        )

    def add(tag, xs):
        out["instances"][tag] = scored(xs, tag in ("ST-MTM", "ST order"))
        s = out["instances"][tag]
        print(
            name,
            tag,
            {
                k: (round(v, 4) if isinstance(v, float) else v)
                for k, v in s.items()
                if k in ("SNS", "TW1", "TD", "ME_errors", "ME_clear", "NPE")
            },
            flush=True,
        )

    if not refresh_st:
        add("reference", Q)
    add("ST-MTM", stored)
    for tag, oo, wt, cap, mf in (
        ("default", orders, (4.0, 1.0, 0.5, 0.1), True, 1.0),
        ("geometry", orders, (0.0, 1.0, 0.5, 0.1), True, 1.0),
        ("uncapped", orders, (0.0, 1.0, 0.5, 0.1), False, 1.0),
        ("ST order", st_orders, (0.0, 1.0, 0.5, 0.1), False, 0.0),
    ):
        if refresh_st and tag != "ST order":
            continue
        t0 = time.perf_counter()
        try:
            add(tag, sequence(d, oo, R, p, wt, cap, mf))
        except Exception as e:
            out["instances"][tag] = dict(failed=repr(e))
            print(name, tag, "FAILED", repr(e), flush=True)
        print(name, tag, f"{time.perf_counter() - t0:.0f}s", flush=True)
    return out


def main():
    parser = argparse.ArgumentParser(description="Score volume layouts for Table S2.")
    parser.add_argument("datasets", nargs="*", metavar="DATASET")
    parser.add_argument(
        "--refresh-st", action="store_true", help="Recompute the two ST-MTM rows"
    )
    args = parser.parse_args()
    names = args.datasets or volumes.FOLDER
    for name in names:
        if name not in volumes.FOLDER:
            parser.error(f"unknown dataset {name!r} (choose Cylinder or Tangaroa)")
    OUT.mkdir(parents=True, exist_ok=True)
    for name in names:
        result = run(name, refresh_st=args.refresh_st)
        (OUT / f"{name.lower()}.json").write_text(json.dumps(result, indent=1) + "\n")


if __name__ == "__main__":
    main()
