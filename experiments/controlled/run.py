"""Enumerate four-leaf orders and verify displacement bounds with independent LPs."""

import argparse
import csv
import json
import sys
import platform
import scipy
from itertools import permutations
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

CODE = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(CODE), str(CODE / "src")]
from anchormtm_paths import OUTPUT
from rgmtm.error_budget import Parameters
from rgmtm.framework import interval_costs, merge_costs
from rgmtm.interval_order import width_certificate, width_order
from rgmtm.order_certificate import certificate, order_errors

OUT = OUTPUT / "controlled"
ORDERS = np.array(list(permutations(range(4))))


def tree(inner=0.4):
    V = np.ones((4, 4))
    np.fill_diagonal(V, 0)
    V[0, 1] = V[1, 0] = V[2, 3] = V[3, 2] = inner
    return V


def references(d):
    return np.array([0.2, 0.5, 0.5 - d, 0.8 - d])


def lp(order, q, w, p):
    """Independent minimax LP: x[4], centre[4], t, no certificate formula."""
    A, b = [], []

    def add(pairs, rhs):
        row = np.zeros(9)
        for k, v in pairs:
            row[k] += v
        A.append(row)
        b.append(rhs)

    for a, c in zip(order[:-1], order[1:]):
        add([(4 + a, 1), (4 + c, -1)], -(w[a] + w[c]) / 2 - p.gap)
    for i in range(4):
        add([(i, 1), (4 + i, -1)], p.rho * w[i] / 2)
        add([(i, -1), (4 + i, 1)], p.rho * w[i] / 2)
        add([(i, 1), (8, -1)], q[i])
        add([(i, -1), (8, -1)], -q[i])
    bounds = [(0, 1)] * 4 + [(v / 2, 1 - v / 2) for v in w] + [(0, None)]
    r = linprog([0] * 8 + [1], A_ub=A, b_ub=b, bounds=bounds, method="highs")
    assert r.success, r.message
    assert np.max(np.array(A) @ r.x - np.array(b)) < 1e-8
    return float(r.fun), r.x[:4].tolist(), r.x[4:8].tolist()


def barrier_check(order, V, delta):
    barriers = [
        min(V[order[i], order[j]] for i in range(k + 1) for j in range(k + 1, 4))
        + delta
        for k in range(3)
    ]
    measured = max(
        abs(max(barriers[i:j]) - V[order[i], order[j]])
        for i in range(4)
        for j in range(i + 1, 4)
    )
    assert abs(measured - delta) < 1e-10
    assert all(
        barriers[k] >= max(V[order[k], order[k]], V[order[k + 1], order[k + 1]])
        for k in range(3)
    )
    return barriers, measured


def run(output):
    output.mkdir(parents=True, exist_ok=True)
    V = tree()
    delta = merge_costs(ORDERS, V)
    exact = delta < 1e-10
    rows, witnesses, checks, enumeration = [], [], [], []
    widths = [0.0, 0.03, 0.06, 0.09]
    for width in widths:
        w = np.full(4, width)
        p = Parameters(canvas=1.0, gap=0.01 if width else 0.0, rho=0.5)
        for d in np.linspace(-0.12, 0.30, 43):
            q = references(d)
            solutions = [lp(o, q, w, p) for o in ORDERS]
            costs = np.array([s[0] for s in solutions])
            enumeration.append(
                dict(
                    overlap=float(d),
                    width=width,
                    references=q.tolist(),
                    orders=ORDERS.tolist(),
                    lp_tau=costs.tolist(),
                    merge_error=delta.tolist(),
                    lp_anchors=[s[1] for s in solutions],
                    lp_centres=[s[2] for s in solutions],
                )
            )
            closed = interval_costs(ORDERS, w, q, p)
            checks.extend(abs(costs - closed).tolist())
            k = int(np.flatnonzero(exact)[np.argmin(costs[exact])])
            kfree = int(np.argmin(costs))
            nu = certificate(V, q, 1.0)
            lower = width_certificate(V, q, w, p)
            search = width_order(V, q, w, p, 0)
            assert abs(nu - max(0.0, d) / 2) < 1e-10
            assert abs(order_errors(ORDERS[exact], q, 1.0).min() - nu) < 1e-10
            assert abs(lower - costs[k]) < 1e-9
            assert abs(search["tau"] - costs[k]) < 1e-9
            expected = max(0.0, (d + (1 - p.rho) * width + p.gap) / 2)
            assert abs(costs[k] - expected) < 1e-9
            H = float(costs[k] - costs[kfree])
            row = dict(
                overlap=float(d),
                width=width,
                gap=p.gap,
                rho=p.rho,
                nu=nu,
                width_lower=lower,
                exact_tau=float(costs[k]),
                free_tau=float(costs[kfree]),
                hierarchy_excess=H,
                search_tau=search["tau"],
                free_merge_error=float(delta[kfree]),
            )
            rows.append(row)
            b, measured = barrier_check(ORDERS[k], V, delta[k])
            witnesses.append(
                dict(
                    **row,
                    order=ORDERS[k].tolist(),
                    x=solutions[k][1],
                    centres=solutions[k][2],
                    barriers=b,
                    measured_merge_error=measured,
                    free_order=ORDERS[kfree].tolist(),
                )
            )
    assert max(checks) < 1e-9
    budgets = []
    q = references(0.24)
    w = np.full(4, 0.06)
    p = Parameters(canvas=1.0, gap=0.01, rho=0.5)
    costs = np.array([lp(o, q, w, p)[0] for o in ORDERS])
    for inner in [0.2, 0.4, 0.6, 0.8]:
        V = tree(inner)
        delta = merge_costs(ORDERS, V)
        for eps in np.linspace(0, 0.45, 91):
            ok = np.flatnonzero(delta <= eps + 1e-10)
            k = int(ok[np.argmin(costs[ok])])
            r = width_order(V, q, w, p, float(eps))
            assert abs(r["tau"] - costs[k]) < 1e-9
            b, measured = barrier_check(ORDERS[k], V, delta[k])
            budgets.append(
                dict(
                    inner_level=inner,
                    parent_gap=1 - inner,
                    epsilon=float(eps),
                    tau=float(costs[k]),
                    delta=float(delta[k]),
                    order=ORDERS[k].tolist(),
                    barriers=b,
                )
            )
    for name, records in [("overlap", rows), ("budget", budgets)]:
        (output / f"{name}.json").write_text(json.dumps(records, indent=2) + "\n")
        with (output / f"{name}.csv").open("w") as f:
            keys = [k for k, v in records[0].items() if not isinstance(v, list)]
            writer = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(records)
    (output / "enumeration.json").write_text(json.dumps(enumeration, indent=2) + "\n")
    (output / "witnesses.json").write_text(json.dumps(witnesses, indent=2) + "\n")
    summary = dict(
        python=platform.python_version(),
        numpy=np.__version__,
        scipy=scipy.__version__,
        overlap_settings=len(rows),
        budget_settings=len(budgets),
        orders_per_setting=24,
        independent_lp_solves=len(checks) + 24,
        max_lp_formula_difference=max(checks),
        max_width_bound_gap=max(r["exact_tau"] - r["width_lower"] for r in rows),
        max_search_gap=max(abs(r["exact_tau"] - r["search_tau"]) for r in rows),
        reference_example=next(
            r
            for r in rows
            if abs(r["width"] - 0.06) < 1e-10 and abs(r["overlap"] - 0.24) < 1e-10
        ),
        scope="Controlled four-leaf trees; no claim of general width-bound tightness or physical time evolution.",
    )
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    run(parser.parse_args().output)
