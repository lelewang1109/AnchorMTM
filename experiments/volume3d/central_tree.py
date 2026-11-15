"""Run the Cylinder hierarchy-guided variant from full native TMTM profiles."""

import argparse
import json
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import osqp
from scipy import sparse
from scipy.optimize import linprog

import anchormtm_paths as paths
from experiments.grid.run import Frame
from experiments.volume3d.compact_tmtm import barrier_matrix
from rgmtm.error_budget import Parameters, constraints


def clean(a):
    if isinstance(a, dict):
        return {str(k): clean(v) for k, v in a.items()}
    if isinstance(a, (list, tuple)):
        return [clean(v) for v in a]
    if isinstance(a, np.ndarray):
        return a.tolist()
    if isinstance(a, np.generic):
        return a.item()
    return a


def dump(p, a):
    p.write_text(json.dumps(clean(a), indent=2, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "prepared", type=Path, help="Cylinder fields, augmented trees and TMTM maps"
    )
    parser.add_argument("--output", type=Path, default=paths.OUTPUT / "cylinder_guided")
    args = parser.parse_args()
    data = args.prepared.resolve()
    O = args.output.resolve()
    O.mkdir(parents=True, exist_ok=True)
    dataset = paths.RESULTS / "volume_inputs/cylinder"
    F = json.loads((dataset / "features.json").read_text())
    oldp = json.loads((dataset / "parameters.json").read_text())
    coords = np.load(data / "coordinates.npz", allow_pickle=False)
    L, origin = oldp["canvas"], oldp["canvas_origin"]
    H, epsilon = 32768, 0.01
    weights = (1.0, 4.0, 16.0)
    T = len(F["frames"])
    started = time.time()
    p = Parameters(
        canvas=1,
        canvas_origin=-0.5,
        width_scale=oldp["width_scale"] / L,
        gap=oldp["gap"] / L,
        rho=0.95,
        extra_budget=oldp["extra_budget"] / L,
        motion_weight=0.5,
        eccentricity_weight=0.1,
        reference_weight=1,
        geometry_weight=1,
        normalize_terms=True,
    )
    # Fix one global true-X cap, based only on feasibility of the fixed tree orders.
    inputs = []
    maxtau = 0
    for t in range(T):
        ids = np.array(F["frames"][t]["leaf_ids"])
        N = int(np.prod(np.load(data / "source_fields.npy", mmap_mode="r").shape[1:]))
        sample_ids = np.load(data / f"maps/TMTM_ids_{t:03d}.npy")
        assert len(sample_ids) == N and np.all(
            np.bincount(sample_ids, minlength=N) == 1
        )
        inv = np.empty(N, int)
        inv[sample_ids] = np.arange(N)
        guide = inv[ids] / (N - 1) - 0.5
        order = np.argsort(guide)
        q = (np.array(F["frames"][t]["extrema"])[:, 0] - origin) / L - 0.5
        w = p.width_scale * np.array(F["frames"][t]["support_counts"])
        A, b, bounds = constraints(w, order, p)
        n = len(ids)
        E = np.c_[np.eye(n), np.zeros((n, n)), -np.ones(n)]
        G = np.c_[-np.eye(n), np.zeros((n, n)), -np.ones(n)]
        lp = linprog(
            np.r_[np.zeros(2 * n), 1.0],
            A_ub=np.vstack([np.c_[-A, np.zeros(len(A))], E, G]),
            b_ub=np.r_[-b, q, -q],
            bounds=bounds + [(0, None)],
            method="highs",
        )
        assert lp.success
        maxtau = max(maxtau, float(lp.fun))
        inputs.append(
            dict(ids=ids, guide=guide, order=order, q=q, w=w, tau=float(lp.fun))
        )
    cap = max(0.5, maxtau + 0.01)
    dump(
        O / "parameters.json",
        dict(
            persistence_fraction=0.02,
            epsilon=epsilon,
            canvas=[-0.5, 0.5],
            physical_X_reference="(X-original_domain_min)/original_domain_length - 0.5",
            original_domain=[origin, origin + L],
            tree_guidance_weights=weights,
            parameters=p.__dict__,
            global_X_cap_fraction=cap,
            max_order_feasibility_tau=maxtau,
            cap_rule="max(0.5, max fixed-order feasibility tau + 0.01), fixed across sequence",
            tree_order="cached full-sample compact TMTM, temporal focus 0",
            raster_length=H,
            interval_boundary_padding_fraction=1 / H,
            scope="Fixed TMTM leaf order with a spatial position cap.",
        ),
    )
    print("FEASIBILITY", cap, maxtau, flush=True)
    results = {
        g: dict(rows=[], maps=np.empty((H, T), np.float32), checks=[]) for g in weights
    }
    for t, inp in enumerate(inputs):
        ids = inp["ids"]
        guide = inp["guide"]
        order = inp["order"]
        q = inp["q"]
        w = inp["w"]
        n = len(ids)
        A, b, bounds = constraints(w, order, p)
        fr = Frame(t, data / "trees", coords["x"], coords["y"], coords["z"])
        V = np.array(F["frames"][t]["source_lca_values"])
        values = np.load(data / f"maps/TMTM_{t:03d}.npy").astype(float)
        sample_ids = np.load(data / f"maps/TMTM_ids_{t:03d}.npy")
        N = len(values)
        inv = np.empty(N, int)
        inv[sample_ids] = np.arange(N)
        leafpos = inv[ids]
        err0 = float(abs(barrier_matrix(values, leafpos) - V).max())
        assert err0 < 1e-8
        source_intervals = []
        for leaf in ids:
            pp = np.r_[inv[leaf], inv[fr.arc_for_child(leaf).regular_vertices]]
            lo0 = int(pp.min())
            hi0 = int(pp.max())
            assert hi0 - lo0 + 1 == len(pp)
            source_intervals.append((lo0, hi0))
        R = np.c_[np.eye(n), np.zeros((n, n))]
        D = np.c_[-np.eye(n), np.eye(n)]
        pairs = list(combinations(order, 2))
        B = np.zeros((len(pairs), 2 * n))
        centers = np.array(F["frames"][t]["centers"]) / L
        dist = np.array([np.linalg.norm(centers[j] - centers[i]) for i, j in pairs])
        for k, (i, j) in enumerate(pairs):
            B[k, j] = 1
            B[k, i] = -1
        bounds = bounds[:n] + [(a + 1 / H, b0 - 1 / H) for a, b0 in bounds[n:]]
        G = np.vstack([A, R, np.eye(2 * n)])
        lower = np.r_[b, q - cap, [a for a, b0 in bounds]]
        upper = np.r_[np.full(len(A), np.inf), q + cap, [b0 for a, b0 in bounds]]
        for g in weights:
            res = results[g]
            design = [
                np.sqrt(g / n) * R,
                np.sqrt(1 / n) * R,
                np.sqrt(0.1 / n) * D,
                np.sqrt(1 / max(1, len(pairs))) * B,
            ]
            target = [
                np.sqrt(g / n) * guide,
                np.sqrt(1 / n) * q,
                np.zeros(n),
                np.sqrt(1 / max(1, len(pairs))) * dist,
            ]
            matched = []
            old = {}
            if t:
                old = {k: i for i, k in enumerate(F["tracks"][t - 1])}
                matched = [i for i, k in enumerate(F["tracks"][t]) if k in old]
                if matched:
                    oldi = [old[F["tracks"][t][i]] for i in matched]
                    goal = (
                        np.array(res["rows"][-1]["x"])[oldi]
                        + q[matched]
                        - inputs[t - 1]["q"][oldi]
                    )
                    design.append(np.sqrt(0.5 / len(matched)) * R[matched])
                    target.append(np.sqrt(0.5 / len(matched)) * goal)
            M = np.vstack(design)
            target = np.concatenate(target)
            P = sparse.csc_matrix(2 * np.einsum("ki,kj->ij", M, M, optimize=False))
            lin = -2 * np.einsum("ki,k->i", M, target, optimize=False)
            assert np.isfinite(P.data).all() and np.isfinite(lin).all()
            solver = osqp.OSQP()
            solver.setup(
                P=P,
                q=lin,
                A=sparse.csc_matrix(G),
                l=lower,
                u=upper,
                eps_abs=1e-9,
                eps_rel=1e-9,
                max_iter=100000,
                polishing=True,
                verbose=False,
            )
            sol = solver.solve()
            assert sol.info.status_val in [1, 2], (t, g, sol.info.status)
            v = sol.x
            x = v[:n]
            z = v[n:]
            slack = float(
                min(
                    np.min(np.einsum("ij,j->i", G, v, optimize=False) - lower),
                    np.min(upper - np.einsum("ij,j->i", G, v, optimize=False)),
                )
            )
            assert slack > -1e-7, (t, g, slack)
            assert np.max(abs(x - q)) <= cap + 1e-7
            # Monotone warp of all original samples through interval boundaries and anchors.
            cp = [(-1.0, -0.5), (float(N), 0.5)]
            for i in order:
                s, e = source_intervals[i]
                cp.extend(
                    [
                        (s - 0.49, z[i] - w[i] / 2),
                        (float(leafpos[i]), x[i]),
                        (e + 0.49, z[i] + w[i] / 2),
                    ]
                )
            cp = sorted(cp)
            ss = np.array([a for a, b0 in cp])
            dd = np.array([b0 for a, b0 in cp])
            assert np.all(np.diff(ss) > 0) and np.all(np.diff(dd) > 0), (
                t,
                g,
                np.min(np.diff(dd)),
            )
            control_error = float(abs(np.interp(leafpos, ss, dd) - x).max())
            assert control_error < 1e-10
            invgrid = np.interp(np.linspace(-0.5, 0.5, H), dd, ss)
            col = values[np.clip(np.rint(invgrid).astype(int), 0, N - 1)].copy()
            pix = np.rint((x + 0.5) * (H - 1)).astype(int)
            assert len(set(pix)) == n
            col[pix] = np.diag(V)
            for i, j in zip(order[:-1], order[1:]):
                seg = col[pix[i] + 1 : pix[j]]
                assert len(seg) > 0
                assert seg.max() <= V[i, j] + 1e-6, (t, g, float(seg.max()), V[i, j])
                seg[np.argmax(seg)] = V[i, j]
            err = float(abs(barrier_matrix(col, pix) - V).max())
            assert err <= epsilon + 1e-12, (t, g, err)
            spurious = 0
            for i, j in zip(order[:-1], order[1:]):
                seg = col[pix[i] : pix[j] + 1]
                peak = int(np.argmax(seg))
                spurious += int(
                    np.any(np.diff(seg[: peak + 1]) < -1e-7)
                    or np.any(np.diff(seg[peak:]) > 1e-7)
                )
            assert spurious == 0, (t, g, spurious)
            pixel = pix / (H - 1) - 0.5
            me = 0.0
            if matched:
                oldi = [old[F["tracks"][t][i]] for i in matched]
                motion = (
                    pixel[matched]
                    - np.array(res["rows"][-1]["pixel_x"])[oldi]
                    - (q[matched] - inputs[t - 1]["q"][oldi])
                )
                me = float(np.sqrt(np.mean(motion**2)))
            row = dict(
                frame=t,
                x=x,
                z=z,
                w=w,
                order=order,
                reference_X=q,
                tree_guide=guide,
                pixel_x=pixel,
                position_mean=float(np.mean(abs(pixel - q))),
                position_max=float(np.max(abs(pixel - q))),
                tree_guide_RMS=float(np.sqrt(np.mean((pixel - guide) ** 2))),
                ME=me,
                NMLEmax=err,
                min_constraint_slack=slack,
                qp_status=sol.info.status,
            )
            res["rows"].append(clean(row))
            res["maps"][:, t] = col
            res["checks"].append(
                dict(
                    frame=t,
                    source_NMLEmax=err0,
                    raster_NMLEmax=err,
                    anchor_control_error=control_error,
                    original_samples=N,
                    spurious_inter_anchor_excursions=spurious,
                    continuous_position_budget_ok=True,
                    raster_position_cap_with_rounding_ok=bool(
                        np.max(abs(pixel - q)) <= cap + 0.5 / (H - 1) + 1e-7
                    ),
                    raster_position_rounding_tolerance=0.5 / (H - 1),
                )
            )
        if t % 10 == 0:
            print(
                "FRAME",
                t,
                "/",
                T - 1,
                "seconds",
                round(time.time() - started, 1),
                flush=True,
            )
    for g, res in results.items():
        d = O / f"guide_{g:g}"
        d.mkdir(exist_ok=True)
        dump(d / "layouts.json", res["rows"])
        np.savez_compressed(d / "map.npz", map=res["maps"])
        dump(d / "verification.json", res["checks"])
        sm = dict(
            guide_weight=g,
            frames=T,
            NPE_mean=float(np.mean([r["position_mean"] for r in res["rows"]])),
            NPE_max=float(max(r["position_max"] for r in res["rows"])),
            ME_mean=float(np.mean([r["ME"] for r in res["rows"]][1:])),
            NMLEmax=float(max(r["NMLEmax"] for r in res["rows"])),
            initial_anchor_mean=float(np.mean(res["rows"][0]["pixel_x"])),
            position_cap=cap,
            all_frames_verified=True,
        )
        dump(d / "summary.json", sm)
        print("SUMMARY", json.dumps(sm), flush=True)
    dump(
        O / "verification.json",
        dict(
            passed=True,
            frames=T,
            weights=weights,
            epsilon=epsilon,
            hard_X_cap=cap,
            seconds=time.time() - started,
            all_original_samples_used_in_continuous_tree_profile=True,
            uniform_raster_includes_exact_leaf_values_and_separating_LCA_levels=True,
            scope="Fixed-tree-order central guidance prototype; continuous interval widths and true-X constraints; raster certificates concern retained leaf labels, not every augmented node.",
        ),
    )
    print("COMPLETE", flush=True)


if __name__ == "__main__":
    main()
