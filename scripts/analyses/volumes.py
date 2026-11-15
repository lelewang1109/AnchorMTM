"""Volume rows of Table 1: reference axes, Exact and Balanced AnchorMTM, TMTM, ST-MTM and 3-D Hilbert/Morton.

Everything uses the experiment code unchanged:
  layouts   experiments/order_layout/run.py (DATASETS, capped_sequence) with budget_order + width_order, i.e. the
            default configuration of the volume rows of Table 1 (same path as metrics_n1c and volume_gallery_data);
  metrics   experiments/stmtm_completion/run.py geometry() (SNS, TW with k = 1) on the support centroids;
  curves    src/rgmtm/baselines/space_filling.py inverse_order() (canonical 3-D Hilbert/Morton, x fastest);
  bounds    rgmtm.order_certificate.certificate (nu*, strict).
NPE = mean over frames of max |x - q| / L. Curve positions are rank / (N - 1) on the axis (continuous; no raster rounding).
TMTM and ST-MTM keep their stored layouts; only the reference q changes with the axis, and each is reflected once if that
lowers its pooled displacement (as for Table 1). Hilbert/Morton keep their fixed orientation.
Axis candidates: +x, +y, +z of the grid and, for Tangaroa, the stored motion direction (Eq. 2). For an axis e, q = extrema @ e,
the canvas is the projected domain, and every length parameter (width scale, gap, slack) is scaled with the canvas, so the
default rules (width 0.5 L / measure, gap and slack proportional to L) are unchanged.
Balanced: the planar rule -- budgets 1..10% in 1% steps on the first ceil(T/4) frames, J = 0.5 NPE_norm + 0.5 SNS_norm over the
candidates whose attained merge error stays within their budget; within 0.01 the smaller budget, then the smaller attained error.
Usage: python scripts/analyses/volumes.py validate | axes | run NAME AXIS   -> output/analyses/volumes/*.json
"""
import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
import anchormtm_paths as AP
from experiments.order_layout import run as order_layout
from experiments.stmtm_completion.run import geometry
from rgmtm.baselines.space_filling import inverse_order
from rgmtm.error_budget import Parameters
from rgmtm.interval_order import width_order
from rgmtm.order_certificate import budget_order, certificate

OUT = AP.OUTPUT / 'analyses' / 'volumes'
read = lambda p: json.loads(Path(p).read_text())
FOLDER = {name: AP.RESULTS / 'volume_inputs' / name.lower() for name in ('Tangaroa', 'Cylinder')}
SHAPE = {'Tangaroa': (120, 180, 300), 'Cylinder': (48, 64, 135)}                           # z, y, x (x fastest; checked)
BOX = {'Tangaroa': [(-.35, .65), (-.3, .3), (-.5, -.3)], 'Cylinder': [(-2.450261780104711, 20.), (-4., 4.), (0., 6.)]}


def axis_vector(name, axis):
    if axis == 'motion': return np.asarray(read(FOLDER[name] / 'reference.json')['direction'], float)
    return np.eye(3)[{'x': 0, 'y': 1, 'z': 2}[axis]]


def load(name, axis):
    """Dataset of order_layout/run.py with its reference replaced by the projection onto `axis`."""
    p0 = read(FOLDER[name] / 'parameters.json'); e = axis_vector(name, axis)
    corners = np.array([[a, b, c] for a in BOX[name][0] for b in BOX[name][1] for c in BOX[name][2]]) @ e
    o, L = float(corners.min()), float(np.ptp(corners)); ratio = L / p0['canvas']
    saved = read(FOLDER[name] / 'features.json'); feats = saved['frames']
    qs = [list(np.asarray(f['extrema'], float) @ e) for f in feats]
    par = dict(p0, canvas=L, canvas_origin=o, width_scale=p0['width_scale'] * ratio, gap=p0['gap'] * ratio,
               extra_budget=p0['extra_budget'] * ratio)
    measure = 'volumes' if name == 'Tangaroa' else 'support_counts'
    d = dict(centers=[f['centers'] for f in feats], measures=[f[measure] for f in feats],
             hierarchies=[f['hierarchy'] for f in feats], merge=[f['source_lca_values'] for f in feats],
             tracks=saved['tracks'], scalar_range=2.0 if name == 'Tangaroa' else 1.0,
             references=qs, parameters=par)
    return d, Parameters(**{k: v for k, v in par.items() if k in order_layout.PARAM_FIELDS}), feats


def subset(d, n):
    return dict(d, **{k: d[k][:n] for k in ('centers', 'measures', 'hierarchies', 'references', 'merge', 'tracks')})


def nu_mean(d, p):
    R = d['scalar_range']
    return 100 * float(np.mean([certificate(np.asarray(V, float), np.asarray(q, float), p.canvas, 0, 0.)   # strict nu*, as in Fig. 6
                                for V, q in zip(d['merge'], d['references'])]))


def score(d, p, xs):
    L, o = p.canvas, p.canvas_origin; qs = [np.asarray(q, float) for q in d['references']]
    u = [(np.asarray(x, float) - o) / L for x in xs]
    g, _ = geometry(d['centers'], u, d['tracks'], True)
    return dict(SNS=float(g['SNS']), TW=float(g['TW1']), NPE=float(np.mean([np.max(np.abs(np.asarray(x) - q)) / L for x, q in zip(xs, qs)])))


def anchor(d, p, eps):
    R = d['scalar_range']; L = p.canvas; orders = []
    for t, (V, q) in enumerate(zip(d['merge'], d['references'])):
        w = p.width_scale * np.asarray(d['measures'][t], float)
        orders.append(width_order(V, q, w, p, eps * R, order_layout.TIE * R, point=budget_order(V, q, L, eps * R, order_layout.TIE * R)))
    rows = order_layout.capped_sequence(d, orders, R, eps, p)
    s = score(d, p, [r['x'] for r in rows]); s['merge_max'] = 100 * float(max(r['delta'] for r in rows)) / R
    return s, [list(map(float, r['x'])) for r in rows]


def native(name, d, p, feats):
    """TMTM and ST-MTM: stored layouts (their own [0,1] coordinate), scored against this axis."""
    P = read(FOLDER[name] / 'parameters.json'); o0, L0 = P['canvas_origin'], P['canvas']; ids = [f['leaf_ids'] for f in feats]
    T = sorted((r for r in read(FOLDER[name] / 'tmtm.json') if r['method'] == 'TMTM'), key=lambda r: r['index'])
    S = sorted(read(FOLDER[name] / 'stmtm.json'), key=lambda r: r['index'])
    qn = [(np.asarray(q, float) - p.canvas_origin) / p.canvas for q in d['references']]; out = {}
    for meth, rows in (('TMTM', T), ('ST-MTM', S)):
        us = []
        for r, leaves in zip(rows, ids):
            a = (np.asarray(r['positions_unit'], float) if 'positions_unit' in r else
                 np.asarray(r['anchors'], float) / (r['raster'] - 1) if 'raster' in r else (np.asarray(r['positions'], float) - o0) / L0)
            if 'leaf_ids' in r: assert list(r['leaf_ids']) == list(leaves)
            us.append(a)
        e0 = np.mean(np.concatenate([np.abs(u - q) for u, q in zip(us, qn)])); e1 = np.mean(np.concatenate([np.abs(1 - u - q) for u, q in zip(us, qn)]))
        flip = bool(e1 < e0); us = [1 - u if flip else u for u in us]
        g, _ = geometry(d['centers'], us, d['tracks'], True)
        out[meth] = dict(SNS=float(g['SNS']), TW=float(g['TW1']), NPE=float(np.mean([np.max(np.abs(u - q)) for u, q in zip(us, qn)])), reflected=flip)
    return out


def curves(name, d, p, feats):
    qn = [(np.asarray(q, float) - p.canvas_origin) / p.canvas for q in d['references']]; out = {}
    for c in ('Hilbert', 'Morton'):
        t0 = time.perf_counter(); inv = inverse_order(SHAPE[name], c); N = len(inv)
        us = [inv[np.asarray(f['leaf_ids'], int)] / (N - 1) for f in feats]
        g, _ = geometry(d['centers'], us, d['tracks'], True)
        out[c] = dict(SNS=float(g['SNS']), TW=float(g['TW1']), NPE=float(np.mean([np.max(np.abs(u - q)) for u, q in zip(us, qn)])),
                      seconds=time.perf_counter() - t0)
        print(name, c, out[c], flush=True)
    return out


def balanced(name, d, p):
    n = math.ceil(len(d['references']) / 4); ds = subset(d, n); cand = []
    for k in range(1, 11):
        e = k / 100; s, _ = anchor(ds, p, e); s['epsilon'] = e; cand.append(s)
        print(name, 'selection', f'{k}%', {a: round(b, 4) for a, b in s.items()}, flush=True)
    legal = [c for c in cand if c['merge_max'] <= 100 * c['epsilon'] + 1e-6]
    for key in ('NPE', 'SNS'):
        lo, hi = min(c[key] for c in legal), max(c[key] for c in legal)
        for c in legal: c[key + '_norm'] = 0. if hi - lo < 1e-15 else (c[key] - lo) / (hi - lo)
    for c in legal: c['J'] = .5 * c['NPE_norm'] + .5 * c['SNS_norm']
    best = min(c['J'] for c in legal)
    chosen = sorted((c for c in legal if c['J'] <= best + .01), key=lambda c: (c['epsilon'], c['merge_max']))[0]
    return chosen['epsilon'], dict(frames=n, candidates=cand, chosen=chosen['epsilon'])


def main():
    parser = argparse.ArgumentParser(description="Recompute paper layouts on the volume data.")
    commands = parser.add_subparsers(dest="mode", required=True)
    commands.add_parser("validate", help="Check Exact and ten-percent layouts on both volumes")
    commands.add_parser("axes", help="Compute the order bounds for each reference axis")
    for mode, help_text in (
        ("run", "Run the baseline comparison and select the Balanced budget"),
        ("figs", "Recompute data for the volume figure panels"),
    ):
        command = commands.add_parser(mode, help=help_text)
        command.add_argument("dataset", choices=tuple(FOLDER))
        command.add_argument("axis", choices=("x", "y", "z", "motion"))
    args = parser.parse_args()
    if getattr(args, "axis", None) == "motion" and args.dataset != "Tangaroa":
        parser.error("the motion axis is available for Tangaroa")

    OUT.mkdir(parents=True, exist_ok=True)
    mode = args.mode
    if mode == 'validate':   # must reproduce the stored volume runs along the default references
        for name, axis in (('Cylinder', 'x'), ('Tangaroa', 'motion')):
            d, p, feats = load(name, axis)
            if name == 'Tangaroa':
                stored = read(FOLDER[name] / 'reference.json')['qs']
                assert max(np.max(np.abs(np.asarray(a) - np.asarray(b))) for a, b in zip(d['references'], stored)) < 1e-9
            for e in (0, .1):
                t0 = time.perf_counter(); s, _ = anchor(d, p, e)
                print(name, axis, f'eps={e:g}', {k: round(v, 4) for k, v in s.items()}, f'{time.perf_counter() - t0:.0f}s', flush=True)
            print(name, 'native', native(name, d, p, feats), flush=True)
    elif mode == 'axes':
        res = {}
        for name in ('Tangaroa', 'Cylinder'):
            for axis in (['x', 'y', 'z', 'motion'] if name == 'Tangaroa' else ['x', 'y', 'z']):
                d, p, _ = load(name, axis); res.setdefault(name, {})[axis] = dict(nu_mean=nu_mean(d, p), canvas=p.canvas)
                print(name, axis, res[name][axis], flush=True)
        (OUT / 'axes.json').write_text(json.dumps(res, indent=1))
    elif mode == 'figs':      # data of Figs. 6, 7 and 8 for one volume on one axis (same functions as plane_data.py etc.)
        from rgmtm.order_certificate import chain_certificate, order_error
        from rgmtm.framework import merge_costs
        name, axis = args.dataset, args.axis
        d, p, feats = load(name, axis); R = d['scalar_range']; L = p.canvas
        frames = [(np.asarray(V, float), np.asarray(q, float)) for V, q in zip(d['merge'], d['references'])]
        # Fig. 6: bound envelope (strict, as plane_data.py), budgeted point-model orders, projection, baseline order error
        BOUND_EPS = sorted({0., 1e-8, 1e-7, 1e-6, 1e-5} | set(np.geomspace(1e-4, 1, 121).tolist()))
        BUILD_EPS = sorted({0., 1e-7, 1e-6, 1e-5, 1e-4} | set(np.geomspace(1e-3, 1, 31).tolist()) | {.005, .01, .02, .05, .1, .2, .5, 1.})
        lower = np.array([[max(certificate(V, q, L, 2 * e * R), chain_certificate(V, q, L, e * R, 0.)) for V, q in frames] for e in BOUND_EPS])
        mean = lower.mean(1); env = np.maximum.accumulate(mean[::-1])[::-1]; build = []
        for e in BUILD_EPS:
            rows = [budget_order(V, q, L, e * R, 0.) for V, q in frames]
            build.append(dict(epsilon=e, NROE_mean=float(np.mean([r['nroe'] for r in rows])), NMLE_max=float(max(r['delta'] for r in rows) / R)))
        proj = max(float(merge_costs(np.array([np.argsort(q, kind='stable')]), V)[0]) for V, q in frames) / R
        nat = {}
        P = read(FOLDER[name] / 'parameters.json'); o0, L0 = P['canvas_origin'], P['canvas']
        T = sorted((r for r in read(FOLDER[name] / 'tmtm.json') if r['method'] == 'TMTM'), key=lambda r: r['index'])
        S = sorted(read(FOLDER[name] / 'stmtm.json'), key=lambda r: r['index'])
        for meth, rows in (('TMTM', T), ('ST-MTM', S)):
            us = [np.asarray(r['positions_unit'], float) if 'positions_unit' in r else np.asarray(r['anchors'], float) / (r['raster'] - 1)
                  if 'raster' in r else (np.asarray(r['positions'], float) - o0) / L0 for r in rows]
            vals = [float(np.mean([order_error(np.argsort(s * u, kind='stable'), q, L) for u, (V, q) in zip(us, frames)])) for s in (1, -1)]
            nat[meth] = dict(as_stored=vals[0], reversed=vals[1], best=min(vals), positions=[list(map(float, u)) for u in us])
        plane = dict(dataset=name, axis=axis, frames=len(frames), extent=L, scalar_range=R, nu_strict=float(mean[0]),
                     bound=dict(epsilon=BOUND_EPS, mean=mean.tolist(), envelope=env.tolist()), construction=build,
                     projection_NMLE_max=proj, baselines_NROE={m: {k: v for k, v in x.items() if k != 'positions'} for m, x in nat.items()})
        (OUT / f'{name.lower()}_{axis}_plane.json').write_text(json.dumps(plane))
        print(name, axis, 'nu*', round(100 * mean[0], 3), 'baseline NROE', {m: round(100 * x['best'], 2) for m, x in nat.items()}, 'proj', round(100 * proj, 1), flush=True)
        # Figs. 7 and 8: AnchorMTM at five budgets, per-frame worst displacement and merge error
        qs = [np.asarray(q, float) for q in d['references']]; gallery = dict(dataset=name, axis=axis, canvas=L, origin=p.canvas_origin,
                                                                           q=[list(map(float, q)) for q in qs], maps={})
        budget = dict(epsilon=[], median=[], p90=[], max=[], mean=[], used=[])
        for e in (0., .02, .05, .1, .2):
            R_ = R; orders = []
            for t, (V, q) in enumerate(frames):
                w = p.width_scale * np.asarray(d['measures'][t], float)
                orders.append(width_order(V, q, w, p, e * R_, order_layout.TIE * R_, point=budget_order(V, q, L, e * R_, order_layout.TIE * R_)))
            rows = order_layout.capped_sequence(d, orders, R_, e, p)
            E = np.array([100 * np.max(np.abs(np.asarray(r['x']) - q)) / L for r, q in zip(rows, qs)]); dm = np.array([r['delta'] / R_ for r in rows])
            for k, v in (('epsilon', e), ('median', float(np.median(E))), ('p90', float(np.percentile(E, 90))), ('max', float(E.max())),
                         ('mean', float(E.mean())), ('used', float(100 * np.mean(dm > 1e-5)))): budget[k].append(v)
            gallery['maps'][f'eps={e:g}'] = dict(x=[list(map(float, r['x'])) for r in rows], NMLE_max=float(dm.max()))
            print(name, axis, f'eps={e:g}', 'mean E', round(E.mean(), 3), 'used', round(budget['used'][-1], 1), 'merge', round(100 * dm.max(), 2), flush=True)
        for meth, x in nat.items():   # stored layouts, oriented to lower pooled displacement (as volume_gallery_data.native)
            qn = [(q - p.canvas_origin) / L for q in qs]; us = [np.asarray(u) for u in x['positions']]
            flip = np.mean(np.concatenate([np.abs(1 - u - q) for u, q in zip(us, qn)])) < np.mean(np.concatenate([np.abs(u - q) for u, q in zip(us, qn)]))
            gallery['maps'][meth] = dict(x=[list(map(float, p.canvas_origin + (1 - u if flip else u) * L)) for u in us], NMLE_max=0., mirrored=bool(flip))
        (OUT / f'{name.lower()}_{axis}_gallery.json').write_text(json.dumps(gallery))
        (OUT / f'{name.lower()}_{axis}_budget.json').write_text(json.dumps(budget))
        print('saved figs', name, axis)
    elif mode == 'run':
        name, axis = args.dataset, args.axis
        d, p, feats = load(name, axis)
        out = dict(dataset=name, axis=axis, canvas=p.canvas, origin=p.canvas_origin, nu_mean=nu_mean(d, p), maps={})
        out['baselines'] = native(name, d, p, feats); out['curves'] = curves(name, d, p, feats)
        eb, sel = balanced(name, d, p); out['balanced_selection'] = sel
        for tag, e in (('Exact', 0.), ('Balanced', eb), ('10%', .1)):
            if tag == '10%' and abs(eb - .1) < 1e-12: out['maps']['10%'] = out['maps']['Balanced']; continue
            t0 = time.perf_counter(); s, xs = anchor(d, p, e); s['epsilon'] = e; out['maps'][tag] = dict(s, x=xs)
            print(name, axis, tag, {k: round(v, 4) for k, v in s.items()}, f'{time.perf_counter() - t0:.0f}s', flush=True)
        out['references'] = [list(map(float, q)) for q in d['references']]
        (OUT / f'{name.lower()}_{axis}.json').write_text(json.dumps(out))
        print('saved', name, axis)


if __name__ == '__main__':
    main()
