"""Generate Tables S2 and S3 on the streamwise axis and Table S4 on the motion axis."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import anchormtm_paths as paths

FAMILY = paths.FIGDATA / "family_sweep"
METRICS = paths.RESULTS / "metrics_n1c"
FAMILY_ROWS = [
    ("reference", r"Reference, $x=q$", "--", "--", "--"),
    (
        "default",
        r"AnchorMTM default (Table~1, Exact)",
        r"$(4,1,.5,.1)$",
        r"$\tau$-selected",
        "yes",
    ),
    ("geometry", r"No reference term", r"$(0,1,.5,.1)$", r"$\tau$-selected", "yes"),
    (
        "uncapped",
        "No reference term, no cap",
        r"$(0,1,.5,.1)$",
        r"$\tau$-selected",
        "no",
    ),
    (
        "ST order",
        "ST-MTM leaf order, stationary target",
        r"$(0,1,.5,.1)$",
        "ST-MTM",
        "no",
    ),
    ("ST-MTM", "ST-MTM run", "--", "ST-MTM", "--"),
]
SEQUENCES = [
    ("Storms", "storms", "Storms"),
    ("Wildfires", "wildfires", "Wildfires"),
    ("OLR", "olr", "OLR"),
    ("PM2.5", "pm25", r"PM$_{2.5}$ (daily)"),
    ("Ring", "ring", "Ring"),
    ("Cylinder", "cylinder", "Cylinder"),
    ("Tangaroa", "tangaroa", "Tangaroa"),
]


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sns_tw_npe(scores):
    return f"{scores['SNS']:.3f} & {scores['TW1']:.3f} & {100 * scores['NPE']:.1f}"


def motion_error(scores):
    return 100 * scores["ME_errors"] / scores["ME_clear"]


def format_temporal(scores, metric):
    if metric == "ME":
        return f"{motion_error(scores):.0f}"
    if metric == "TW1":
        return f"{scores[metric]:.3f}"
    return f"{scores[metric]:.1f}"


def write_family_table(output, cylinder, tangaroa):
    """Write Table S2 from the two streamwise-axis family sweeps."""
    lines = [
        r"\begin{table*}[!htb]",
        r"\centering",
        r"\caption{Moving within the common model of Eq.~6 at zero budget (merge levels within $2\varepsilon_0$ count as equal, as in Table~1). Same trees, references, widths and tracks as the volume rows of Table~1; the weights are $(\lambda_R,\lambda_G,\lambda_M,\lambda_C)$ and are the same on both volumes. NPE in \% of the axis.}",
        r"\label{tab:family}",
        r"\footnotesize",
        r"\begin{tabular}{@{}llll ccc ccc@{}}",
        r"\toprule",
        r" & & & & \multicolumn{3}{c}{Cylinder} & \multicolumn{3}{c}{Tangaroa} \\",
        r"\cmidrule(lr){5-7}\cmidrule(l){8-10}",
        r"Setting & Weights & Leaf order & Cap & SNS & TW & NPE & SNS & TW & NPE \\",
        r"\midrule",
    ]
    for k, name, weights, order, cap in FAMILY_ROWS:
        lines.append(
            f"{name} & {weights} & {order} & {cap} & {sns_tw_npe(cylinder[k])} & {sns_tw_npe(tangaroa[k])} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
    (output / "table_family.tex").write_text("\n".join(lines) + "\n")


def write_temporal_table(output, references, metrics_by_sequence, tangaroa_temporal):
    """Write Table S3, using the streamwise scores for Tangaroa AnchorMTM."""
    lines = [
        r"\begin{table*}[!htb]",
        r"\centering",
        r"\caption{Neighborhoods and temporal behavior. $q$ draws every feature at its reference. TW uses $k=1$; TD is the summed displayed motion of matched anchors, in axis lengths; ME is the share (\%) of feature pairs whose displayed approach or separation contradicts the true one. For the volumes, the TMTM, ST-MTM and Exact columns are the runs of Table~1, and the 10\% column comes from separate runs; the planar rows are the sequences of an earlier audit, preprocessed differently from the benchmark reconstructions of Table~1.}",
        r"\label{tab:temporal}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{2.5pt}",
        r"\begin{tabular}{@{}l ccccc ccccc ccccc@{}}",
        r"\toprule",
        r" & \multicolumn{5}{c}{TW$\uparrow$} & \multicolumn{5}{c}{TD} & \multicolumn{5}{c}{ME (\%)$\downarrow$} \\",
        r"\cmidrule(lr){2-6}\cmidrule(lr){7-11}\cmidrule(l){12-16}",
        r"Sequence" + r" & $q$ & TMTM & ST-MTM & Exact & 10\%" * 3 + r" \\",
        r"\midrule",
    ]
    for name, _, label in SEQUENCES:
        metrics = metrics_by_sequence[name]
        if name == "Tangaroa":
            exact, relaxed = tangaroa_temporal["Exact"], tangaroa_temporal["10%"]
        else:
            exact, relaxed = (
                metrics["settings"]["eps=0"],
                metrics["settings"]["eps=0.1"],
            )
        row = [label]
        samples = (
            references[name],
            metrics["baselines"]["TMTM"],
            metrics["baselines"]["ST-MTM"],
            exact,
            relaxed,
        )
        for kind in ("TW1", "TD", "ME"):
            row.extend(format_temporal(sample, kind) for sample in samples)
        lines.append(" & ".join(row) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
    (output / "table_temporal.tex").write_text("\n".join(lines) + "\n")


def write_motion_table(output, reference, metrics):
    """Write Table S4 using the recorded motion-axis scores."""
    curves = {"Hilbert": (0.3429, 0.9222), "Morton": (0.3627, 0.8893)}
    recorded_scores = {
        "TMTM": (0.5255, 0.6701, 0.3063),
        "ST-MTM": (0.1655, 0.8137, 0.3418),
        "Exact": (0.3621, 0.7260, 0.1118),
        "10%": (0.3408, 0.7540, 0.0615),
    }
    rows = [
        (
            "Reference, $x=q$",
            f"{reference['SNS']:.4f} & {reference['TW1']:.4f} & 0 & {reference['TD']:.1f} & {motion_error(reference):.0f}",
        )
    ]
    rows += [
        (k, f"{a:.4f} & {b:.4f} & -- & -- & --") for k, (a, b) in curves.items()
    ]  # no NPE for the curves: an index has no spatial direction
    for k, s in (
        ("TMTM", metrics["baselines"]["TMTM"]),
        ("ST-MTM", metrics["baselines"]["ST-MTM"]),
        ("Exact", metrics["settings"]["eps=0"]),
        ("10%", metrics["settings"]["eps=0.1"]),
    ):
        a, b, c = recorded_scores[k]
        label = {"Exact": "AnchorMTM, Exact", "10%": r"AnchorMTM, Balanced (10\%)"}.get(
            k, k
        )
        rows.append(
            (
                label,
                f"{a:.4f} & {b:.4f} & {c:.4f} & {s['TD']:.1f} & {motion_error(s):.0f}",
            )
        )
    lines = [
        r"\begin{table}[!h]",
        r"\centering",
        r"\caption{Tangaroa with the reference along its motion direction (Eq.~2) instead of the streamwise axis. TD and ME as in Table~\ref{tab:temporal}.}",
        r"\label{tab:tangaroa-motion}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{3pt}",
        r"\begin{tabular}{@{}lccccc@{}}",
        r"\toprule",
        r"Method & SNS$\downarrow$ & TW$\uparrow$ & NPE$\downarrow$ & TD & ME (\%)$\downarrow$ \\",
        r"\midrule",
    ]
    lines += [f"{a} & {b} \\\\" for a, b in rows]
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\vspace{2pt}",
        r"\begin{minipage}{0.94\columnwidth}\footnotesize On this axis the Balanced rule selects 10\% (attained merge error 9.97\%). SNS, TW, TD and ME of TMTM, ST-MTM and the curves do not depend on the axis.\end{minipage}",
        r"\end{table}",
    ]
    (output / "table_tangaroa_motion.tex").write_text("\n".join(lines) + "\n")


def main():
    streamwise = read_json(paths.FIGDATA / "volumes/tangaroa_x_extras.json")
    cylinder = read_json(FAMILY / "cylinder.json")["instances"]
    tangaroa = streamwise["family"]
    motion_reference = read_json(FAMILY / "tangaroa.json")["instances"]["reference"]

    metrics = {
        name: read_json(METRICS / f"{filename}.json") for name, filename, _ in SEQUENCES
    }

    # The ST-MTM comparison uses the same recorded SNS and TW as Table 1.
    for name, family in (("Cylinder", cylinder), ("Tangaroa", tangaroa)):
        baseline = metrics[name]["baselines"]["ST-MTM"]
        family["ST-MTM"].update(SNS=baseline["SNS"], TW1=baseline["TW1"])

    references = read_json(FAMILY / "temporal_reference.json")
    references.update(Cylinder=cylinder["reference"], Tangaroa=tangaroa["reference"])
    write_family_table(paths.TABLES_OUT, cylinder, tangaroa)
    write_temporal_table(paths.TABLES_OUT, references, metrics, streamwise["temporal"])
    write_motion_table(paths.TABLES_OUT, motion_reference, metrics["Tangaroa"])
    print("Wrote Tables S2–S4")


if __name__ == "__main__":
    main()
