"""Compare generated tables and figures with the paper LaTeX directory."""

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anchormtm_paths as paths

SUPPLEMENT = (
    "table_family.tex",
    "table_temporal.tex",
    "table_tangaroa_motion.tex",
    "table_cylinder_exact.tex",
    "table_ablation.tex",
    "table_depth.tex",
    "table_motion.tex",
)


def normalized(text):
    return " ".join(text.split())


def check_tables(submission, generated):
    reference = submission / "tables" if submission else paths.RESULTS / "paper_tables"
    table1 = (reference / "table1_comparison.tex").read_text(encoding="utf-8")
    submitted = [
        normalized(line)
        for line in table1.splitlines()
        if " & Exact & " in line or " & Balanced & " in line
    ]
    actual = (generated / "table1_rows.tex").read_text(encoding="utf-8")
    if len(submitted) != 12 or submitted != [
        normalized(line) for line in actual.splitlines()
    ]:
        raise ValueError("Table 1 data rows or emphasis differ from the submission")

    table_s1 = (reference / "table2_parameters_runtime.tex").read_text(encoding="utf-8")
    labels = dict(
        N_r=r"$N_r$", s_w=r"$s_w$", t_opt=r"$t_{opt}$", t_total=r"$t_{total}$"
    )
    names = ("Ring", "Wildfires", "Storms", "PM2.5")
    with (generated / "table_parameters_runtime.csv").open(newline="") as stream:
        measured = {
            (r["sequence"], r["parameter"]): float(r["value"])
            for r in csv.DictReader(stream)
        }
    for key, label in labels.items():
        match = re.search(r"^" + re.escape(label) + r"\s*&(.+?)\\\\", table_s1, re.M)
        if not match:
            raise ValueError(f"Table S1 is missing {label}")
        values = match[1].split("&")
        if len(values) != len(names):
            raise ValueError(f"Table S1 row {label} has the wrong number of cells")
        for name, cell in zip(names, values, strict=True):
            value = cell.strip()
            digits = len(value.partition(".")[2])
            if f"{measured[name, key]:.{digits}f}" != value:
                raise ValueError(f"Table S1 differs for {name}, {key}")

    for name in SUPPLEMENT:
        expected = (reference / name).read_text(encoding="utf-8").strip()
        actual = (generated / name).read_text(encoding="utf-8").strip()
        if actual != expected:
            raise ValueError(f"{name} differs from the submission")
    print("Table 1 and Tables S1–S8 match the paper at the reported precision.")


def check_figures(submission, generated):
    import numpy as np
    from PIL import Image

    submitted = sorted((submission / "figures").glob("*.png"))
    if not submitted:
        raise ValueError("No submitted PNGs found")
    identical = 0
    for reference in submitted:
        name = {"budget_response_v30.png": "budget_response.png"}.get(
            reference.name, reference.name
        )
        path = generated / name
        with (
            Image.open(reference) as submitted_image,
            Image.open(path) as generated_image,
        ):
            if submitted_image.size != generated_image.size:
                print(
                    f"{path.name}: {generated_image.size} vs {submitted_image.size} (submitted)"
                )
                continue
            expected = np.asarray(submitted_image.convert("RGB"), dtype=np.int16)
            actual = np.asarray(generated_image.convert("RGB"), dtype=np.int16)
            difference = np.abs(actual - expected)
            if difference.any():
                changed = np.mean(np.any(difference != 0, axis=2))
                print(
                    f"{path.name}: {changed:.2%} changed pixels, mean |RGB delta| {difference.mean():.3f}/255"
                )
            else:
                identical += 1
    print(f"{identical}/{len(submitted)} submitted figures match pixel for pixel.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "submission",
        nargs="?",
        type=Path,
        help="Paper source directory for an external comparison",
    )
    parser.add_argument(
        "--figures", action="store_true", help="Also compare the exported PNGs"
    )
    args = parser.parse_args()
    check_tables(args.submission, paths.TABLES_OUT)
    if args.figures:
        if args.submission is None:
            parser.error("--figures requires the paper source directory")
        check_figures(args.submission, paths.FIGURES_OUT / "png")


if __name__ == "__main__":
    main()
