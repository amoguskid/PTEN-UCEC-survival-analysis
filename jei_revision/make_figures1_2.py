#!/usr/bin/env python3
"""Recreate manuscript Figures 1–2 from the rerun cohort and subtype counts.

Example:
    python jei_revision/make_figures1_2.py \
        --os jei_revision_outputs/data_processed/OS_eligible_cohort.csv \
        --subtype-statistics results/subtype_statistics.csv \
        --out jei_revision_outputs/figures

Figure 1 uses all recorded follow-up to fit Kaplan–Meier and unadjusted Cox
models and the log-rank test, but displays 0–120 months. Figure 2 derives
percentages and Wilson 95% intervals from the subtype counts. No survival
coordinates, scientific estimates, or patient counts are supplied by hand.
"""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.stats import chi2_contingency, norm


BLUE = "#0072B2"
ORANGE = "#D55E00"
SUBTYPE_ORDER = [
    "POLE-ultramutated", "MSI-hypermutated", "Copy-number low", "Copy-number high"
]
SUBTYPE_COLORS = ["#8EC5E8", "#5FA8D3", BLUE, "#AAB7C4"]
DISPLAY_MONTHS = 120
RISK_TIMES = np.arange(0, DISPLAY_MONTHS + 1, 24)


def read_os(path):
    frame = pd.read_csv(path)
    required = {"OS_TIME", "OS_EVENT", "Altered"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"OS input is missing columns: {sorted(missing)}")
    if frame.empty:
        raise ValueError("OS input is empty")
    for name in required:
        frame[name] = pd.to_numeric(frame[name], errors="raise")
        if not np.isfinite(frame[name]).all():
            raise ValueError(f"OS input contains missing/nonfinite {name}")
    if not frame["OS_TIME"].ge(0).all():
        raise ValueError("OS durations must be nonnegative and measured in months")
    if not frame["OS_EVENT"].isin([0, 1]).all() or not frame["Altered"].isin([0, 1]).all():
        raise ValueError("OS_EVENT and Altered must be coded 0/1")
    if set(frame["Altered"].unique()) != {0, 1}:
        raise ValueError("OS input must contain both PTEN groups")
    if "patientId" in frame and not frame["patientId"].is_unique:
        raise ValueError("OS input contains repeated patient identifiers")
    frame[["OS_EVENT", "Altered"]] = frame[["OS_EVENT", "Altered"]].astype(int)
    return frame


def save_figure(fig, out, stem, dpi):
    # Fixed canvas and explicit margins preserve the manuscript aspect ratio.
    fig.savefig(out / f"{stem}.png", dpi=dpi, facecolor="white")
    fig.savefig(out / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def make_figure1(frame, out, dpi):
    groups = [
        ("PTEN-altered", frame.loc[frame["Altered"].eq(1)], BLUE),
        ("PTEN-unaltered", frame.loc[frame["Altered"].eq(0)], ORANGE),
    ]
    logrank = logrank_test(
        groups[0][1]["OS_TIME"], groups[1][1]["OS_TIME"],
        event_observed_A=groups[0][1]["OS_EVENT"],
        event_observed_B=groups[1][1]["OS_EVENT"],
    )
    cox = CoxPHFitter()
    cox.fit(frame[["OS_TIME", "OS_EVENT", "Altered"]],
            duration_col="OS_TIME", event_col="OS_EVENT")
    coefficient = cox.summary.loc["Altered"]
    hr = float(coefficient["exp(coef)"])
    low = float(coefficient["exp(coef) lower 95%"])
    high = float(coefficient["exp(coef) upper 95%"])

    # Original manuscript image: 1420 × 957 pixels at approximately 220 dpi.
    fig = plt.figure(figsize=(1420 / 220, 957 / 220), facecolor="white")
    ax = fig.add_axes([0.144, 0.216, 0.829, 0.700])
    curve_tables, risk_rows, group_summary = [], [], []
    for name, group, color in groups:
        km = KaplanMeierFitter(label=f"{name} (n={len(group)})")
        km.fit(group["OS_TIME"], event_observed=group["OS_EVENT"])
        km.plot_survival_function(
            ax=ax, ci_show=True, ci_alpha=0.12, color=color, linewidth=1.55,
            show_censors=True, censor_styles={"ms": 2.0, "marker": "+"},
        )
        curve = km.survival_function_.copy()
        curve.columns = ["survival_probability"]
        curve = curve.join(km.confidence_interval_)
        curve.columns = ["survival_probability", "CI_95_lower", "CI_95_upper"]
        curve = curve.join(km.event_table[["at_risk", "observed", "censored"]])
        curve.index.name = "time_months"
        curve = curve.reset_index()
        curve.insert(0, "group", name)
        curve_tables.append(curve)
        for time in RISK_TIMES:
            risk_rows.append({
                "group": name, "time_months": int(time),
                "at_risk_at_start_of_time": int(group["OS_TIME"].ge(time).sum()),
            })
        median = float(km.median_survival_time_)
        group_summary.append({
            "group": name, "patients_n": int(len(group)),
            "deaths_n": int(group["OS_EVENT"].sum()),
            "censored_n": int(group["OS_EVENT"].eq(0).sum()),
            "maximum_observed_duration_months": float(group["OS_TIME"].max()),
            "KM_median_months": median if np.isfinite(median) else None,
        })

    ax.set(xlim=(0, DISPLAY_MONTHS), ylim=(0, 1.02), xticks=RISK_TIMES,
           yticks=np.arange(0, 1.01, 0.2), xlabel="Overall survival (months)",
           ylabel="Estimated survival probability")
    ax.xaxis.label.set_size(9)
    ax.yaxis.label.set_size(9)
    ax.tick_params(axis="both", labelsize=8, length=3, width=0.6, pad=3)
    ax.set_axisbelow(True)
    ax.grid(color="#D9E5F2", linewidth=0.65)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.6)
    ax.legend(loc="lower left", frameon=False, fontsize=8.2,
              handlelength=1.9, borderpad=0.1, labelspacing=0.4)
    annotation = (
        f"Log-rank p = {logrank.p_value: .3f}\n"
        f"Unadjusted HR = {hr:.2f} (95% CI {low:.2f}–{high:.2f})\n"
        f"Deaths = {int(frame['OS_EVENT'].sum())}"
    )
    ax.text(0.985, 0.960, annotation, transform=ax.transAxes,
            ha="right", va="top", fontsize=7.5,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white",
                  "edgecolor": "#D9E2EC", "linewidth": 0.8})

    ax.text(-0.090, -0.175, "At risk", transform=ax.transAxes,
            ha="left", va="center", fontsize=7.5, fontweight="bold")
    for index, (name, group, color) in enumerate(groups):
        y = -0.225 - index * 0.060
        ax.text(-0.026, y, name, transform=ax.transAxes,
                ha="right", va="center", fontsize=7, color=color)
        for time in RISK_TIMES:
            count = int(group["OS_TIME"].ge(time).sum())
            ax.text(time / DISPLAY_MONTHS, y, str(count), transform=ax.transAxes,
                    ha="center", va="center", fontsize=7)

    save_figure(fig, out, "Figure1_Kaplan_Meier", dpi)
    pd.concat(curve_tables, ignore_index=True).to_csv(out / "Figure1_plotted_survival.csv", index=False)
    pd.DataFrame(risk_rows).to_csv(out / "Figure1_numbers_at_risk.csv", index=False)
    return {
        "groups": group_summary,
        "logrank_statistic": float(logrank.test_statistic),
        "logrank_p_value": float(logrank.p_value),
        "unadjusted_Cox_HR": hr, "CI_95_lower": low, "CI_95_upper": high,
        "unadjusted_Cox_p_value": float(coefficient["p"]),
        "display_window_months": [0, DISPLAY_MONTHS],
        "all_recorded_followup_used": True, "extrapolation_used": False,
        "risk_count_definition": "Observed duration >= displayed time (start-of-time risk set)",
    }


def subtype_statistics(path):
    table = pd.read_csv(path)
    required = {"subtype", "altered_n", "total_n"}
    if not required.issubset(table.columns):
        raise ValueError(f"Subtype input must contain {sorted(required)}")
    aliases = {"Copy-number-low": "Copy-number low", "Copy-number-high": "Copy-number high"}
    table["subtype"] = table["subtype"].replace(aliases)
    if not table["subtype"].is_unique or set(table["subtype"]) != set(SUBTYPE_ORDER):
        raise ValueError("Subtype input must contain the four unique TCGA UCEC molecular subtypes")
    table = table.set_index("subtype").loc[SUBTYPE_ORDER].reset_index()
    for column in ("altered_n", "total_n"):
        values = pd.to_numeric(table[column], errors="raise")
        if not np.isfinite(values).all() or not values.eq(np.floor(values)).all():
            raise ValueError(f"Subtype {column} must contain finite integer counts")
        table[column] = values.astype(int)
    if not table["total_n"].gt(0).all() or not table["altered_n"].between(0, table["total_n"]).all():
        raise ValueError("Subtype counts must satisfy 0 <= altered_n <= total_n and total_n > 0")
    n = table["total_n"].to_numpy(dtype=float)
    proportion = table["altered_n"].to_numpy(dtype=float) / n
    z = float(norm.ppf(0.975))
    denominator = 1 + z * z / n
    center = (proportion + z * z / (2 * n)) / denominator
    half_width = z * np.sqrt(proportion * (1 - proportion) / n + z * z / (4 * n * n)) / denominator
    table["unaltered_n"] = table["total_n"] - table["altered_n"]
    table["altered_percent"] = 100 * proportion
    table["wilson_ci95_lower_percent"] = 100 * (center - half_width)
    table["wilson_ci95_upper_percent"] = 100 * (center + half_width)
    return table


def make_figure2(table, out, dpi):
    # Original manuscript image: 1418 × 717 pixels at approximately 220 dpi.
    fig = plt.figure(figsize=(1418 / 220, 717 / 220), facecolor="white")
    ax = fig.add_axes([0.2073, 0.1227, 0.5973, 0.7797])
    positions = np.arange(len(table))
    percentages = table["altered_percent"].to_numpy()
    lower = table["wilson_ci95_lower_percent"].to_numpy()
    upper = table["wilson_ci95_upper_percent"].to_numpy()
    ax.barh(positions, percentages, height=0.62, color=SUBTYPE_COLORS,
            edgecolor="white", linewidth=0.4, zorder=2)
    ax.errorbar(percentages, positions, xerr=np.vstack([percentages - lower, upper - percentages]),
                fmt="none", ecolor="#243447", elinewidth=0.85, capsize=2.7,
                capthick=0.7, zorder=3)
    ax.set_yticks(positions, table["subtype"])
    # Keep the manuscript's 0–75% window while leaving room for labels when
    # the same script is used to verify a changed input file.
    xmax = max(75, 10 * np.ceil((float(upper.max()) + 10) / 10))
    ax.set(xlim=(0, xmax), ylim=(len(table) - 0.5, -0.5),
           xlabel="Tumors with a PTEN alteration", xticks=np.arange(0, xmax, 10))
    ax.xaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    ax.xaxis.label.set_size(7.2)
    ax.xaxis.label.set_verticalalignment("center")
    ax.xaxis.set_label_coords(0.5, -0.088)
    ax.tick_params(axis="x", labelsize=6, length=3, width=0.6, pad=2.8)
    ax.tick_params(axis="y", labelsize=6.8, length=0, pad=1.2)
    ax.set_axisbelow(True)
    ax.grid(axis="x", color="#E0E7EF", linewidth=0.65)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.6)
    for y, row in enumerate(table.itertuples(index=False)):
        ax.text(row.wilson_ci95_upper_percent + 1.0, y,
                f"{row.altered_percent:.1f}%  ({row.altered_n}/{row.total_n})",
                ha="left", va="center", fontsize=6, color="#243447")
    save_figure(fig, out, "Figure2_Subtype_frequency", dpi)
    table.to_csv(out / "Figure2_plotted_subtype_statistics.csv", index=False)
    test = chi2_contingency(table[["altered_n", "unaltered_n"]].to_numpy(), correction=False)
    return {
        "patients_n": int(table["total_n"].sum()),
        "altered_n": int(table["altered_n"].sum()),
        "chi_square_statistic": float(test.statistic), "degrees_of_freedom": int(test.dof),
        "chi_square_p_value": float(test.pvalue), "confidence_interval": "Wilson 95%",
        "counts_and_percentages": table.to_dict(orient="records"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--os", type=Path, required=True,
                        help="Rerun OS_eligible_cohort.csv (OS_TIME in months)")
    parser.add_argument("--subtype-statistics", "--subtype", dest="subtype", type=Path, required=True,
                        help="Original results/subtype_statistics.csv with subtype, altered_n, total_n")
    parser.add_argument("--out", type=Path, required=True, help="Output figures directory")
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution (default: 300; PDF is vector)")
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    args.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    summary = {
        "inputs": [{"filename": path.name,
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                   for path in [args.os, args.subtype]],
        "Figure1": make_figure1(read_os(args.os), args.out, args.dpi),
        "Figure2": make_figure2(subtype_statistics(args.subtype), args.out, args.dpi),
        "png_dpi": args.dpi,
    }
    (args.out / "Figures1_2_analysis_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"out": str(args.out), "Figure1": summary["Figure1"],
                      "Figure2_patients_n": summary["Figure2"]["patients_n"]}, indent=2))


if __name__ == "__main__":
    main()
