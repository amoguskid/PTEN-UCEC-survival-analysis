#!/usr/bin/env python3
"""Plot existing Cox coefficients for JEI Figure 3; no models are refitted.

Usage:
    python make_figure3.py --coefficients results/all_cox_coefficients.csv --out figures
The saved CSV retains the exact fitted values used for the plotted intervals.
"""
from pathlib import Path
import argparse
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
from matplotlib.lines import Line2D


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coefficients", type=Path, default=Path(__file__).parent / "results/all_cox_coefficients.csv")
    parser.add_argument("--out", type=Path, default=Path(__file__).parent / "figures")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    with args.coefficients.open(newline="") as f:
        coefficients = list(csv.DictReader(f))

    specs = [
        ("PTEN only: fixed complete cases", "Altered", "PTEN status\nAltered vs unaltered", "0.60 (0.35–1.00)", "0.050"),
        ("Age, stage and subtype: fixed complete cases", "Altered", "PTEN status\nAltered vs unaltered", "1.07 (0.59–1.95)", "0.821"),
        ("Age, stage and subtype: fixed complete cases", "AGE_10Y", "Age\nPer 10-year increase", "1.10 (0.86–1.42)", "0.451"),
        ("Age, stage and subtype: fixed complete cases", "ADVANCED_STAGE", "Stage\nAdvanced vs early", "3.82 (2.34–6.24)", r"$9.0\times10^{-8}$"),
        ("Age, stage and subtype: fixed complete cases", "SUBTYPE_UCEC_CN_LOW", "Copy-number-low\nvs copy-number-high", "0.42 (0.20–0.89)", "0.023"),
        ("Age, stage and subtype: fixed complete cases", "SUBTYPE_UCEC_MSI", "MSI-hypermutated\nvs copy-number-high", "0.58 (0.30–1.11)", "0.101"),
        ("Age, stage and subtype: fixed complete cases", "SUBTYPE_UCEC_POLE", "POLE-ultramutated\nvs copy-number-high", "0.10 (0.02–0.46)", "0.003"),
    ]
    selected = []
    for model, variable, label, rounded, p_label in specs:
        match = [r for r in coefficients if r["endpoint"] == "OS" and r["definition"] == "Primary portal definition" and r["model"] == model and r["variable"] == variable]
        if len(match) != 1:
            raise ValueError(f"Expected one saved coefficient for {model}: {variable}; got {len(match)}")
        row = match[0]
        selected.append(dict(row, display_label=label, display_hr_ci=rounded, display_p=p_label))

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9.5,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.labelsize": 9.5, "xtick.labelsize": 9,
    })
    fig = plt.figure(figsize=(6.5, 4.8), facecolor="white")
    # Shared y coordinates keep the plot, labels, and numerical columns aligned.
    bottom, height = 0.145, 0.75
    forest = fig.add_axes([0.36, bottom, 0.30, height])
    forest.set_xscale("log")
    forest.set_xlim(0.02, 8)
    forest.set_ylim(0, 9.8)
    forest.set_yticks([])
    forest.xaxis.set_major_locator(FixedLocator([0.02, 0.1, 0.5, 1, 5]))
    forest.xaxis.set_major_formatter(FixedFormatter(["0.02", "0.1", "0.5", "1", "5"]))
    forest.xaxis.set_minor_locator(NullLocator())
    forest.grid(axis="x", color="#d9e5f2", linewidth=0.65, zorder=0)
    forest.axvline(1, color="#6c7b8e", linestyle="--", linewidth=1.0, zorder=1)
    for side in ("top", "left", "right"):
        forest.spines[side].set_visible(False)
    forest.spines["bottom"].set_linewidth(0.7)
    forest.tick_params(axis="x", length=3, width=0.6, pad=3)
    forest.set_xlabel("Hazard ratio (log scale)", labelpad=6)

    labels = fig.add_axes([0.025, bottom, 0.32, height], frameon=False)
    numbers = fig.add_axes([0.69, bottom, 0.29, height], frameon=False)
    for ax in (labels, numbers):
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 9.8)
        ax.axis("off")

    fig.text(0.025, 0.97, "Overall survival", fontsize=11, fontweight="bold", va="top")
    fig.text(0.025, 0.925, "459 patients; 77 observed deaths", fontsize=9.5, va="top")
    labels.text(0, 9.55, "Variable / comparison", fontsize=9, fontweight="bold", va="center")
    numbers.text(0, 9.55, "HR (95% CI)", fontsize=9, fontweight="bold", va="center")
    numbers.text(1, 9.55, "p value", fontsize=9, fontweight="bold", ha="right", va="center")
    labels.text(0, 8.92, "PTEN-only model", fontsize=9.5, fontweight="bold", va="center")
    labels.text(0, 7.20, "Fully adjusted model", fontsize=9.5, fontweight="bold", va="center")
    y_positions = [8.20, 6.45, 5.43, 4.41, 3.39, 2.37, 1.35]
    blue = "#0082be"
    for row, y in zip(selected, y_positions):
        hr = float(row["exp(coef)"])
        lo = float(row["exp(coef) lower 95%"])
        hi = float(row["exp(coef) upper 95%"])
        forest.errorbar(hr, y, xerr=[[hr-lo], [hi-hr]], fmt="o", markersize=4.5,
                        color=blue, ecolor=blue, elinewidth=1.4, capsize=3, capthick=1.1,
                        zorder=3)
        labels.text(0, y, row["display_label"], fontsize=9.5, va="center", linespacing=1.25)
        numbers.text(0, y, row["display_hr_ci"], fontsize=9, va="center")
        numbers.text(1, y, row["display_p"], fontsize=9, va="center", ha="right")

    # A subtle separator differentiates the one-variable comparison from Model 2.
    separator_y = bottom + height * 7.65 / 9.8
    fig.add_artist(Line2D([0.025, 0.98], [separator_y, separator_y],
                          transform=fig.transFigure, color="#c8c8c8", linewidth=0.6))
    fig.text(0.025, 0.010, "Early stage: I–II. Advanced stage: III–IV. MSI: microsatellite instability.",
             fontsize=8.5, va="bottom")
    fig.savefig(args.out / "Figure3_Cox_models.png", dpi=400, facecolor="white")
    fig.savefig(args.out / "Figure3_Cox_models.pdf", facecolor="white")
    plt.close(fig)
    with (args.out / "Figure3_plotted_coefficients.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(selected[0]))
        writer.writeheader()
        writer.writerows(selected)


if __name__ == "__main__":
    main()
