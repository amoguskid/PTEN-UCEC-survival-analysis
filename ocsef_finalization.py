"""Finalize the OCSEF analysis from saved inputs and Phase 1--3 outputs.

This module deliberately does not fit statistical models. It reconstructs the
portal PTEN label, assembles the final findings table, and renders the three
final figures from the saved Phase 1–3 datasets and estimates.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter


S1_NAME = "Supplementary_Table_S1_PTEN_altered_unaltered_sample_matrix.tsv"
S2_NAME = "Supplementary_Table_S2_PTEN_discrete_CNA_table.tsv"
S3_NAME = "Supplementary_Table_S3_PTEN_mutation_table.tsv"

FINAL_OUTPUTS = (
    "data_processed/PTEN_status_reconstruction.csv",
    "results/PTEN_status_reconstruction_validation.csv",
    "results/PTEN_status_reconstruction_summary.csv",
    "results/final_key_findings.csv",
    "figures/final_01_kaplan_meier.png",
    "figures/final_01_kaplan_meier.pdf",
    "figures/final_02_subtype_frequency.png",
    "figures/final_02_subtype_frequency.pdf",
    "figures/final_03_adjusted_cox_forest.png",
    "figures/final_03_adjusted_cox_forest.pdf",
)


def _column(frame, alternatives, description):
    """Return the first supported column, retaining documented portal aliases."""
    for name in alternatives:
        if name in frame.columns:
            return name
    raise ValueError(
        f"Missing {description}; expected one of {list(alternatives)}, "
        f"found {list(frame.columns)}."
    )


def _sample_ids(frame):
    """Read sample identifiers from supplied or standard cBioPortal columns."""
    name = _column(
        frame,
        ("Sample ID", "SAMPLE_ID", "sampleId", "Tumor_Sample_Barcode"),
        "sample identifier column",
    )
    return _normalize_sample_ids(frame[name])


def _normalize_sample_ids(series):
    """Remove an optional study prefix and normalize TCGA sample identifiers."""
    return (
        series.astype("string")
        .str.split(":", n=1)
        .str[-1]
        .str.strip()
        .str.upper()
    )


def _read_table(path):
    return pd.read_csv(path, sep="\t", low_memory=False)


def _as_binary(series, name):
    values = pd.to_numeric(series, errors="coerce")
    if values.isna().any() or not set(values.unique()).issubset({0, 1}):
        raise ValueError(f"{name} must contain only 0 and 1.")
    return values.astype(int)


def reconstruct_pten(s1, s2, s3):
    """Reconstruct mutation-or-high-level-CNA status for every S1 sample."""
    # S1 commonly prefixes the sample with the study ID.
    if "studyID:sampleId" in s1.columns:
        s1_sample = _normalize_sample_ids(s1["studyID:sampleId"])
    else:
        s1_sample = _sample_ids(s1)
    portal_col = _column(s1, ("Altered", "ALTERED", "PORTAL_ALTERED"), "portal label")

    cna_sample = _sample_ids(s2)
    cna_col = _column(s2, ("PTEN", "CNA", "DISCRETE_CNA"), "PTEN CNA column")
    cna = pd.DataFrame({
        "SAMPLE_ID": cna_sample,
        "PTEN_CNA": pd.to_numeric(s2[cna_col], errors="coerce"),
    })
    if cna["SAMPLE_ID"].duplicated().any():
        raise ValueError("Duplicate sample identifiers in the discrete CNA table.")
    if set(s1_sample.dropna()) != set(cna["SAMPLE_ID"].dropna()):
        raise ValueError("S1 and S2 do not contain the same sample identifiers.")

    gene_col = _column(
        s3, ("Gene", "Hugo_Symbol", "HUGO_SYMBOL", "GENE"), "gene column"
    )
    mutation_samples = set(
        _sample_ids(s3).loc[
            s3[gene_col].astype("string").str.strip().str.upper().eq("PTEN")
        ].dropna()
    )

    output = pd.DataFrame({
        "SAMPLE_ID": s1_sample,
        "PORTAL_ALTERED": _as_binary(s1[portal_col], "Portal label"),
    }).merge(cna, on="SAMPLE_ID", how="left", validate="one_to_one")
    if output["SAMPLE_ID"].duplicated().any():
        raise ValueError("Duplicate sample identifiers in the supplied PTEN matrix.")
    # The supplied CNA table contains ten expected NP calls. They remain missing;
    # only numeric -2 and +2 calls count as high-level alterations.
    if int(output["PTEN_CNA"].isna().sum()) != 10:
        raise ValueError(
            "Expected exactly 10 unavailable (NP) PTEN CNA calls, found "
            f"{int(output['PTEN_CNA'].isna().sum())}."
        )
    output["MUTATION_PRESENT"] = output["SAMPLE_ID"].isin(mutation_samples)
    output["HIGH_LEVEL_CNA_PRESENT"] = output["PTEN_CNA"].isin((-2, 2))
    output["RECONSTRUCTED_ALTERED"] = (
        output["MUTATION_PRESENT"] | output["HIGH_LEVEL_CNA_PRESENT"]
    ).astype(int)
    output["ALTERATION_SOURCE"] = "none"
    output.loc[output["MUTATION_PRESENT"], "ALTERATION_SOURCE"] = "mutation only"
    output.loc[output["HIGH_LEVEL_CNA_PRESENT"], "ALTERATION_SOURCE"] = "CNA only"
    output.loc[
        output["MUTATION_PRESENT"] & output["HIGH_LEVEL_CNA_PRESENT"],
        "ALTERATION_SOURCE",
    ] = "mutation and CNA"
    output["MATCHES_PORTAL_LABEL"] = (
        output["RECONSTRUCTED_ALTERED"] == output["PORTAL_ALTERED"]
    )
    return output


def reconstruction_summary(reconstruction):
    measures = {
        "total_samples": len(reconstruction),
        "mutation_present": int(reconstruction["MUTATION_PRESENT"].sum()),
        "high_level_cna_present": int(reconstruction["HIGH_LEVEL_CNA_PRESENT"].sum()),
        "mutation_and_high_level_cna": int(
            (reconstruction["MUTATION_PRESENT"] & reconstruction["HIGH_LEVEL_CNA_PRESENT"]).sum()
        ),
        "reconstructed_altered": int(reconstruction["RECONSTRUCTED_ALTERED"].sum()),
        "portal_altered": int(reconstruction["PORTAL_ALTERED"].sum()),
        "exact_matches": int(reconstruction["MATCHES_PORTAL_LABEL"].sum()),
        "mismatches": int((~reconstruction["MATCHES_PORTAL_LABEL"]).sum()),
    }
    return pd.DataFrame(measures.items(), columns=["measure", "count"])


def _validate_reconstruction(summary):
    counts = summary.set_index("measure")["count"]
    expected = {
        "total_samples": 549,
        "mutation_present": 161,
        "high_level_cna_present": 26,
        "mutation_and_high_level_cna": 6,
        "reconstructed_altered": 181,
        "exact_matches": 549,
    }
    validation = pd.DataFrame(
        [
            {
                "check": key,
                "expected": expected_value,
                "observed": int(counts[key]),
                "matches_expected": int(counts[key]) == expected_value,
            }
            for key, expected_value in expected.items()
        ]
    )
    return validation


def _save_figure(figure, figures_dir, stem):
    figure.savefig(figures_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    figure.savefig(figures_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(figure)


def make_final_figures(processed_dir, results_dir, figures_dir):
    """Render the three final figures from the saved Phase 1–3 artifacts."""
    survival = pd.read_csv(processed_dir / "survival_analysis_dataset.csv")
    subtype_stats = pd.read_csv(results_dir / "subtype_statistics.csv")
    comparison = pd.read_csv(results_dir / "PTEN_adjustment_comparison.csv")

    figure, axis = plt.subplots(figsize=(8, 6))
    for value, label, color in (
        (1, "PTEN-altered", "#0072B2"),
        (0, "PTEN-unaltered", "#D55E00"),
    ):
        group = survival.loc[survival["Altered"].eq(value)]
        KaplanMeierFitter(label=f"{label} (n={len(group)})").fit(
            group["OS_MONTHS_NUMERIC"], event_observed=group["OS_EVENT"]
        ).plot_survival_function(
            ax=axis,
            ci_show=True,
            color=color,
            linewidth=2,
            show_censors=True,
        )
    axis.set(
        xlabel="Overall survival (months)",
        ylabel="Estimated survival probability",
        title="Overall survival by PTEN alteration status",
    )
    axis.grid(alpha=0.2)
    figure.tight_layout()
    _save_figure(figure, figures_dir, "final_01_kaplan_meier")

    labels = subtype_stats["subtype"]
    positions = np.arange(len(subtype_stats))
    percentages = subtype_stats["altered_percent"]
    lower = percentages - subtype_stats["wilson_ci95_lower_percent"]
    upper = subtype_stats["wilson_ci95_upper_percent"] - percentages
    figure, axis = plt.subplots(figsize=(9, 6))
    bars = axis.bar(
        positions,
        percentages,
        yerr=np.vstack([lower, upper]),
        capsize=5,
        color="#56B4E9",
    )
    axis.set_xticks(positions, labels, rotation=15, ha="right")
    axis.set(
        xlabel="TCGA-UCEC molecular subtype",
        ylabel="PTEN-altered tumors (%)",
        title="PTEN alteration frequency by molecular subtype",
    )
    axis.grid(axis="y", alpha=0.2)
    for bar, row in zip(bars, subtype_stats.itertuples()):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 1,
            f"{row.altered_n}/{row.total_n}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    figure.tight_layout()
    _save_figure(figure, figures_dir, "final_02_subtype_frequency")

    positions = np.arange(len(comparison))
    figure, axis = plt.subplots(figsize=(8, 5))
    axis.errorbar(
        comparison["PTEN_HR"],
        positions,
        xerr=np.vstack(
            [
                comparison["PTEN_HR"] - comparison["CI_95_low"],
                comparison["CI_95_high"] - comparison["PTEN_HR"],
            ]
        ),
        fmt="o",
        color="#0072B2",
        capsize=4,
    )
    axis.axvline(1, color="black", linestyle="--", linewidth=1)
    axis.set_yticks(positions, comparison["model"])
    axis.set(
        xlabel="PTEN hazard ratio (95% CI)",
        title="PTEN association across Cox models",
    )
    axis.invert_yaxis()
    figure.tight_layout()
    _save_figure(figure, figures_dir, "final_03_adjusted_cox_forest")


def build_final_key_findings(results_dir, summary):
    """Build the narrative table strictly from saved result CSV files."""
    survival = pd.read_csv(results_dir / "survival_statistics.csv").iloc[0]
    subtype = pd.read_csv(results_dir / "subtype_chi_square_test.csv").iloc[0]
    adjusted = pd.read_csv(results_dir / "PTEN_adjustment_comparison.csv")
    model = adjusted.loc[adjusted["model"].eq("model_2_subtype_adjusted")].iloc[0]
    counts = summary.set_index("measure")["count"]
    return pd.DataFrame([
        {"finding": "Direct PTEN reconstruction", "analysis_n": int(counts["total_samples"]),
         "result": f'{int(counts["exact_matches"])}/{int(counts["total_samples"])} exact matches',
         "interpretation": "Mutation or high-level discrete CNA reproduces the supplied portal label."},
        {"finding": "Independent survival rerun", "analysis_n": int(survival["total_n"]),
         "result": f'Log-rank p={survival["logrank_p_value"]:.6g}; unadjusted HR={survival["cox_hr_altered_vs_unaltered"]:.4f}',
         "interpretation": "Exported-data rerun; not an exact cBioPortal cohort replication."},
        {"finding": "PTEN status by molecular subtype", "analysis_n": int(subtype["analysis_n"]),
         "result": f'Chi-square={subtype["chi_square"]:.4f}, df={int(subtype["degrees_of_freedom"])}, p={subtype["p_value"]:.6g}, Cramer\'s V={subtype["cramers_v"]:.4f}',
         "interpretation": "Unadjusted association; it does not establish causation or prognosis."},
        {"finding": "Fully adjusted Cox result", "analysis_n": int(model["patients_n"]),
         "result": f'PTEN HR={model["PTEN_HR"]:.4f} (95% CI {model["CI_95_low"]:.4f}-{model["CI_95_high"]:.4f}); p={model["PTEN_p_value"]:.4f}',
         "interpretation": "PTEN was not independently associated with survival after adjustment."},
    ])


def main(project_root=None):
    root = Path(project_root) if project_root else Path(__file__).resolve().parent
    raw, processed, results, figures = (root / name for name in
        ("data_raw", "data_processed", "results", "figures"))
    for directory in (processed, results, figures):
        directory.mkdir(parents=True, exist_ok=True)

    # Remove only outputs owned by this phase, so a failed rerun cannot leave stale results.
    for relative in FINAL_OUTPUTS:
        (root / relative).unlink(missing_ok=True)

    paths = [raw / name for name in (S1_NAME, S2_NAME, S3_NAME)]
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required Phase 4 inputs: {missing}")

    reconstruction = reconstruct_pten(*(_read_table(path) for path in paths))
    summary = reconstruction_summary(reconstruction)
    validation = _validate_reconstruction(summary)
    reconstruction.to_csv(processed / "PTEN_status_reconstruction.csv", index=False)
    summary.to_csv(results / "PTEN_status_reconstruction_summary.csv", index=False)
    validation.to_csv(
        results / "PTEN_status_reconstruction_validation.csv", index=False
    )
    if not validation["matches_expected"].all():
        failures = validation.loc[~validation["matches_expected"]].to_dict("records")
        raise ValueError(f"PTEN reconstruction count validation failed: {failures}")
    build_final_key_findings(results, summary).to_csv(
        results / "final_key_findings.csv", index=False
    )
    make_final_figures(processed, results, figures)

    matched = int(reconstruction["MATCHES_PORTAL_LABEL"].sum())
    print(f"Phase 4 PTEN reconstruction: {matched}/{len(reconstruction)} labels matched.")
    print("Final figures and final_key_findings.csv were written successfully.")


if __name__ == "__main__":
    main()
