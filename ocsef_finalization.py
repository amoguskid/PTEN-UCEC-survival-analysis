"""Phase 4: validate PTEN reconstruction and create the final project outputs.

This file is deliberately executable both from the master notebook and from the
command line.  It consumes only the raw inputs and the outputs written by Phases
1--3; no state from the notebook kernel is required.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter


PROJECT_ROOT = Path(__file__).resolve().parent
RAW_DIR = PROJECT_ROOT / "data_raw"
PROCESSED_DIR = PROJECT_ROOT / "data_processed"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"

for directory in (PROCESSED_DIR, RESULTS_DIR, FIGURES_DIR):
    directory.mkdir(parents=True, exist_ok=True)


def require(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required input: {path}")
    return path


def sample_id_column(frame, source):
    candidates = ("sampleId", "SAMPLE_ID", "Sample ID", "studyID:sampleId")
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
    raise ValueError(f"Could not identify the sample-ID column in {source}.")


def normalized_sample_ids(series):
    return series.astype("string").str.strip().str.split(":").str[-1].str.upper()


def mutation_sample_ids(frame):
    sample_column = sample_id_column(frame, "the mutation table")
    gene_columns = [column for column in ("Hugo_Symbol", "GENE", "Gene") if column in frame]
    if gene_columns:
        frame = frame.loc[frame[gene_columns[0]].astype("string").str.upper().eq("PTEN")]
    elif "PTEN" in frame:
        annotation = frame["PTEN"].astype("string").str.strip().str.lower()
        frame = frame.loc[~annotation.isin(("", "0", "nan", "na", "none"))]
    else:
        raise ValueError("The mutation table has neither a gene column nor a PTEN column.")
    return set(normalized_sample_ids(frame[sample_column]).dropna())


def save_figure(fig, stem):
    fig.savefig(FIGURES_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


s1_path = require(RAW_DIR / "Supplementary_Table_S1_PTEN_altered_unaltered_sample_matrix.tsv")
s2_path = require(RAW_DIR / "Supplementary_Table_S2_PTEN_discrete_CNA_table.tsv")
s3_path = require(RAW_DIR / "Supplementary_Table_S3_PTEN_mutation_table.tsv")
survival_path = require(PROCESSED_DIR / "survival_analysis_dataset.csv")
subtype_path = require(RESULTS_DIR / "subtype_statistics.csv")
adjustment_path = require(RESULTS_DIR / "PTEN_adjustment_comparison.csv")

s1 = pd.read_csv(s1_path, sep="\t")
s2 = pd.read_csv(s2_path, sep="\t")
s3 = pd.read_csv(s3_path, sep="\t")

s1_id = sample_id_column(s1, "the PTEN status matrix")
s2_id = sample_id_column(s2, "the discrete CNA table")
if "Altered" not in s1 or "PTEN" not in s2:
    raise ValueError("Expected Altered in S1 and PTEN in S2.")

portal = pd.DataFrame({
    "SAMPLE_ID": normalized_sample_ids(s1[s1_id]),
    "PORTAL_ALTERED": pd.to_numeric(s1["Altered"], errors="raise").astype(int),
})
cna = pd.DataFrame({
    "SAMPLE_ID": normalized_sample_ids(s2[s2_id]),
    "PTEN_CNA": pd.to_numeric(s2["PTEN"], errors="coerce"),
})
if portal["SAMPLE_ID"].duplicated().any() or cna["SAMPLE_ID"].duplicated().any():
    raise ValueError("Duplicate sample identifiers prevent one-to-one reconstruction.")

reconstruction = portal.merge(cna, on="SAMPLE_ID", how="left", validate="one_to_one")
mutated = mutation_sample_ids(s3)
reconstruction["MUTATION_PRESENT"] = reconstruction["SAMPLE_ID"].isin(mutated)
reconstruction["HIGH_LEVEL_CNA_PRESENT"] = reconstruction["PTEN_CNA"].isin((-2, 2))
reconstruction["RECONSTRUCTED_ALTERED"] = (
    reconstruction["MUTATION_PRESENT"] | reconstruction["HIGH_LEVEL_CNA_PRESENT"]
).astype(int)
reconstruction["ALTERATION_SOURCE"] = np.select(
    [
        reconstruction["MUTATION_PRESENT"] & reconstruction["HIGH_LEVEL_CNA_PRESENT"],
        reconstruction["MUTATION_PRESENT"],
        reconstruction["HIGH_LEVEL_CNA_PRESENT"],
    ],
    ["Mutation and high-level CNA", "Mutation only", "High-level CNA only"],
    default="None",
)
reconstruction["MATCHES_PORTAL_LABEL"] = (
    reconstruction["RECONSTRUCTED_ALTERED"] == reconstruction["PORTAL_ALTERED"]
)
reconstruction.to_csv(PROCESSED_DIR / "PTEN_status_reconstruction.csv", index=False)

summary_values = {
    "samples_reconstructed": len(reconstruction),
    "mutation_present": int(reconstruction["MUTATION_PRESENT"].sum()),
    "high_level_cna_present": int(reconstruction["HIGH_LEVEL_CNA_PRESENT"].sum()),
    "mutation_and_high_level_cna": int((reconstruction["MUTATION_PRESENT"] & reconstruction["HIGH_LEVEL_CNA_PRESENT"]).sum()),
    "reconstructed_altered": int(reconstruction["RECONSTRUCTED_ALTERED"].sum()),
    "portal_altered": int(reconstruction["PORTAL_ALTERED"].sum()),
    "exact_matches": int(reconstruction["MATCHES_PORTAL_LABEL"].sum()),
    "mismatches": int((~reconstruction["MATCHES_PORTAL_LABEL"]).sum()),
}
pd.DataFrame(summary_values.items(), columns=["measure", "value"]).to_csv(
    RESULTS_DIR / "PTEN_status_reconstruction_summary.csv", index=False
)
validation = pd.DataFrame([
    {"check": "All 549 samples reconstructed", "passed": summary_values["samples_reconstructed"] == 549},
    {"check": "181 altered samples reconstructed", "passed": summary_values["reconstructed_altered"] == 181},
    {"check": "Reconstruction matches every portal label", "passed": summary_values["mismatches"] == 0},
])
validation.to_csv(RESULTS_DIR / "PTEN_status_reconstruction_validation.csv", index=False)
if not validation["passed"].all():
    raise RuntimeError("One or more Phase 4 PTEN reconstruction checks failed.")

# Final Kaplan--Meier figure, using the Phase 1 cohort and unchanged grouping.
survival = pd.read_csv(survival_path)
fig, ax = plt.subplots(figsize=(8, 6))
for value, label, color in ((1, "PTEN altered", "#D55E00"), (0, "PTEN unaltered", "#0072B2")):
    group = survival.loc[survival["Altered"].eq(value)]
    KaplanMeierFitter(label=label).fit(
        group["OS_MONTHS_NUMERIC"], event_observed=group["OS_EVENT"]
    ).plot_survival_function(ax=ax, ci_show=True, color=color)
ax.set(xlabel="Overall survival (months)", ylabel="Survival probability", title="Overall survival by PTEN alteration status")
ax.grid(alpha=0.2)
save_figure(fig, "final_01_kaplan_meier")

# Final subtype-frequency figure from the Phase 2 estimates.
subtypes = pd.read_csv(subtype_path)
label_column = "subtype" if "subtype" in subtypes else "SUBTYPE_READABLE"
fig, ax = plt.subplots(figsize=(9, 6))
ax.bar(subtypes[label_column], subtypes["altered_percent"], color="#56B4E9")
ax.set(ylabel="PTEN altered (%)", xlabel="Molecular subtype", title="PTEN alteration frequency by molecular subtype")
ax.tick_params(axis="x", rotation=20)
fig.tight_layout()
save_figure(fig, "final_02_subtype_frequency")

# Final forest figure from the unchanged Phase 3 model estimates.
models = pd.read_csv(adjustment_path)
y = np.arange(len(models))
fig, ax = plt.subplots(figsize=(8, 5))
ax.errorbar(models["PTEN_HR"], y, xerr=[models["PTEN_HR"] - models["CI_95_low"], models["CI_95_high"] - models["PTEN_HR"]], fmt="o", color="#0072B2", capsize=4)
ax.axvline(1, color="black", linestyle="--", linewidth=1)
ax.set(yticks=y, yticklabels=models["model"], xlabel="PTEN hazard ratio (95% CI)", title="PTEN association across Cox models")
ax.invert_yaxis()
fig.tight_layout()
save_figure(fig, "final_03_adjusted_cox_forest")

survival_stats = pd.read_csv(require(RESULTS_DIR / "survival_statistics.csv")).iloc[0]
chi_square = pd.read_csv(require(RESULTS_DIR / "subtype_chi_square_test.csv")).iloc[0]
fully_adjusted = models.loc[models["model"].eq("model_2_subtype_adjusted")].iloc[0]
findings = pd.DataFrame([
    {"analysis": "Direct PTEN reconstruction", "result": f"{summary_values['exact_matches']}/{summary_values['samples_reconstructed']} labels matched"},
    {"analysis": "Independent survival rerun", "result": f"HR {survival_stats['cox_hr_altered_vs_unaltered']:.3f}"},
    {"analysis": "Subtype association", "result": f"Chi-square {chi_square['chi_square']:.2f}"},
    {"analysis": "Fully adjusted Cox model", "result": f"HR {fully_adjusted['PTEN_HR']:.3f}"},
])
findings.to_csv(RESULTS_DIR / "final_key_findings.csv", index=False)

print("Phase 4 completed successfully.")
print(pd.DataFrame(summary_values.items(), columns=["measure", "value"]).to_string(index=False))
