"""Phase 4 reconstruction checks and publication-ready summary figures.

This module deliberately consumes the artifacts created by Phases 1--3; it does
not refit or redefine any of their statistical models.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data_raw"
PROCESSED = ROOT / "data_processed"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
for directory in (PROCESSED, RESULTS, FIGURES):
    directory.mkdir(exist_ok=True)

S1 = RAW / "Supplementary_Table_S1_PTEN_altered_unaltered_sample_matrix.tsv"
S2 = RAW / "Supplementary_Table_S2_PTEN_discrete_CNA_table.tsv"
S3 = RAW / "Supplementary_Table_S3_PTEN_mutation_table.tsv"
for path in (S1, S2, S3):
    if not path.is_file():
        raise FileNotFoundError(f"Missing required Phase 4 input: {path}")


def _sample_ids(frame: pd.DataFrame) -> pd.Series:
    """Return normalized sample IDs from either documented export layout."""
    for column in ("studyID:sampleId", "SAMPLE_ID", "Tumor_Sample_Barcode"):
        if column in frame:
            return frame[column].astype("string").str.split(":", n=1).str[-1].str.strip()
    raise ValueError("Input has no sample identifier column.")


labels = pd.read_csv(S1, sep="\t", dtype={"studyID:sampleId": "string"})
cna = pd.read_csv(S2, sep="\t", dtype={"studyID:sampleId": "string"})
mutations = pd.read_csv(S3, sep="\t", dtype={"studyID:sampleId": "string"})

required_label_columns = {"studyID:sampleId", "Altered", "PTEN"}
missing = required_label_columns - set(labels.columns)
if missing:
    raise ValueError(f"S1 is missing required columns: {sorted(missing)}")
if "PTEN" not in cna:
    raise ValueError("S2 is missing required column: PTEN")

labels = labels.assign(SAMPLE_ID=_sample_ids(labels))
cna = cna.assign(SAMPLE_ID=_sample_ids(cna))
if labels["SAMPLE_ID"].duplicated().any() or cna["SAMPLE_ID"].duplicated().any():
    raise ValueError("S1 and S2 must contain one row per sample.")

portal_altered = pd.to_numeric(labels["Altered"], errors="raise").astype(int)
pten_cna = pd.to_numeric(cna["PTEN"], errors="coerce")
cna_lookup = pd.Series(pten_cna.array, index=cna["SAMPLE_ID"])

mutation_sample_ids = _sample_ids(mutations)
if "PTEN" in mutations:
    values = mutations["PTEN"].astype("string").str.strip()
    present = values.notna() & ~values.str.upper().isin({"", "0", "NA", "NAN"})
    mutation_ids = set(mutation_sample_ids[present].dropna())
elif "Hugo_Symbol" in mutations:
    mutation_ids = set(
        mutation_sample_ids[mutations["Hugo_Symbol"].astype(str).str.upper().eq("PTEN")]
        .dropna()
    )
else:
    # S3 may already be a PTEN-only long-form mutation export.
    mutation_ids = set(mutation_sample_ids.dropna())

reconstruction = pd.DataFrame({"SAMPLE_ID": labels["SAMPLE_ID"]})
reconstruction["PORTAL_ALTERED"] = portal_altered.to_numpy()
reconstruction["PTEN_CNA"] = reconstruction["SAMPLE_ID"].map(cna_lookup)
reconstruction["MUTATION_PRESENT"] = reconstruction["SAMPLE_ID"].isin(mutation_ids)
reconstruction["HIGH_LEVEL_CNA_PRESENT"] = reconstruction["PTEN_CNA"].isin([-2, 2])
reconstruction["RECONSTRUCTED_ALTERED"] = (
    reconstruction["MUTATION_PRESENT"] | reconstruction["HIGH_LEVEL_CNA_PRESENT"]
).astype(int)
reconstruction["ALTERATION_SOURCE"] = np.select(
    [
        reconstruction["MUTATION_PRESENT"] & reconstruction["HIGH_LEVEL_CNA_PRESENT"],
        reconstruction["MUTATION_PRESENT"],
        reconstruction["HIGH_LEVEL_CNA_PRESENT"],
    ],
    ["Both", "Mutation only", "CNA only"],
    default="None",
)
reconstruction["EXACT_PORTAL_LABEL_MATCH"] = (
    reconstruction["RECONSTRUCTED_ALTERED"] == reconstruction["PORTAL_ALTERED"]
)

observed = {
    "Samples reconstructed": len(reconstruction),
    "Mutation-positive samples": int(reconstruction["MUTATION_PRESENT"].sum()),
    "High-level-CNA samples": int(reconstruction["HIGH_LEVEL_CNA_PRESENT"].sum()),
    "Mutation and high-level-CNA samples": int(
        (reconstruction["MUTATION_PRESENT"] & reconstruction["HIGH_LEVEL_CNA_PRESENT"]).sum()
    ),
    "Reconstructed altered samples": int(reconstruction["RECONSTRUCTED_ALTERED"].sum()),
    "Exact portal-label matches": int(reconstruction["EXACT_PORTAL_LABEL_MATCH"].sum()),
}
expected = {
    "Samples reconstructed": 549,
    "Mutation-positive samples": 161,
    "High-level-CNA samples": 26,
    "Mutation and high-level-CNA samples": 6,
    "Reconstructed altered samples": 181,
    "Exact portal-label matches": 549,
}
validation = pd.DataFrame(
    {
        "check": expected.keys(),
        "expected": expected.values(),
        "observed": [observed[key] for key in expected],
    }
)
validation["matches_expected"] = validation["expected"] == validation["observed"]
reconstruction.to_csv(PROCESSED / "PTEN_status_reconstruction.csv", index=False)
validation.to_csv(RESULTS / "PTEN_status_reconstruction_validation.csv", index=False)
validation[["check", "observed"]].rename(columns={"observed": "count"}).to_csv(
    RESULTS / "PTEN_status_reconstruction_summary.csv", index=False
)
if not validation["matches_expected"].all():
    failures = validation.loc[~validation["matches_expected"]].to_dict("records")
    raise ValueError(f"Phase 4 reconstruction count validation failed: {failures}")

# Preserve the established analyses: final figures are copies rendered from their
# saved Phase 1--3 data/results rather than newly specified models or cohorts.
survival = pd.read_csv(PROCESSED / "survival_analysis_dataset.csv")
subtypes = pd.read_csv(PROCESSED / "subtype_analysis_dataset.csv")
comparison = pd.read_csv(RESULTS / "PTEN_adjustment_comparison.csv")

from lifelines import KaplanMeierFitter

fig, ax = plt.subplots(figsize=(8, 6))
for value, label in ((1, "PTEN-altered"), (0, "PTEN-unaltered")):
    group = survival.loc[survival["Altered"] == value]
    KaplanMeierFitter(label=f"{label} (n={len(group)})").fit(
        group["OS_MONTHS_NUMERIC"], group["OS_EVENT"]
    ).plot_survival_function(ax=ax, ci_show=True)
ax.set(xlabel="Overall survival (months)", ylabel="Estimated survival probability")
fig.tight_layout()
fig.savefig(FIGURES / "final_01_kaplan_meier.png", dpi=300)
fig.savefig(FIGURES / "final_01_kaplan_meier.pdf")
plt.close(fig)

status_column = "Altered" if "Altered" in subtypes else "PTEN_ALTERED"
subtype_column = "SUBTYPE_LABEL" if "SUBTYPE_LABEL" in subtypes else "SUBTYPE"
frequency = subtypes.groupby(subtype_column)[status_column].agg(["sum", "count"])
fig, ax = plt.subplots(figsize=(8, 5))
(100 * frequency["sum"] / frequency["count"]).plot.bar(ax=ax)
ax.set(xlabel="Molecular subtype", ylabel="PTEN altered (%)")
fig.tight_layout()
fig.savefig(FIGURES / "final_02_subtype_frequency.png", dpi=300)
fig.savefig(FIGURES / "final_02_subtype_frequency.pdf")
plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 4.5))
y = np.arange(len(comparison))
ax.errorbar(
    comparison["PTEN_HR"], y,
    xerr=[comparison["PTEN_HR"] - comparison["CI_95_low"], comparison["CI_95_high"] - comparison["PTEN_HR"]],
    fmt="o", capsize=3,
)
ax.axvline(1, color="black", linestyle="--", linewidth=1)
ax.set(yticks=y, yticklabels=comparison["model"], xlabel="PTEN hazard ratio (95% CI)")
fig.tight_layout()
fig.savefig(FIGURES / "final_03_adjusted_cox_forest.png", dpi=300)
fig.savefig(FIGURES / "final_03_adjusted_cox_forest.pdf")
plt.close(fig)

print("Phase 4 direct PTEN reconstruction: 549/549 exact label matches.")
print("Final figures created successfully from the established Phase 1-3 outputs.")
