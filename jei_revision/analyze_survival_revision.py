#!/usr/bin/env python3
"""Reproduce published OS models and add review-requested PFS/PFI models.

Run with the bundled primary Python. Optional --raw, --cdr, and --out flags
allow rerunning on the archived sources. No penalizer or multiple imputation.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.metadata as metadata
import json
import platform
import sys
import warnings

SCRIPT_DIR = Path(__file__).resolve().parent
if (SCRIPT_DIR / "dependencies").is_dir():
    sys.path.insert(0, str(SCRIPT_DIR / "dependencies"))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test

parser = argparse.ArgumentParser()
parser.add_argument("--raw", type=Path, default=SCRIPT_DIR / "data_raw")
parser.add_argument("--cdr", type=Path, default=SCRIPT_DIR / "reference_inputs/TCGA-CDR-SupplementalTableS1.xlsx")
parser.add_argument("--mutations", type=Path, default=SCRIPT_DIR / "reference_inputs/Supplementary_Table_S3_PTEN_mutation_table.tsv")
parser.add_argument("--out", type=Path, default=SCRIPT_DIR)
args = parser.parse_args()
OUT = args.out
for name in ["results", "data_processed", "figures"]:
    (OUT / name).mkdir(parents=True, exist_ok=True)

VERSIONS = {"Python": platform.python_version()}
for name in ["pandas", "numpy", "scipy", "matplotlib", "lifelines", "autograd", "formulaic", "openpyxl"]:
    VERSIONS[name] = metadata.version(name)

S1 = args.raw / "Supplementary_Table_S1_PTEN_altered_unaltered_sample_matrix.tsv"
S2 = args.raw / "Supplementary_Table_S2_PTEN_discrete_CNA_table.tsv"
S4 = args.raw / "Supplementary_Table_S4_cBioPortal_clinical_survival_table.tsv"
PATIENT = args.raw / "data_clinical_patient.txt"
hash_rows = [{"name": p.name, "path": str(p), "size_bytes": p.stat().st_size,
              "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
             for p in [S1, S2, S4, PATIENT, args.mutations, args.cdr]]
pd.DataFrame(hash_rows).to_csv(OUT / "results/input_checksums.csv", index=False)
pd.DataFrame(list(VERSIONS.items()), columns=["software", "version"]).to_csv(
    OUT / "results/software_versions.csv", index=False)

s1 = pd.read_csv(S1, sep="\t")
s2 = pd.read_csv(S2, sep="\t")
s4 = pd.read_csv(S4, sep="\t")
patient = pd.read_csv(PATIENT, sep="\t", comment="#")
cdr = pd.read_excel(args.cdr, sheet_name="TCGA-CDR")
extra = pd.read_excel(args.cdr, sheet_name="ExtraEndpoints")
mutations = pd.read_csv(args.mutations, sep="\t")
s1["sampleId"] = s1["studyID:sampleId"].str.split(":", n=1).str[-1]
assert len(s1) == 549 and s1["sampleId"].is_unique and int(s1["Altered"].sum()) == 181
assert s4["sampleId"].is_unique and s4["patientId"].is_unique
assert patient["PATIENT_ID"].is_unique
assert cdr["bcr_patient_barcode"].is_unique and extra["bcr_patient_barcode"].is_unique

# Exact original inclusion definition: missing, unrecognized, or negative OS
# values are excluded. Zero times are not excluded.
joined = s4.merge(s1[["sampleId", "Altered"]], on="sampleId", validate="one_to_one")
joined["OS_TIME"] = pd.to_numeric(joined["OS_MONTHS"], errors="coerce")
joined["OS_EVENT"] = joined["OS_STATUS"].astype("string").str.split(":").str[0].map({"0": 0, "1": 1})
joined["original_os_exclusion"] = ""
joined.loc[joined["OS_TIME"].isna(), "original_os_exclusion"] = "Missing or nonnumeric OS time"
joined.loc[joined["OS_EVENT"].isna(), "original_os_exclusion"] = "Missing or unrecognized OS status"
joined.loc[joined["OS_TIME"].lt(0), "original_os_exclusion"] = "Negative OS time"
os = joined[joined["original_os_exclusion"].eq("")].copy()
os["OS_EVENT"] = os["OS_EVENT"].astype(int)
os["Altered"] = os["Altered"].astype(int)
assert (len(os), int(os["OS_EVENT"].sum()), int(os["Altered"].sum())) == (498, 85, 146)
os = os.merge(patient[["PATIENT_ID", "SUBTYPE"]], left_on="patientId", right_on="PATIENT_ID",
              how="left", validate="one_to_one")
os["AGE_10Y"] = (pd.to_numeric(os["AGE"], errors="coerce") - 64) / 10
stage_text = os["CLINICAL_STAGE"].astype("string").str.strip()
stage = stage_text.str.extract(r"^Stage\s+(IV|III|II|I)", expand=False)
os["ADVANCED_STAGE"] = stage.map({"I": 0, "II": 0, "III": 1, "IV": 1})
SUBTYPES = ["UCEC_CN_HIGH", "UCEC_CN_LOW", "UCEC_MSI", "UCEC_POLE"]
os["SUBTYPE"] = os["SUBTYPE"].where(os["SUBTYPE"].isin(SUBTYPES))
os["complete_case"] = os[["AGE_10Y", "ADVANCED_STAGE", "SUBTYPE"]].notna().all(axis=1)
assert (int(os["complete_case"].sum()), int(os.loc[os["complete_case"], "OS_EVENT"].sum()),
        int(os.loc[os["complete_case"], "Altered"].sum())) == (459, 77, 139)

# Reconstruct both PTEN definitions from raw mutation and CNA records.
# The other two GISTIC+2 tumors also contain PTEN mutations and remain altered.
mutation_samples = set(mutations.loc[mutations["Gene"].eq("PTEN"), "Sample ID"])
assert len(mutations) == 255 and len(mutation_samples) == 161
reconstruction = s1[["sampleId", "Altered"]].merge(
    s2[["SAMPLE_ID", "PTEN"]], left_on="sampleId", right_on="SAMPLE_ID", validate="one_to_one")
reconstruction["PTEN_CNA"] = pd.to_numeric(reconstruction["PTEN"], errors="coerce")
reconstruction["MUTATION_PRESENT"] = reconstruction["sampleId"].isin(mutation_samples)
reconstruction["RECONSTRUCTED_ALTERED"] = (reconstruction["MUTATION_PRESENT"] |
    reconstruction["PTEN_CNA"].isin([-2, 2])).astype(int)
reconstruction["PTEN_NO_AMP"] = (reconstruction["MUTATION_PRESENT"] |
    reconstruction["PTEN_CNA"].eq(-2)).astype(int)
assert reconstruction["RECONSTRUCTED_ALTERED"].eq(reconstruction["Altered"]).all()
assert int(reconstruction["PTEN_NO_AMP"].sum()) == 178
AMPLIFICATION_ONLY = sorted(reconstruction.loc[
    reconstruction["RECONSTRUCTED_ALTERED"].ne(reconstruction["PTEN_NO_AMP"]), "sampleId"].tolist())
assert AMPLIFICATION_ONLY == ["TCGA-AJ-A2QK-01", "TCGA-AP-A05J-01", "TCGA-EO-A3B1-01"]
reconstruction.to_csv(OUT / "data_processed/PTEN_status_reconstruction.csv", index=False)
os = os.merge(reconstruction[["sampleId", "PTEN_NO_AMP"]], on="sampleId", how="left", validate="one_to_one")
assert int(os["PTEN_NO_AMP"].sum()) == 143
assert int(os.loc[os["complete_case"], "PTEN_NO_AMP"].sum()) == 136

os = os.merge(cdr[["bcr_patient_barcode", "PFI", "PFI.time"]], left_on="patientId",
              right_on="bcr_patient_barcode", how="left", validate="one_to_one")
os = os.rename(columns={"PFI": "PFI_EVENT", "PFI.time": "PFI_DAYS",
                        "bcr_patient_barcode": "PFI_MATCH_ID"})
os = os.merge(extra[["bcr_patient_barcode", "PFS", "PFS.time"]], left_on="patientId",
              right_on="bcr_patient_barcode", how="left", validate="one_to_one")
os = os.rename(columns={"PFS": "PFS_EVENT", "PFS.time": "PFS_DAYS",
                        "bcr_patient_barcode": "PFS_MATCH_ID"})
for endpoint in ["PFS", "PFI"]:
    os[f"{endpoint}_TIME"] = pd.to_numeric(os[f"{endpoint}_DAYS"], errors="coerce") / (365.25 / 12)
    os[f"{endpoint}_EVENT"] = pd.to_numeric(os[f"{endpoint}_EVENT"], errors="coerce")
    os[f"{endpoint}_exclusion"] = ""
    os.loc[os[f"{endpoint}_MATCH_ID"].isna(), f"{endpoint}_exclusion"] = "No exact TCGA-CDR patient match"
    os.loc[os[f"{endpoint}_TIME"].isna(), f"{endpoint}_exclusion"] = "Missing endpoint time"
    os.loc[~os[f"{endpoint}_EVENT"].isin([0, 1]), f"{endpoint}_exclusion"] = "Missing or invalid event indicator"
    os.loc[os[f"{endpoint}_TIME"].lt(0), f"{endpoint}_exclusion"] = "Negative endpoint time"
os.to_csv(OUT / "data_processed/all_498_patients_endpoint_and_covariate_audit.csv", index=False)
joined.to_csv(OUT / "data_processed/original_os_inclusion_audit.csv", index=False)
all_original = s1[["sampleId", "Altered"]].merge(
    joined[["sampleId", "patientId", "OS_TIME", "OS_EVENT", "original_os_exclusion"]],
    on="sampleId", how="left", validate="one_to_one")
all_original.loc[all_original["patientId"].isna(), "original_os_exclusion"] = "No matching clinical export"
all_original["patientId"] = all_original["patientId"].fillna(all_original["sampleId"].str[:12])
all_original.to_csv(OUT / "data_processed/all_549_original_os_cohort_audit.csv", index=False)

# Document the cBioPortal field-label discrepancy by comparison with the CDR.
cb = patient[["PATIENT_ID", "PFS_STATUS", "PFS_MONTHS"]].merge(
    cdr[["bcr_patient_barcode", "PFI", "PFI.time"]], left_on="PATIENT_ID", right_on="bcr_patient_barcode",
    validate="one_to_one").merge(extra[["bcr_patient_barcode", "PFS", "PFS.time"]],
                                on="bcr_patient_barcode", validate="one_to_one")
cb["CBIO_EVENT"] = cb["PFS_STATUS"].astype("string").str.split(":").str[0].astype(int)
cb["CBIO_DAYS"] = pd.to_numeric(cb["PFS_MONTHS"], errors="coerce") * 30.417
cb.to_csv(OUT / "data_processed/cbio_PFS_label_endpoint_comparison.csv", index=False)
assert cb["CBIO_EVENT"].eq(cb["PFI"]).all()
valid_time = cb["CBIO_DAYS"].notna() & cb["PFI.time"].notna()
assert np.allclose(cb.loc[valid_time, "CBIO_DAYS"], cb.loc[valid_time, "PFI.time"], atol=0.001, rtol=0)

model_rows, coefficient_rows, ph_rows, warning_rows, cohort_rows, logrank_rows = [], [], [], [], [], []
km_models = {}

def encode(frame):
    dummies = pd.get_dummies(pd.Categorical(frame["SUBTYPE"], categories=SUBTYPES),
                            prefix="SUBTYPE", dtype=int)
    dummies.index = frame.index
    return pd.concat([frame, dummies.drop(columns="SUBTYPE_UCEC_CN_HIGH")], axis=1)

def fit(frame, endpoint, exposure, model_name, covariates, strata=None):
    duration_column = f"{endpoint}_TIME" if endpoint == "OS" else f"{endpoint}_DAYS"
    cols = [duration_column, f"{endpoint}_EVENT", *covariates, *(strata or [])]
    data = frame[list(dict.fromkeys(cols))].copy()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = CoxPHFitter()
        model.fit(data, duration_col=duration_column, event_col=f"{endpoint}_EVENT",
                  strata=strata, show_progress=False)
        ph = proportional_hazard_test(model, data, time_transform="rank")
    for warning in caught:
        warning_rows.append({"endpoint": endpoint, "exposure": exposure, "model": model_name,
                             "category": warning.category.__name__, "message": str(warning.message)})
    row = model.summary.loc[exposure]
    definition = "Primary portal definition" if exposure == "Altered" else "Amplification excluded"
    info = {"endpoint": endpoint, "definition": definition, "model": model_name,
            "patients_n": len(data), "events_n": int(data[f"{endpoint}_EVENT"].sum()),
            "altered_n": int(frame[exposure].sum()), "PTEN_HR": float(row["exp(coef)"]),
            "CI_95_low": float(row["exp(coef) lower 95%"]),
            "CI_95_high": float(row["exp(coef) upper 95%"]), "PTEN_p_value": float(row["p"]),
            "PTEN_PH_p_value": float(ph.summary.loc[exposure, "p"]),
            "concordance_index": float(model.concordance_index_), "penalizer": 0,
            "strata": ",".join(strata or []), "model_time_unit": "months" if endpoint == "OS" else "days"}
    model_rows.append(info)
    coefs = model.summary.reset_index().rename(columns={"covariate": "variable"})
    coefs["endpoint"], coefs["definition"], coefs["model"] = endpoint, definition, model_name
    coefficient_rows.append(coefs)
    diagnostics = ph.summary.reset_index().rename(columns={"index": "variable"})
    diagnostics["endpoint"], diagnostics["definition"], diagnostics["model"] = endpoint, definition, model_name
    ph_rows.append(diagnostics)
    return info

for endpoint in ["OS", "PFS", "PFI"]:
    cohort = os.copy() if endpoint == "OS" else os[os[f"{endpoint}_exclusion"].eq("")].copy()
    cohort[f"{endpoint}_EVENT"] = cohort[f"{endpoint}_EVENT"].astype(int)
    complete = cohort[cohort["complete_case"]].copy()
    cohort.to_csv(OUT / f"data_processed/{endpoint}_eligible_cohort.csv", index=False)
    complete.to_csv(OUT / f"data_processed/{endpoint}_complete_case_cohort.csv", index=False)
    for label, frame in [("Endpoint-eligible", cohort), ("Fixed complete case", complete)]:
        cohort_rows.append({"endpoint": endpoint, "cohort": label, "patients_n": len(frame),
                            "events_n": int(frame[f"{endpoint}_EVENT"].sum()),
                            "altered_n": int(frame["Altered"].sum()),
                            "unaltered_n": int(frame["Altered"].eq(0).sum()),
                            "altered_events_n": int(frame.loc[frame["Altered"].eq(1), f"{endpoint}_EVENT"].sum()),
                            "unaltered_events_n": int(frame.loc[frame["Altered"].eq(0), f"{endpoint}_EVENT"].sum()),
                            "zero_times_n": int(frame[f"{endpoint}_TIME"].eq(0).sum())})
    for exposure in ["Altered", "PTEN_NO_AMP"]:
        fit(cohort, endpoint, exposure, "Unadjusted: all eligible", [exposure])
        encoded = encode(complete)
        fit(encoded, endpoint, exposure, "PTEN only: fixed complete cases", [exposure])
        fit(encoded, endpoint, exposure, "Age and stage: fixed complete cases", [exposure, "AGE_10Y", "ADVANCED_STAGE"])
        fit(encoded, endpoint, exposure, "Age, stage and subtype: fixed complete cases",
            [exposure, "AGE_10Y", "ADVANCED_STAGE", "SUBTYPE_UCEC_CN_LOW", "SUBTYPE_UCEC_MSI", "SUBTYPE_UCEC_POLE"])
        fit(complete, endpoint, exposure, "Stage/subtype stratified: fixed complete cases",
            [exposure, "AGE_10Y"], strata=["ADVANCED_STAGE", "SUBTYPE"])
        altered = cohort[exposure].eq(1)
        duration_column = f"{endpoint}_TIME" if endpoint == "OS" else f"{endpoint}_DAYS"
        lr = logrank_test(cohort.loc[altered, duration_column], cohort.loc[~altered, duration_column],
                          event_observed_A=cohort.loc[altered, f"{endpoint}_EVENT"],
                          event_observed_B=cohort.loc[~altered, f"{endpoint}_EVENT"])
        logrank_rows.append({"endpoint": endpoint, "exposure": exposure, "patients_n": len(cohort),
                            "chi_square": float(lr.test_statistic), "p_value": float(lr.p_value)})
    km_models[endpoint] = []
    for value, label in [(1, "PTEN-altered"), (0, "PTEN-unaltered")]:
        group = cohort[cohort["Altered"].eq(value)]
        km = KaplanMeierFitter(label=f"{label} (n={len(group)})").fit(
            group[f"{endpoint}_TIME"], event_observed=group[f"{endpoint}_EVENT"])
        km_models[endpoint].append((km, group))

models = pd.DataFrame(model_rows)
for name, expected in [("Unadjusted: all eligible", 0.593361),
                       ("PTEN only: fixed complete cases", 0.595209),
                       ("Age and stage: fixed complete cases", 0.825825),
                       ("Age, stage and subtype: fixed complete cases", 1.071310),
                       ("Stage/subtype stratified: fixed complete cases", 1.138520)]:
    observed = models.loc[models["endpoint"].eq("OS") & models["definition"].eq("Primary portal definition") & models["model"].eq(name), "PTEN_HR"].iloc[0]
    assert np.isclose(observed, expected, atol=2e-6, rtol=0), (name, observed, expected)

models.to_csv(OUT / "results/model_comparison.csv", index=False)
pd.concat(coefficient_rows, ignore_index=True).to_csv(OUT / "results/all_cox_coefficients.csv", index=False)
pd.concat(ph_rows, ignore_index=True).to_csv(OUT / "results/proportional_hazards_rank_tests.csv", index=False)
pd.DataFrame(cohort_rows).to_csv(OUT / "results/cohort_counts.csv", index=False)
pd.DataFrame(logrank_rows).to_csv(OUT / "results/logrank_tests.csv", index=False)
pd.DataFrame(warning_rows, columns=["endpoint", "exposure", "model", "category", "message"]).to_csv(
    OUT / "results/model_warnings.csv", index=False)

# Publication figure: PFS requested by reviewer and progression-focused PFI.
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42,
                     "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(2, 1, figsize=(6.5, 6.3), sharey=True)
colors = ["#0072B2", "#D55E00"]
risk_rows = []
for ax, endpoint, letter in zip(axes, ["PFS", "PFI"], ["A", "B"]):
    # Display the same 120-month window as original Figure 1. All follow-up
    # remains in the estimators and tests, including times beyond 120 months.
    xmax = 120
    ticks = np.arange(0, xmax + 1, 24)
    for (km, group), color in zip(km_models[endpoint], colors):
        km.plot_survival_function(ax=ax, ci_show=True, color=color, linewidth=1.8,
                                  show_censors=True, censor_styles={"ms": 2.5, "marker": "+"})
        at_risk = [int(group[f"{endpoint}_TIME"].ge(t).sum()) for t in ticks]
        risk_rows.extend({"endpoint": endpoint, "group": km._label, "time_months": float(t),
                          "at_risk_at_start_of_time": n} for t, n in zip(ticks, at_risk))
    lr = next(row for row in logrank_rows if row["endpoint"] == endpoint and row["exposure"] == "Altered")
    ax.set(xlim=(0, xmax), ylim=(0, 1.02), xticks=ticks, xlabel="Time from diagnosis (months)")
    ax.set_title(f"{letter}  {'Progression-free survival' if endpoint == 'PFS' else 'Progression-free interval'}", loc="left", fontsize=11)
    ax.text(0.97, 0.95, f"Log-rank p = {lr['p_value']:.3f}", ha="right", va="top", transform=ax.transAxes)
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    ax.grid(alpha=0.18)
    ax.text(0, -0.34, "Number at risk", transform=ax.transAxes, ha="left", fontsize=9)
    for j, ((km, group), color) in enumerate(zip(km_models[endpoint], colors)):
        y = -0.45 - j * 0.12
        ax.text(-0.03, y, "Altered" if j == 0 else "Unaltered", transform=ax.transAxes,
                ha="right", va="center", fontsize=9, color=color)
        for t in ticks:
            count = int(group[f"{endpoint}_TIME"].ge(t).sum())
            ax.text(t / xmax, y, str(count), transform=ax.transAxes,
                    ha="center", va="center", fontsize=9)
for ax in axes:
    ax.set_ylabel("Event-free probability", fontsize=9)
fig.subplots_adjust(left=0.19, right=0.96, bottom=0.19, top=0.95, hspace=1.00)
fig.savefig(OUT / "figures/Figure4_PFS_PFI_Kaplan_Meier.png", dpi=300)
fig.savefig(OUT / "figures/Figure4_PFS_PFI_Kaplan_Meier.pdf")
plt.close(fig)
pd.DataFrame(risk_rows).to_csv(OUT / "results/Figure4_numbers_at_risk.csv", index=False)

source_counts = []
for frame, endpoint, time_column, sheet in [(cdr, "PFI", "PFI.time", "TCGA-CDR"),
                                          (extra, "PFS", "PFS.time", "ExtraEndpoints")]:
    source = frame[frame["type"].eq("UCEC")]
    usable = source[endpoint].isin([0, 1]) & source[time_column].notna() & source[time_column].ge(0)
    source_counts.append({"source_sheet": sheet, "endpoint": endpoint, "source_UCEC_n": len(source),
                          "source_usable_endpoint_n": int(usable.sum()),
                          "source_usable_events_n": int(source.loc[usable, endpoint].sum()),
                          "source_zero_time_n": int(source.loc[usable, time_column].eq(0).sum()),
                          "matched_original_OS_cohort_n": int(os[f"{endpoint}_MATCH_ID"].notna().sum()),
                          "eligible_primary_endpoint_cohort_n": int(os[f"{endpoint}_exclusion"].eq("").sum())})
pd.DataFrame(source_counts).to_csv(OUT / "results/source_endpoint_coverage.csv", index=False)

summary = {"software": VERSIONS, "original_os_reproduced": True,
           "input_sha256": hash_rows, "cohort_counts": cohort_rows,
           "models": model_rows, "logrank": logrank_rows,
           "warnings": warning_rows,
           "convergence_warning_count": sum("convergence" in row["category"].lower() for row in warning_rows),
           "source_endpoint_counts": source_counts, "amplification_only_sample_ids": AMPLIFICATION_ONLY,
           "endpoint_exclusions": {e: os[f"{e}_exclusion"].value_counts().to_dict() for e in ["PFS", "PFI"]},
           "covariate_missing": {"age_n": int(os["AGE_10Y"].isna().sum()),
                                  "stage_n": int(os["ADVANCED_STAGE"].isna().sum()),
                                  "subtype_n": int(os["SUBTYPE"].isna().sum())},
           "endpoint_label_verification": {"patients_compared": len(cb),
                                          "cbio_events_matching_PFI": int(cb["CBIO_EVENT"].eq(cb["PFI"]).sum()),
                                          "cbio_events_differing_from_true_PFS": int(cb["CBIO_EVENT"].ne(cb["PFS"]).sum()),
                                          "nonmissing_time_matches": int(valid_time.sum()), "historical_cbio_days_per_month": 30.417},
           "new_endpoint_month_conversion": "For KM display only: days divided by 365.25/12 (30.4375)",
           "cox_and_logrank_units": {"OS": "original cBioPortal months", "PFS": "original TCGA-CDR days", "PFI": "original TCGA-CDR days"},
           "age_covariate_definition": "(AGE - 64)/10; hazard ratio per 10-year increase",
           "figure_display_window_months": [0, 120], "all_followup_used_in_models": True,
           "analysis_status": "Additional analysis requested during review; exploratory secondary endpoints",
           "limitations": ["Retrospective observational association; residual confounding remains.",
                           "No treatment covariates or adjudicated causes of death were added.",
                           "PFS includes progression or death from any cause; PFI includes new tumor events or death with tumor.",
                           "PFI is not a confirmed cancer-specific mortality endpoint.",
                           "Complete-case analysis excludes 39 patients with missing age and/or subtype.",
                           "P-values are unadjusted for multiple secondary and sensitivity analyses.",
                           "Stage/subtype stratification uses sparse strata; its concordance index is not directly comparable with unstratified models."]}
(OUT / "results/analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(pd.DataFrame(cohort_rows).to_string(index=False))
print(models.loc[models["definition"].eq("Primary portal definition")].to_string(index=False))
print("Warnings captured:", len(warning_rows))
