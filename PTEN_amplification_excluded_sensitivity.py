from pathlib import Path

import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test


project_root = Path(globals().get("PROJECT_ROOT", Path.cwd()))
processed_dir = project_root / "data_processed"
results_dir = project_root / "results"
results_dir.mkdir(parents=True, exist_ok=True)

reconstruction_path = processed_dir / "PTEN_status_reconstruction.csv"
survival_path = processed_dir / "survival_analysis_dataset.csv"
complete_path = processed_dir / "adjusted_cox_complete_case_dataset.csv"

for required_path in [reconstruction_path, survival_path, complete_path]:
    if not required_path.exists():
        raise FileNotFoundError(
            f"Missing {required_path}. Run the master notebook through Phase 3 first."
        )

# Reconstruct the narrower sensitivity definition:
# PTEN mutation OR GISTIC -2 deep deletion; GISTIC +2 alone is excluded.
reconstruction = pd.read_csv(reconstruction_path)
mutation_present = (
    reconstruction["MUTATION_PRESENT"]
    .astype("string")
    .str.strip()
    .str.lower()
    .eq("true")
)
pten_cna = pd.to_numeric(reconstruction["PTEN_CNA"], errors="coerce")
reconstruction["PTEN_NO_AMP"] = (
    mutation_present | pten_cna.eq(-2)
).astype(int)
reconstruction["PRIMARY_ALT"] = pd.to_numeric(
    reconstruction["RECONSTRUCTED_ALTERED"], errors="raise"
).astype(int)
reconstruction["RECLASSIFIED"] = (
    reconstruction["PRIMARY_ALT"] != reconstruction["PTEN_NO_AMP"]
)

reclassified = reconstruction.loc[
    reconstruction["RECLASSIFIED"],
    [
        "SAMPLE_ID",
        "PTEN_CNA",
        "MUTATION_PRESENT",
        "PRIMARY_ALT",
        "PTEN_NO_AMP",
    ],
].copy()

if len(reconstruction) != 549:
    raise ValueError(f"Expected 549 reconstructed tumors, found {len(reconstruction)}.")
if int(reconstruction["PRIMARY_ALT"].sum()) != 181:
    raise ValueError("Primary altered count does not reproduce 181.")
if int(reconstruction["PTEN_NO_AMP"].sum()) != 178:
    raise ValueError("Expected 178 altered tumors after excluding amplification alone.")
if len(reclassified) != 3:
    raise ValueError(f"Expected exactly 3 reclassified tumors, found {len(reclassified)}.")

# Add the sensitivity classification to the existing 498-patient survival cohort.
status_lookup = reconstruction[["SAMPLE_ID", "PTEN_NO_AMP"]].copy()
survival = pd.read_csv(survival_path)
survival = survival.merge(
    status_lookup,
    left_on="sampleId",
    right_on="SAMPLE_ID",
    how="left",
    validate="one_to_one",
)
if survival["PTEN_NO_AMP"].isna().any():
    raise ValueError("Some survival patients did not receive a sensitivity PTEN status.")
survival["PTEN_NO_AMP"] = survival["PTEN_NO_AMP"].astype(int)
survival["Altered"] = survival["Altered"].astype(int)

if len(survival) != 498 or int(survival["OS_EVENT"].sum()) != 85:
    raise ValueError("The survival cohort or event count changed unexpectedly.")
if int(survival["PTEN_NO_AMP"].sum()) != 143:
    raise ValueError("Expected 143 altered patients in the sensitivity survival cohort.")

# Add the same classification to the fixed 459-patient complete-case cohort.
patient_status = survival[["patientId", "PTEN_NO_AMP"]].drop_duplicates()
if patient_status["patientId"].duplicated().any():
    raise ValueError("Conflicting sensitivity labels were found for a patient.")

complete = pd.read_csv(complete_path).merge(
    patient_status,
    on="patientId",
    how="left",
    validate="one_to_one",
)
if complete["PTEN_NO_AMP"].isna().any():
    raise ValueError("Some complete-case patients did not receive a sensitivity PTEN status.")
complete["PTEN_NO_AMP"] = complete["PTEN_NO_AMP"].astype(int)
complete["Altered"] = complete["Altered"].astype(int)
complete["OS_EVENT"] = complete["OS_EVENT"].astype(int)
complete["ADVANCED_STAGE"] = complete["ADVANCED_STAGE"].astype(int)

if len(complete) != 459 or int(complete["OS_EVENT"].sum()) != 77:
    raise ValueError("The complete-case cohort or event count changed unexpectedly.")
if int(complete["PTEN_NO_AMP"].sum()) != 136:
    raise ValueError("Expected 136 altered patients in the sensitivity complete-case cohort.")

# Preserve string subtype labels for the stratified models.
complete_stratified = complete.copy()

# CN-high remains the reference subtype in the fully adjusted models.
subtype_order = ["UCEC_CN_HIGH", "UCEC_CN_LOW", "UCEC_MSI", "UCEC_POLE"]
complete["SUBTYPE"] = pd.Categorical(complete["SUBTYPE"], categories=subtype_order)
subtype_dummies = pd.get_dummies(
    complete["SUBTYPE"], prefix="SUBTYPE", dtype=int
).drop(columns=["SUBTYPE_UCEC_CN_HIGH"])
complete_encoded = pd.concat(
    [complete.drop(columns=["SUBTYPE"]), subtype_dummies], axis=1
)


def fit_and_extract(dataframe, exposure, model_name, covariates, strata=None):
    model_columns = ["OS_MONTHS_NUMERIC", "OS_EVENT", *covariates]
    if strata:
        model_columns.extend([column for column in strata if column not in model_columns])
    model_data = dataframe[model_columns].copy()

    model = CoxPHFitter()
    model.fit(
        model_data,
        duration_col="OS_MONTHS_NUMERIC",
        event_col="OS_EVENT",
        strata=strata,
        show_progress=False,
    )
    row = model.summary.loc[exposure]
    ph = proportional_hazard_test(model, model_data, time_transform="rank")
    ph_p = float(ph.summary.loc[exposure, "p"])

    result = {
        "model": model_name,
        "definition": "Primary portal definition" if exposure == "Altered" else "Amplification excluded",
        "patients_n": len(model_data),
        "deaths_n": int(model_data["OS_EVENT"].sum()),
        "altered_n": int(model_data[exposure].sum()),
        "PTEN_HR": float(row["exp(coef)"]),
        "CI_95_low": float(row["exp(coef) lower 95%"]),
        "CI_95_high": float(row["exp(coef) upper 95%"]),
        "PTEN_p_value": float(row["p"]),
        "PTEN_PH_p_value": ph_p,
        "concordance_index": float(model.concordance_index_),
    }
    return model, result


results = []

# Unadjusted models in the full 498-patient survival cohort.
for exposure in ["Altered", "PTEN_NO_AMP"]:
    _, result = fit_and_extract(
        survival,
        exposure,
        "Unadjusted Cox (498 patients)",
        [exposure],
    )
    results.append(result)

# PTEN-only models on the fixed 459-patient complete-case cohort.
for exposure in ["Altered", "PTEN_NO_AMP"]:
    _, result = fit_and_extract(
        complete_encoded,
        exposure,
        "PTEN-only Cox (459 complete cases)",
        [exposure],
    )
    results.append(result)

# Fully adjusted models on those same 459 patients.
for exposure in ["Altered", "PTEN_NO_AMP"]:
    adjusted_covariates = [
        exposure,
        "AGE_10Y",
        "ADVANCED_STAGE",
        "SUBTYPE_UCEC_CN_LOW",
        "SUBTYPE_UCEC_MSI",
        "SUBTYPE_UCEC_POLE",
    ]
    _, result = fit_and_extract(
        complete_encoded,
        exposure,
        "Fully adjusted Cox (459 complete cases)",
        adjusted_covariates,
    )
    results.append(result)

# Stage- and subtype-stratified sensitivity models.
for exposure in ["Altered", "PTEN_NO_AMP"]:
    _, result = fit_and_extract(
        complete_stratified,
        exposure,
        "Stage/subtype-stratified Cox (459 complete cases)",
        [exposure, "AGE_10Y"],
        strata=["ADVANCED_STAGE", "SUBTYPE"],
    )
    results.append(result)

comparison = pd.DataFrame(results)

# Log-rank comparison for the narrower definition in the full survival cohort.
sens_altered = survival["PTEN_NO_AMP"].eq(1)
logrank = logrank_test(
    survival.loc[sens_altered, "OS_MONTHS_NUMERIC"],
    survival.loc[~sens_altered, "OS_MONTHS_NUMERIC"],
    event_observed_A=survival.loc[sens_altered, "OS_EVENT"],
    event_observed_B=survival.loc[~sens_altered, "OS_EVENT"],
)
logrank_summary = pd.DataFrame(
    {
        "definition": ["Amplification excluded"],
        "patients_n": [len(survival)],
        "altered_n": [int(sens_altered.sum())],
        "unaltered_n": [int((~sens_altered).sum())],
        "deaths_n": [int(survival["OS_EVENT"].sum())],
        "logrank_chi_square": [float(logrank.test_statistic)],
        "logrank_p_value": [float(logrank.p_value)],
    }
)

reclassified.to_csv(
    results_dir / "amplification_excluded_reclassified_tumors.csv", index=False
)
comparison.to_csv(
    results_dir / "amplification_excluded_cox_comparison.csv", index=False
)
logrank_summary.to_csv(
    results_dir / "amplification_excluded_logrank.csv", index=False
)

print("Amplification-excluded validation counts")
print(
    pd.DataFrame(
        {
            "cohort": ["Reconstruction", "Survival", "Complete case"],
            "total_n": [len(reconstruction), len(survival), len(complete)],
            "primary_altered_n": [
                int(reconstruction["PRIMARY_ALT"].sum()),
                int(survival["Altered"].sum()),
                int(complete["Altered"].sum()),
            ],
            "no_amplification_altered_n": [
                int(reconstruction["PTEN_NO_AMP"].sum()),
                int(survival["PTEN_NO_AMP"].sum()),
                int(complete["PTEN_NO_AMP"].sum()),
            ],
        }
    ).to_string(index=False)
)
print("\nCox comparison")
print(comparison.to_string(index=False))
print("\nAmplification-excluded log-rank result")
print(logrank_summary.to_string(index=False))

try:
    display(comparison)
    display(logrank_summary)
except NameError:
    pass
