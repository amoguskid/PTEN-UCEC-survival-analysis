"""Compare the primary PTEN definition with one that excludes amplification alone."""

from pathlib import Path

import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test


RECONSTRUCTION_N = 549
SURVIVAL_N = 498
COMPLETE_CASE_N = 459
SUBTYPE_ORDER = ["UCEC_CN_HIGH", "UCEC_CN_LOW", "UCEC_MSI", "UCEC_POLE"]


def require_files(paths):
    """Fail with a useful message when the notebook outputs are unavailable."""
    missing = [path for path in paths if not path.exists()]
    if missing:
        formatted = ", ".join(str(path) for path in missing)
        raise FileNotFoundError(
            f"Missing required notebook output(s): {formatted}. "
            "Run the master notebook through Phase 3 first."
        )


def expect_count(actual, expected, label):
    """Keep cohort drift from silently changing the reported analysis."""
    if actual != expected:
        raise ValueError(f"Expected {expected} {label}, found {actual}.")


def build_sensitivity_status(reconstruction_path):
    reconstruction = pd.read_csv(reconstruction_path)
    mutation_present = (
        reconstruction["MUTATION_PRESENT"]
        .astype("string")
        .str.strip()
        .str.lower()
        .eq("true")
    )
    pten_cna = pd.to_numeric(reconstruction["PTEN_CNA"], errors="coerce")

    reconstruction["PTEN_NO_AMP"] = (mutation_present | pten_cna.eq(-2)).astype(int)
    reconstruction["PRIMARY_ALT"] = pd.to_numeric(
        reconstruction["RECONSTRUCTED_ALTERED"], errors="raise"
    ).astype(int)
    reconstruction["RECLASSIFIED"] = (
        reconstruction["PRIMARY_ALT"] != reconstruction["PTEN_NO_AMP"]
    )

    reclassified = reconstruction.loc[
        reconstruction["RECLASSIFIED"],
        ["SAMPLE_ID", "PTEN_CNA", "MUTATION_PRESENT", "PRIMARY_ALT", "PTEN_NO_AMP"],
    ].copy()

    expect_count(len(reconstruction), RECONSTRUCTION_N, "reconstructed tumors")
    expect_count(int(reconstruction["PRIMARY_ALT"].sum()), 181, "primary altered tumors")
    expect_count(
        int(reconstruction["PTEN_NO_AMP"].sum()),
        178,
        "altered tumors after excluding amplification alone",
    )
    expect_count(len(reclassified), 3, "reclassified tumors")
    return reconstruction, reclassified


def build_analysis_cohorts(reconstruction, survival_path, complete_path):
    status_lookup = reconstruction[["SAMPLE_ID", "PTEN_NO_AMP"]]
    survival = pd.read_csv(survival_path).merge(
        status_lookup,
        left_on="sampleId",
        right_on="SAMPLE_ID",
        how="left",
        validate="one_to_one",
    )
    if survival["PTEN_NO_AMP"].isna().any():
        raise ValueError("Some survival patients lack a sensitivity PTEN status.")
    survival[["PTEN_NO_AMP", "Altered"]] = survival[
        ["PTEN_NO_AMP", "Altered"]
    ].astype(int)

    expect_count(len(survival), SURVIVAL_N, "survival patients")
    expect_count(int(survival["OS_EVENT"].sum()), 85, "survival deaths")
    expect_count(int(survival["PTEN_NO_AMP"].sum()), 143, "sensitivity altered patients")

    patient_status = survival[["patientId", "PTEN_NO_AMP"]].drop_duplicates()
    if patient_status["patientId"].duplicated().any():
        raise ValueError("Conflicting sensitivity labels were found for a patient.")

    complete = pd.read_csv(complete_path).merge(
        patient_status, on="patientId", how="left", validate="one_to_one"
    )
    if complete["PTEN_NO_AMP"].isna().any():
        raise ValueError("Some complete-case patients lack a sensitivity PTEN status.")
    integer_columns = ["PTEN_NO_AMP", "Altered", "OS_EVENT", "ADVANCED_STAGE"]
    complete[integer_columns] = complete[integer_columns].astype(int)

    expect_count(len(complete), COMPLETE_CASE_N, "complete-case patients")
    expect_count(int(complete["OS_EVENT"].sum()), 77, "complete-case deaths")
    expect_count(
        int(complete["PTEN_NO_AMP"].sum()), 136, "complete-case sensitivity patients"
    )
    return survival, complete


def encode_subtypes(complete):
    encoded = complete.copy()
    encoded["SUBTYPE"] = pd.Categorical(encoded["SUBTYPE"], categories=SUBTYPE_ORDER)
    subtype_dummies = pd.get_dummies(
        encoded["SUBTYPE"], prefix="SUBTYPE", dtype=int
    ).drop(columns=["SUBTYPE_UCEC_CN_HIGH"])
    return pd.concat([encoded.drop(columns=["SUBTYPE"]), subtype_dummies], axis=1)


def fit_and_extract(dataframe, exposure, model_name, covariates, strata=None):
    model_columns = ["OS_MONTHS_NUMERIC", "OS_EVENT", *covariates]
    if strata:
        model_columns.extend(column for column in strata if column not in model_columns)
    model_data = dataframe[model_columns].copy()

    model = CoxPHFitter()
    model.fit(
        model_data,
        duration_col="OS_MONTHS_NUMERIC",
        event_col="OS_EVENT",
        strata=strata,
        show_progress=False,
    )
    coefficient = model.summary.loc[exposure]
    ph_test = proportional_hazard_test(model, model_data, time_transform="rank")

    return {
        "model": model_name,
        "definition": (
            "Primary portal definition" if exposure == "Altered" else "Amplification excluded"
        ),
        "patients_n": len(model_data),
        "deaths_n": int(model_data["OS_EVENT"].sum()),
        "altered_n": int(model_data[exposure].sum()),
        "PTEN_HR": float(coefficient["exp(coef)"]),
        "CI_95_low": float(coefficient["exp(coef) lower 95%"]),
        "CI_95_high": float(coefficient["exp(coef) upper 95%"]),
        "PTEN_p_value": float(coefficient["p"]),
        "PTEN_PH_p_value": float(ph_test.summary.loc[exposure, "p"]),
        "concordance_index": float(model.concordance_index_),
    }


def compare_cox_models(survival, complete):
    complete_encoded = encode_subtypes(complete)
    model_specs = [
        ("Unadjusted Cox (498 patients)", survival, lambda exposure: [exposure], None),
        (
            "PTEN-only Cox (459 complete cases)",
            complete_encoded,
            lambda exposure: [exposure],
            None,
        ),
        (
            "Fully adjusted Cox (459 complete cases)",
            complete_encoded,
            lambda exposure: [
                exposure,
                "AGE_10Y",
                "ADVANCED_STAGE",
                "SUBTYPE_UCEC_CN_LOW",
                "SUBTYPE_UCEC_MSI",
                "SUBTYPE_UCEC_POLE",
            ],
            None,
        ),
        (
            "Stage/subtype-stratified Cox (459 complete cases)",
            complete,
            lambda exposure: [exposure, "AGE_10Y"],
            ["ADVANCED_STAGE", "SUBTYPE"],
        ),
    ]
    results = []
    for model_name, cohort, covariates, strata in model_specs:
        for exposure in ["Altered", "PTEN_NO_AMP"]:
            results.append(
                fit_and_extract(
                    cohort, exposure, model_name, covariates(exposure), strata=strata
                )
            )
    return pd.DataFrame(results)


def run_analysis(project_root=None):
    project_root = Path(project_root or globals().get("PROJECT_ROOT", Path.cwd()))
    processed_dir = project_root / "data_processed"
    results_dir = project_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    reconstruction_path = processed_dir / "PTEN_status_reconstruction.csv"
    survival_path = processed_dir / "survival_analysis_dataset.csv"
    complete_path = processed_dir / "adjusted_cox_complete_case_dataset.csv"
    require_files([reconstruction_path, survival_path, complete_path])

    reconstruction, reclassified = build_sensitivity_status(reconstruction_path)
    survival, complete = build_analysis_cohorts(
        reconstruction, survival_path, complete_path
    )
    comparison = compare_cox_models(survival, complete)

    sensitivity_group = survival["PTEN_NO_AMP"].eq(1)
    logrank = logrank_test(
        survival.loc[sensitivity_group, "OS_MONTHS_NUMERIC"],
        survival.loc[~sensitivity_group, "OS_MONTHS_NUMERIC"],
        event_observed_A=survival.loc[sensitivity_group, "OS_EVENT"],
        event_observed_B=survival.loc[~sensitivity_group, "OS_EVENT"],
    )
    logrank_summary = pd.DataFrame(
        {
            "definition": ["Amplification excluded"],
            "patients_n": [len(survival)],
            "altered_n": [int(sensitivity_group.sum())],
            "unaltered_n": [int((~sensitivity_group).sum())],
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

    counts = pd.DataFrame(
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
    )
    print("Amplification-excluded validation counts")
    print(counts.to_string(index=False))
    print("\nCox comparison")
    print(comparison.to_string(index=False))
    print("\nAmplification-excluded log-rank result")
    print(logrank_summary.to_string(index=False))
    return comparison, logrank_summary


if __name__ == "__main__":
    run_analysis()
