#!/usr/bin/env python3
"""Describe the fixed OS cohorts for JEI table/figure revisions.

OS_TIME is the observed time to death or censoring in months. Its descriptive
median is not a Kaplan-Meier survival median or reverse-KM median follow-up.
Only descriptive statistics and product-limit KM steps are calculated; no Cox
or other regression model is fitted.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent


def numeric_summary(series):
    values = pd.to_numeric(series, errors="coerce")
    valid = values[np.isfinite(values)].copy()
    quartiles = valid.quantile([0.25, 0.5, 0.75], interpolation="linear")
    return {
        "available_n": int(len(valid)),
        "missing_or_nonnumeric_n": int(len(values) - len(valid)),
        "minimum": float(valid.min()),
        "q1": float(quartiles.loc[0.25]),
        "median": float(quartiles.loc[0.5]),
        "q3": float(quartiles.loc[0.75]),
        "iqr_width": float(quartiles.loc[0.75] - quartiles.loc[0.25]),
        "maximum": float(valid.max()),
    }


def km_steps(group):
    """Product-limit survival at each observed time, using pre-time risk sets."""
    table = group.groupby("OS_TIME")["OS_EVENT"].agg(
        removed_n="size", deaths_n="sum").reset_index().sort_values("OS_TIME")
    table["censored_n"] = table["removed_n"] - table["deaths_n"]
    table["at_risk_before_time_n"] = len(group) - table["removed_n"].cumsum().shift(fill_value=0)
    table["conditional_survival"] = 1 - table["deaths_n"] / table["at_risk_before_time_n"]
    table["KM_survival"] = table["conditional_survival"].cumprod()
    table["KM_survival_before_time"] = table["KM_survival"].shift(fill_value=1.0)
    assert table["at_risk_before_time_n"].gt(0).all()
    assert table["KM_survival"].between(0, 1).all()
    assert table["KM_survival"].diff().fillna(0).le(1e-12).all()
    assert int(table["removed_n"].sum()) == len(group)
    return table


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--os", type=Path,
                        default=ROOT / "data_processed/OS_eligible_cohort.csv")
    parser.add_argument("--audit", type=Path,
                        default=ROOT / "data_processed/all_498_patients_endpoint_and_covariate_audit.csv")
    parser.add_argument("--clinical", type=Path,
                        default=ROOT / "data_raw/Supplementary_Table_S4_cBioPortal_clinical_survival_table.tsv")
    parser.add_argument("--out", type=Path, default=ROOT / "results")
    args = parser.parse_args()

    os = pd.read_csv(args.os)
    audit = pd.read_csv(args.audit)
    raw = pd.read_csv(args.clinical, sep="\t")
    assert (len(os), int(os["OS_EVENT"].sum()), int(os["Altered"].sum())) == (498, 85, 146)
    assert os["patientId"].is_unique and os["sampleId"].is_unique
    assert os["OS_TIME"].notna().all() and os["OS_TIME"].ge(0).all()
    assert os["OS_EVENT"].isin([0, 1]).all() and os["Altered"].isin([0, 1]).all()
    assert audit["patientId"].is_unique and set(audit["patientId"]) == set(os["patientId"])
    columns = ["patientId", "sampleId", "AGE", "OS_TIME", "OS_EVENT", "Altered", "complete_case"]
    pd.testing.assert_frame_equal(
        os[columns].sort_values("patientId").reset_index(drop=True),
        audit[columns].sort_values("patientId").reset_index(drop=True),
        check_dtype=False)

    # Validate the used ages/durations against the archived Firehose export.
    assert raw["patientId"].is_unique and raw["sampleId"].is_unique
    raw_compare = os[["patientId", "sampleId", "AGE", "OS_TIME", "OS_EVENT"]].merge(
        raw[["patientId", "sampleId", "AGE", "OS_MONTHS", "OS_STATUS"]],
        on=["patientId", "sampleId"], how="left", suffixes=("", "_raw"),
        validate="one_to_one", indicator=True)
    assert raw_compare["_merge"].eq("both").all()
    np.testing.assert_allclose(pd.to_numeric(raw_compare["AGE"], errors="coerce"),
                               pd.to_numeric(raw_compare["AGE_raw"], errors="coerce"),
                               equal_nan=True, atol=0, rtol=0)
    np.testing.assert_allclose(raw_compare["OS_TIME"],
                               pd.to_numeric(raw_compare["OS_MONTHS"], errors="coerce"),
                               equal_nan=True, atol=0, rtol=0)
    raw_events = raw_compare["OS_STATUS"].astype("string").str.split(":").str[0].map({"0": 0, "1": 1})
    assert raw_events.eq(raw_compare["OS_EVENT"]).all()

    complete = os[os["complete_case"].eq(True)].copy()
    assert (len(complete), int(complete["OS_EVENT"].sum()), int(complete["Altered"].sum())) == (459, 77, 139)
    descriptions, km_summaries, step_tables = [], [], []
    for cohort_name, frame in [("498-patient OS cohort", os), ("459 complete cases", complete)]:
        for group_name, group in [("Overall", frame),
                                  ("PTEN-altered", frame[frame["Altered"].eq(1)]),
                                  ("PTEN-unaltered", frame[frame["Altered"].eq(0)])]:
            age = numeric_summary(group["AGE"])
            duration = numeric_summary(group["OS_TIME"])
            row = {"cohort": cohort_name, "group": group_name, "patients_n": len(group),
                   "deaths_n": int(group["OS_EVENT"].sum()),
                   "censored_n": int(group["OS_EVENT"].eq(0).sum())}
            row.update({"age_years_" + key: value for key, value in age.items()})
            row.update({"observed_duration_months_" + key: value for key, value in duration.items()})
            descriptions.append(row)
            table = km_steps(group)
            crossing = table[table["KM_survival"].le(0.5)]
            first = crossing.iloc[0] if not crossing.empty else None
            km_row = {
                "cohort": cohort_name, "group": group_name, "patients_n": len(group),
                "deaths_n": int(group["OS_EVENT"].sum()),
                "censored_n": int(group["OS_EVENT"].eq(0).sum()),
                "minimum_observed_duration_months": float(group["OS_TIME"].min()),
                "maximum_observed_duration_months": float(group["OS_TIME"].max()),
                "last_observed_death_months": float(group.loc[group["OS_EVENT"].eq(1), "OS_TIME"].max()),
                "KM_final_survival": float(table["KM_survival"].iloc[-1]),
                "KM_minimum_survival_including_initial_1": float(table["KM_survival"].min()),
                "KM_maximum_survival_including_initial_1": 1.0,
                "KM_ever_at_or_below_half": not crossing.empty,
                "KM_median_months": None if first is None else float(first["OS_TIME"]),
                "risk_set_at_median_n": None if first is None else int(first["at_risk_before_time_n"]),
                "deaths_at_median_n": None if first is None else int(first["deaths_n"]),
                "extrapolation_used": False,
            }
            km_summaries.append(km_row)
            table.insert(0, "group", group_name)
            table.insert(0, "cohort", cohort_name)
            step_tables.append(table)
            if cohort_name == "498-patient OS cohort" and group_name == "PTEN-altered":
                assert len(group) == 146 and int(group["OS_EVENT"].sum()) == 20
                assert not km_row["KM_ever_at_or_below_half"]
                assert km_row["KM_final_survival"] > 0.5
                assert km_row["maximum_observed_duration_months"] == 129.7
            if cohort_name == "498-patient OS cohort" and group_name == "PTEN-unaltered":
                assert (km_row["KM_median_months"], km_row["risk_set_at_median_n"],
                        km_row["deaths_at_median_n"]) == (110.55, 10, 1)
    assert descriptions[0]["age_years_available_n"] == 496
    assert descriptions[0]["age_years_missing_or_nonnumeric_n"] == 2
    assert descriptions[0]["observed_duration_months_minimum"] == 0.07
    assert descriptions[0]["observed_duration_months_maximum"] == 225.33

    args.out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(descriptions).to_csv(args.out / "OS_age_and_observed_duration_summary.csv", index=False)
    pd.DataFrame(km_summaries).to_csv(args.out / "OS_Kaplan_Meier_summary.csv", index=False)
    pd.concat(step_tables, ignore_index=True).to_csv(args.out / "OS_Kaplan_Meier_steps.csv", index=False)
    summary = {
        "definitions": {
            "cohort": "Fixed 498-patient original OS cohort; 459 complete cases are a subset.",
            "age": "AGE from the archived cBioPortal Firehose Legacy clinical-survival export, in years.",
            "observed_duration": "OS_TIME / OS_MONTHS: observed time to death or last follow-up/censoring, in months.",
            "quantiles": "pandas Series.quantile with linear interpolation; IQR bounds are Q1 and Q3.",
            "iqr_width": "Q3 minus Q1; distinguish the width from the interval Q1–Q3.",
            "KM_median": "First observed time with product-limit survival <=0.5; not reached is JSON null.",
            "KM_survival_range": "Includes the initial survival probability of 1 at time zero.",
            "reverse_KM_follow_up": "Not estimated. Observed-duration median is not reverse-KM median follow-up.",
            "models": "No Cox or other regression model was fitted.",
        },
        "source_checks": {
            "fixed_cohort_and_audit_agree": True,
            "all_498_match_archived_clinical_export": True,
            "age_OS_time_and_event_agree_with_archived_export": True,
        },
        "inputs": [{"name": path.name, "path": str(path), "bytes": path.stat().st_size,
                    "sha256": sha256(path)} for path in [args.os, args.audit, args.clinical]],
        "software": {"Python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__},
        "descriptive_statistics": descriptions,
        "Kaplan_Meier_summary": km_summaries,
    }
    (args.out / "OS_descriptive_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"out": str(args.out), "descriptive_statistics": descriptions,
                      "Kaplan_Meier_summary": km_summaries[:3]}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
