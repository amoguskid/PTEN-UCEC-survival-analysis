#!/usr/bin/env python3
"""Describe archived mutation records and audit the existing OS Kaplan-Meier median.

No Cox models are fitted. Mutation counts do not establish loss of function.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--mutations", type=Path, default=ROOT / "reference_inputs/Supplementary_Table_S3_PTEN_mutation_table.tsv")
parser.add_argument("--labels", type=Path, default=ROOT / "data_raw/Supplementary_Table_S1_PTEN_altered_unaltered_sample_matrix.tsv")
parser.add_argument("--os", type=Path, default=ROOT / "data_processed/OS_eligible_cohort.csv")
parser.add_argument("--out", type=Path, default=ROOT / "results")
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=True)

mutations = pd.read_csv(args.mutations, sep="\t")
labels = pd.read_csv(args.labels, sep="\t")
os = pd.read_csv(args.os)
labels["sampleId"] = labels["studyID:sampleId"].str.split(":", n=1).str[-1]
assert len(labels) == 549 and labels["sampleId"].is_unique
assert len(mutations) == 255 and mutations["Gene"].eq("PTEN").all()
assert mutations["Sample ID"].nunique() == 161
assert set(mutations["Sample ID"]).issubset(set(labels["sampleId"]))
assert (len(os), int(os["OS_EVENT"].sum()), int(os["Altered"].sum())) == (498, 85, 146)
assert os["sampleId"].is_unique and os["patientId"].is_unique
assert os["OS_TIME"].notna().all() and os["OS_TIME"].ge(0).all()
assert os["OS_EVENT"].isin([0, 1]).all()

complete = os[os["complete_case"].eq(True)].copy()
assert (len(complete), int(complete["OS_EVENT"].sum()), int(complete["Altered"].sum())) == (459, 77, 139)
cohorts = [("All 549 source tumors", labels["sampleId"], 255, 161),
           ("498-patient OS cohort", os["sampleId"], 201, 128),
           ("459 complete cases", complete["sampleId"], 190, 122)]
types = ["Missense_Mutation", "Frame_Shift_Del", "Nonsense_Mutation", "Frame_Shift_Ins",
         "Splice_Site", "In_Frame_Del", "In_Frame_Ins"]
assert set(mutations["Mutation Type"]) == set(types)

mutation_rows, annotation_rows, cohort_summaries = [], [], []
for name, samples, expected_records, expected_tumors in cohorts:
    records = mutations[mutations["Sample ID"].isin(set(samples))].copy()
    assert (len(records), records["Sample ID"].nunique()) == (expected_records, expected_tumors)
    type_counts = {}
    for mutation_type in types:
        selected = records[records["Mutation Type"].eq(mutation_type)]
        row = {"cohort": name, "cohort_tumors_n": len(samples), "mutation_type": mutation_type,
               "mutation_records_n": len(selected), "distinct_tumors_n": selected["Sample ID"].nunique()}
        mutation_rows.append(row)
        type_counts[mutation_type] = {"records_n": row["mutation_records_n"],
                                      "distinct_tumors_n": row["distinct_tumors_n"]}
    onco = records["Annotation"].astype("string").str.extract(r"OncoKB:\s*([^,;]+)", expand=False).fillna("No annotation")
    for category, count in onco.value_counts().items():
        annotation_rows.append({"cohort": name, "annotation_source": "Archived OncoKB",
                                "category": category, "mutation_records_n": int(count)})
    functional = records["Functional Impact"].fillna("").str.strip()
    unavailable = functional.eq("MutationAssessor: NA;SIFT: NA;Polyphen-2: NA;AlphaMissense: NA") | functional.eq("")
    for category, count in [("Prediction values unavailable", unavailable.sum()),
                            ("Prediction values present", (~unavailable).sum())]:
        annotation_rows.append({"cohort": name, "annotation_source": "Archived functional predictions",
                                "category": category, "mutation_records_n": int(count)})
    for category, count in records["ClinVar"].fillna("No annotation").value_counts().items():
        annotation_rows.append({"cohort": name, "annotation_source": "Archived ClinVar",
                                "category": category, "mutation_records_n": int(count)})
    if len(samples) == 549:
        assert onco.value_counts().to_dict() == {"Likely Oncogenic": 159, "Oncogenic": 71, "Unknown": 25}
        assert int(unavailable.sum()) == 151
        assert records["Mutation Type"].value_counts().to_dict() == {
            "Missense_Mutation": 104, "Frame_Shift_Del": 54, "Nonsense_Mutation": 53,
            "Frame_Shift_Ins": 19, "Splice_Site": 16, "In_Frame_Del": 8, "In_Frame_Ins": 1}
    cohort_summaries.append({"cohort": name, "cohort_tumors_n": len(samples),
                             "mutation_records_n": len(records), "mutation_bearing_tumors_n": records["Sample ID"].nunique(),
                             "tumors_with_multiple_records_n": int(records["Sample ID"].value_counts().gt(1).sum()),
                             "mutation_type_counts": type_counts,
                             "OncoKB_record_counts": {str(k): int(v) for k, v in onco.value_counts().items()},
                             "functional_prediction_records_unavailable_n": int(unavailable.sum())})

pd.DataFrame(mutation_rows).to_csv(args.out / "supplement_mutation_type_counts.csv", index=False)
pd.DataFrame(annotation_rows).to_csv(args.out / "supplement_annotation_counts.csv", index=False)

# Compute the product-limit steps directly from the existing OS data.
# Censoring removes a patient from future risk sets; it does not count as a death.
km_rows, km_tables = [], []
for value, name in [(0, "PTEN-unaltered"), (1, "PTEN-altered")]:
    group = os[os["Altered"].eq(value)].copy()
    table = group.groupby("OS_TIME")["OS_EVENT"].agg(removed_n="size", deaths_n="sum").reset_index().sort_values("OS_TIME")
    table["censored_n"] = table["removed_n"] - table["deaths_n"]
    table["at_risk_before_time_n"] = len(group) - table["removed_n"].cumsum().shift(fill_value=0)
    table["conditional_survival"] = 1 - table["deaths_n"] / table["at_risk_before_time_n"]
    table["KM_survival"] = table["conditional_survival"].cumprod()
    table["KM_survival_before_time"] = table["KM_survival"].shift(fill_value=1.0)
    table.insert(0, "group", name)
    km_tables.append(table)
    crossing = table[table["KM_survival"].le(0.5)]
    row = {"group": name, "patients_n": len(group), "deaths_n": int(group["OS_EVENT"].sum()),
           "censored_n": int(group["OS_EVENT"].eq(0).sum()),
           "crude_death_fraction": float(group["OS_EVENT"].mean()),
           "maximum_observed_followup_months": float(group["OS_TIME"].max()),
           "last_observed_death_months": float(group.loc[group["OS_EVENT"].eq(1), "OS_TIME"].max()),
           "KM_final_survival": float(table["KM_survival"].iloc[-1]),
           "median_reached": not crossing.empty, "KM_median_months": None,
           "risk_set_at_median_n": None, "deaths_at_median_n": None,
           "KM_survival_before_median": None, "KM_survival_at_median": None,
           "median_is_observed_death_time": None, "extrapolation_used": False}
    if not crossing.empty:
        first = crossing.iloc[0]
        median = float(first["OS_TIME"])
        row.update({"KM_median_months": median, "risk_set_at_median_n": int(first["at_risk_before_time_n"]),
                    "deaths_at_median_n": int(first["deaths_n"]),
                    "KM_survival_before_median": float(first["KM_survival_before_time"]),
                    "KM_survival_at_median": float(first["KM_survival"]),
                    "deaths_before_median_n": int(group.loc[group["OS_TIME"].lt(median), "OS_EVENT"].sum()),
                    "censors_before_median_n": int((group["OS_TIME"].lt(median) & group["OS_EVENT"].eq(0)).sum()),
                    "median_is_observed_death_time": bool((group["OS_TIME"].eq(median) & group["OS_EVENT"].eq(1)).any())})
        assert row["median_is_observed_death_time"] and median <= row["last_observed_death_months"]
    if value == 0:
        assert (row["patients_n"], row["deaths_n"], row["risk_set_at_median_n"]) == (352, 65, 10)
        assert np.isclose(row["KM_median_months"], 110.55)
        assert np.isclose(row["KM_survival_before_median"], 0.5427723592619427)
        assert np.isclose(row["KM_survival_at_median"], 0.4884951233357485)
    else:
        assert (row["patients_n"], row["deaths_n"]) == (146, 20) and not row["median_reached"]
    km_rows.append(row)

pd.DataFrame(km_rows).to_csv(args.out / "supplement_kaplan_meier_audit.csv", index=False)
pd.concat(km_tables, ignore_index=True).to_csv(args.out / "supplement_kaplan_meier_steps.csv", index=False)
summary = {"analysis": "Descriptive mutation and existing OS Kaplan-Meier audits; no Cox models rerun",
           "mutation_source_columns": ["Gene", "Sample ID", "Mutation Type", "Annotation", "Functional Impact", "ClinVar"],
           "OS_source_columns": ["sampleId", "patientId", "Altered", "complete_case", "OS_TIME", "OS_EVENT"],
           "mutation_cohorts": cohort_summaries, "Kaplan_Meier_groups": km_rows,
           "interpretation_notes": [
               "Mutation records and mutation-bearing tumors use different denominators.",
               "A tumor may contribute multiple records or occur in multiple mutation-type categories, so distinct tumor counts by type overlap.",
               "Archived oncogenicity/pathogenicity and prediction annotations do not establish that every mutation causes loss of PTEN function.",
               "The unaltered median is an observed Kaplan-Meier crossing, not an extrapolated value.",
               "Fewer than half the original patients can have observed deaths while Kaplan-Meier survival crosses 0.5 because censoring reduces later risk sets.",
               "Only 10 unaltered patients remained at risk at the median, so the late estimate has limited precision.",
               "Kaplan-Meier estimation relies on the assumption that censoring is noninformative."]}
(args.out / "supplement_descriptive_audits.json").write_text(json.dumps(summary, indent=2) + "\n")
print(pd.DataFrame(mutation_rows).to_string(index=False))
print(pd.DataFrame(km_rows).to_string(index=False))
print("Descriptive audit assertions passed; no Cox models were fitted.")
