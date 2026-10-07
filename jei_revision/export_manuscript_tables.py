#!/usr/bin/env python3
"""Export manuscript Tables 1–3 from the JEI revision analysis outputs.

Usage:
    python jei_revision/export_manuscript_tables.py --input revision_outputs \
        --out revision_outputs

The input root contains ``data_processed/`` and ``results/``. The output root
receives ``tables/`` with formatted CSVs and standalone HTML tables. Values
come from the analysis outputs; this script does not fit models or extrapolate
survival. Decimal half-up rounding preserves the manuscript's displayed values.
Only Python's standard library is required.
"""

import argparse
import csv
from decimal import Decimal, ROUND_HALF_UP
from html import escape
import json
from pathlib import Path


PRIMARY = "Primary portal definition"
FULL_MODEL = "Age, stage and subtype: fixed complete cases"
MODELS = [
    ("Unadjusted: all eligible", "Unadjusted all eligible"),
    ("PTEN only: fixed complete cases", "PTEN only complete cases"),
    ("Age and stage: fixed complete cases", "PTEN plus age and stage"),
    (FULL_MODEL, "PTEN plus age stage and subtype"),
    ("Stage/subtype stratified: fixed complete cases", "Stage and subtype stratified"),
]
VARIABLES = [
    ("Altered", "PTEN status", "Altered vs unaltered"),
    ("AGE_10Y", "Age", "Per 10-year increase"),
    ("ADVANCED_STAGE", "Stage", "Advanced vs early"),
    ("SUBTYPE_UCEC_CN_LOW", "Molecular subtype", "Copy-number-low vs copy-number-high"),
    ("SUBTYPE_UCEC_MSI", "Molecular subtype", "MSI-hypermutated vs copy-number-high"),
    ("SUBTYPE_UCEC_POLE", "Molecular subtype", "POLE-ultramutated vs copy-number-high"),
]
GROUPS = ("Overall", "PTEN-altered", "PTEN-unaltered")


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def select_one(rows, **criteria):
    selected = [r for r in rows if all(r.get(k) == v for k, v in criteria.items())]
    if len(selected) != 1:
        raise ValueError(f"Expected one row for {criteria}, found {len(selected)}")
    return selected[0]


def integer(value):
    number = Decimal(str(value))
    if not number.is_finite() or number != number.to_integral_value():
        raise ValueError(f"Expected a finite integer, got {value!r}")
    return int(number)


def rounded(value, places, trim=False):
    number = Decimal(str(value))
    if not number.is_finite():
        raise ValueError(f"Expected a finite value, got {value!r}")
    result = format(number.quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP), "f")
    if trim and "." in result:
        result = result.rstrip("0").rstrip(".")
    return result


def p_value(value):
    """Three decimals, or two significant figures for very small p values."""
    number = Decimal(str(value))
    if not number.is_finite() or not 0 <= number <= 1:
        raise ValueError(f"Invalid p value {value!r}")
    if number == 0:
        return "<0.001"
    if number >= Decimal("0.001"):
        return rounded(number, 3)
    exponent = number.adjusted()
    coefficient = number.scaleb(-exponent).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP)
    if coefficient == 10:
        coefficient = Decimal("1.0")
        exponent += 1
    superscript = str(exponent).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{coefficient:.1f} × 10{superscript}"


def boolean(value):
    if str(value).lower() in {"true", "1"}:
        return True
    if str(value).lower() in {"false", "0"}:
        return False
    raise ValueError(f"Expected a Boolean, got {value!r}")


def save_csv(path, headers, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)


def html_table(headers, rows, caption, section_rows=()):
    header = "".join(f'<th scope="col">{escape(str(x))}</th>' for x in headers)
    body = []
    for index, row in enumerate(rows):
        if index in section_rows:
            body.append(f'<tr class="section"><th colspan="{len(headers)}" scope="colgroup">'
                        f"{escape(str(row[0]))}</th></tr>")
        else:
            cells = [f'<th scope="row">{escape(str(row[0]))}</th>']
            cells.extend(f"<td>{escape(str(x))}</td>" for x in row[1:])
            body.append("<tr>" + "".join(cells) + "</tr>")
    return (f"<table><caption>{escape(caption)}</caption>"
            f"<thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table>")


def save_html(path, title, tables, notes):
    """Readable black-and-white layout suitable for a Colab HTML display."""
    style = """
.jei-manuscript-table {font-family: Arial, sans-serif; color: #000; background: #fff;
      margin: 24px; line-height: 1.4; max-width: 1200px;}
.jei-manuscript-table h1 {font-size: 20px; margin-bottom: 16px;}
.jei-manuscript-table table {border-collapse: collapse; width: 100%; margin: 16px 0; font-size: 14px;}
.jei-manuscript-table caption {font-weight: bold; text-align: left; margin-bottom: 8px;}
.jei-manuscript-table th, .jei-manuscript-table td {border: 1px solid #999; padding: 7px 9px; text-align: left; vertical-align: top;}
.jei-manuscript-table thead th, .jei-manuscript-table .section th {background: #eee;}
.jei-manuscript-table tbody th {font-weight: normal;}
.jei-manuscript-table p {font-size: 14px; margin: 8px 0;}
@media print {.jei-manuscript-table {margin: 0;} .jei-manuscript-table thead {display: table-header-group;}
              .jei-manuscript-table tr {break-inside: avoid;} }
"""
    paragraphs = "".join(f"<p>{escape(note)}</p>" for note in notes)
    document = ("<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
                '<meta name="viewport" content="width=device-width, initial-scale=1">'
                f"<title>{escape(title)}</title><style>{style}</style></head><body>"
                f'<main class="jei-manuscript-table"><h1>{escape(title)}</h1>'
                f"{''.join(tables)}{paragraphs}</main></body></html>\n")
    path.write_text(document, encoding="utf-8")


def make_table1(input_root, table_root):
    reconstruction = read_csv(input_root / "data_processed/PTEN_status_reconstruction.csv")
    descriptions = read_csv(input_root / "results/OS_age_and_observed_duration_summary.csv")
    medians = read_csv(input_root / "results/OS_Kaplan_Meier_summary.csv")
    counts = read_csv(input_root / "results/cohort_counts.csv")
    if len({r["sampleId"] for r in reconstruction}) != len(reconstruction):
        raise ValueError("PTEN reconstruction contains duplicate sample IDs")
    mutation = [boolean(r["MUTATION_PRESENT"]) for r in reconstruction]
    high_cna = [r["PTEN_CNA"] != "" and Decimal(r["PTEN_CNA"]) in {-2, 2}
                for r in reconstruction]
    altered = [integer(r["RECONSTRUCTED_ALTERED"]) for r in reconstruction]
    if any(a not in {0, 1} for a in altered):
        raise ValueError("PTEN reconstruction must use 0/1 status")
    if any(a != int(m or c) for a, m, c in zip(altered, mutation, high_cna)):
        raise ValueError("Reconstructed PTEN status disagrees with mutation/CNA rule")
    matches = sum(a == integer(r["Altered"]) for a, r in zip(altered, reconstruction))
    rows_a = [
        ["Samples in PTEN matrix", str(len(reconstruction))],
        ["Reconstructed PTEN-altered samples", str(sum(altered))],
        ["Reconstructed PTEN-unaltered samples", str(len(altered) - sum(altered))],
        ["Mutation only", str(sum(m and not c for m, c in zip(mutation, high_cna)))],
        ["High-level copy-number alteration only", str(sum(c and not m for m, c in zip(mutation, high_cna)))],
        ["Mutation and high-level copy-number alteration", str(sum(m and c for m, c in zip(mutation, high_cna)))],
        ["Exact matches to supplied portal labels", str(matches)],
        ["Label mismatches", str(len(reconstruction) - matches)],
    ]
    # Derive the named cohort from the endpoint audit instead of fixing a count.
    audit = select_one(counts, endpoint="OS", cohort="Endpoint-eligible")
    cohort_name = f"{integer(audit['patients_n'])}-patient OS cohort"
    groups = {name: select_one(descriptions, cohort=cohort_name, group=name) for name in GROUPS}
    km = {name: select_one(medians, cohort=cohort_name, group=name) for name in GROUPS}
    for name in GROUPS:
        row = groups[name]
        if integer(row["patients_n"]) != integer(row["deaths_n"]) + integer(row["censored_n"]):
            raise ValueError(f"Inconsistent survival counts for {name}")
        if integer(row["patients_n"]) != integer(row["age_years_available_n"]) + integer(row["age_years_missing_or_nonnumeric_n"]):
            raise ValueError(f"Inconsistent age availability for {name}")
        if boolean(km[name]["extrapolation_used"]):
            raise ValueError("Table 1 requires observed Kaplan–Meier estimates without extrapolation")
    for field, name in [("patients_n", "Overall"), ("altered_n", "PTEN-altered"), ("unaltered_n", "PTEN-unaltered")]:
        if integer(audit[field]) != integer(groups[name]["patients_n"]):
            raise ValueError(f"OS cohort audit disagrees for {name}")
    rows_b = []
    rows_b.append(["Analyzable patients"] + [str(integer(groups[n]["patients_n"])) for n in GROUPS])
    rows_b.append(["Age, median (IQR), years"] + [
        f"{rounded(groups[n]['age_years_median'], 1, trim=True)} "
        f"({rounded(groups[n]['age_years_q1'], 1, trim=True)}–"
        f"{rounded(groups[n]['age_years_q3'], 1, trim=True)})" for n in GROUPS])
    rows_b.append(["Age available / missing"] + [
        f"{integer(groups[n]['age_years_available_n'])} / "
        f"{integer(groups[n]['age_years_missing_or_nonnumeric_n'])}" for n in GROUPS])
    rows_b.append(["Deaths"] + [str(integer(groups[n]["deaths_n"])) for n in GROUPS])
    rows_b.append(["Censored"] + [str(integer(groups[n]["censored_n"])) for n in GROUPS])
    rows_b.append(["Median overall survival (months)", ""] + [
        rounded(km[n]["KM_median_months"], 1) if km[n]["KM_median_months"] else "Not reached"
        for n in GROUPS[1:]])
    headers_b = ["Measure", *GROUPS]
    combined = [["A. PTEN status reconstruction", "", "", ""]]
    combined.extend([r[0], r[1], "", ""] for r in rows_a)
    combined.append(["B. Patient-level survival cohort", "", "", ""])
    combined.extend(rows_b)
    save_csv(table_root / "Table1.csv", headers_b, combined)
    save_csv(table_root / "Table1A_PTEN_reconstruction.csv", ["Measure", "Count"], rows_a)
    save_csv(table_root / "Table1B_survival_cohort.csv", headers_b, rows_b)
    notes = [
        "Tumors were classified as PTEN-altered if they had a reported PTEN mutation, a GISTIC −2 deep deletion, or a GISTIC +2 amplification. High-level copy-number alteration includes GISTIC −2 or +2.",
        "Censored counts equal patients minus deaths. Median survival is the first recorded time at which the Kaplan–Meier estimate reaches 50% or below. No extrapolation was used. Censoring is accounted for in the risk sets, so the Kaplan–Meier estimate can reach 50% even when fewer than half of the original patients have recorded deaths.",
        "Age is reported as the median and interquartile range (IQR, 25th–75th percentiles) among patients with available age data. Available and missing age counts are shown separately. The overall median-survival cell is left blank to match the manuscript table.",
    ]
    unaltered = km["PTEN-unaltered"]
    if unaltered["KM_median_months"]:
        notes.insert(2, f"The unaltered median depended on {integer(unaltered['risk_set_at_median_n'])} patients still at risk.")
    if not km["PTEN-altered"]["KM_median_months"]:
        notes.insert(2, "The altered curve remained above 50% throughout the recorded follow-up, so its median was not reached.")
    title = "Table 1. PTEN reconstruction and survival-cohort characteristics"
    save_html(table_root / "Table1.html", title, [
        html_table(["Measure", "Count"], rows_a, "A. PTEN status reconstruction"),
        html_table(headers_b, rows_b, "B. Patient-level survival cohort")], notes)


def make_table2(input_root, table_root, models, ph_tests):
    coefficients = read_csv(input_root / "results/all_cox_coefficients.csv")
    selected = [r for r in coefficients if r["endpoint"] == "OS" and
                r["definition"] == PRIMARY and r["model"] == FULL_MODEL]
    if set(r["variable"] for r in selected) != set(v[0] for v in VARIABLES) or len(selected) != len(VARIABLES):
        raise ValueError("Fully adjusted OS model must have the six manuscript covariates")
    rows = []
    for code, variable, comparison in VARIABLES:
        row = select_one(selected, variable=code)
        rows.append([variable, comparison, rounded(row["exp(coef)"], 2),
                     rounded(row["exp(coef) lower 95%"], 2),
                     rounded(row["exp(coef) upper 95%"], 2), p_value(row["p"])])
    headers = ["Variable", "Comparison", "Hazard ratio", "95% CI lower", "95% CI upper", "p value"]
    model = select_one(models, endpoint="OS", definition=PRIMARY, model=FULL_MODEL)
    notes = [
        f"The model included {integer(model['patients_n'])} patients and {integer(model['events_n'])} deaths. Reference groups were PTEN-unaltered, early stage, and copy-number-high. Age was measured per 10-year increase.",
        "Hazard ratios above 1 indicate a higher estimated rate of death relative to the stated reference group; values below 1 indicate a lower estimated rate. HR, hazard ratio; CI, confidence interval; MSI, microsatellite instability.",
    ]
    relevant_ph = [r for r in ph_tests if r["endpoint"] == "OS" and
                   r["definition"] == PRIMARY and r["model"] == FULL_MODEL]
    if set(r["variable"] for r in relevant_ph) != set(v[0] for v in VARIABLES) or len(relevant_ph) != len(VARIABLES):
        raise ValueError("Fully adjusted OS model lacks complete proportional-hazards diagnostics")
    violations = [r for r in relevant_ph if Decimal(r["p"]) < Decimal("0.05")]
    labels = {v[0]: v[2] for v in VARIABLES}
    if violations:
        terms = "; ".join(labels[r["variable"]] for r in violations)
        notes.append(f"Ranked-time Schoenfeld tests detected proportional-hazards violations for: {terms}. The stage- and subtype-stratified sensitivity model is reported in the Results.")
    else:
        notes.append("No proportional-hazards violations were detected in this model using ranked-time Schoenfeld tests.")
    title = "Table 2. Fully adjusted Cox proportional-hazards model"
    save_csv(table_root / "Table2.csv", headers, rows)
    save_html(table_root / "Table2.html", title, [html_table(headers, rows, title)], notes)


def make_table3(table_root, models, ph_tests):
    rows = []
    selected_models = []
    for endpoint in ("PFS", "PFI"):
        for model_name, display_name in MODELS:
            row = select_one(models, endpoint=endpoint, definition=PRIMARY, model=model_name)
            select_one(ph_tests, endpoint=endpoint, definition=PRIMARY,
                       model=model_name, variable="Altered")
            selected_models.append(row)
            rows.append([endpoint, display_name, str(integer(row["patients_n"])),
                         str(integer(row["events_n"])), rounded(row["PTEN_HR"], 2),
                         f"{rounded(row['CI_95_low'], 2)}–{rounded(row['CI_95_high'], 2)}",
                         p_value(row["PTEN_p_value"])])
    complete_counts = {integer(r["patients_n"]) for r in selected_models if r["model"] != MODELS[0][0]}
    all_counts = {integer(r["patients_n"]) for r in selected_models if r["model"] == MODELS[0][0]}
    if len(complete_counts) != 1 or len(all_counts) != 1:
        raise ValueError("Manuscript progression models must use fixed common patient counts")
    headers = ["Endpoint", "Model", "Patients", "Events", "HR", "95% CI", "p value"]
    notes = [
        f"The sequential and stratified models used the same {next(iter(complete_counts))} complete cases; the all-eligible unadjusted models used {next(iter(all_counts))} patients. PTEN-unaltered tumors were the reference.",
        "Adjusted models used age per 10 years, early stage as the stage reference, and copy-number-high as the subtype reference. Stratified models estimated PTEN and age effects while allowing baseline hazards to differ by stage and subtype.",
        "PFS events included progression, recurrence, distant metastasis, a new primary tumor at any site, or death from any cause. PFI included the same tumor events or death with tumor; deaths without tumor were censored. PFS therefore does not remove all non-cancer deaths, and death with tumor does not establish the cause of death.",
        "All secondary p values are nominal and were not adjusted for multiple comparisons. CI, confidence interval; PFS, progression-free survival; PFI, progression-free interval.",
    ]
    violations = [r for r in ph_tests if r["endpoint"] in {"PFS", "PFI"}
                  and r["definition"] == PRIMARY and Decimal(r["p"]) < Decimal("0.05")]
    intermediate = {MODELS[1][0], MODELS[2][0]}
    expected_intermediate = {(ep, model, "Altered") for ep in ("PFS", "PFI") for model in intermediate}
    actual = {(r["endpoint"], r["model"], r["variable"]) for r in violations}
    if actual == expected_intermediate:
        notes.insert(3, "PTEN violated the proportional-hazards assumption in the complete-case PTEN-only and age-and-stage models for both endpoints. These intermediate HRs summarize the observed association over follow-up and should not be interpreted as constant effects. No violations were detected in the fully adjusted or stratified models.")
    elif violations:
        display_models = dict(MODELS)
        terms = "; ".join(f"{r['endpoint']}, {display_models.get(r['model'], r['model'])}, {r['variable']}"
                          for r in violations)
        notes.insert(3, f"Ranked-time Schoenfeld tests detected proportional-hazards violations for: {terms}. Affected HRs should be interpreted as summary estimates rather than constant effects over follow-up.")
    else:
        notes.insert(3, "No proportional-hazards violations were detected in these models using ranked-time Schoenfeld tests.")
    title = "Table 3. PTEN hazard ratios for secondary progression endpoints"
    save_csv(table_root / "Table3.csv", headers, rows)
    save_html(table_root / "Table3.html", title, [html_table(headers, rows, title)], notes)


def export_tables(input_root, output_root):
    input_root = Path(input_root).resolve()
    table_root = Path(output_root).resolve() / "tables"
    table_root.mkdir(parents=True, exist_ok=True)
    models = read_csv(input_root / "results/model_comparison.csv")
    ph_tests = read_csv(input_root / "results/proportional_hazards_rank_tests.csv")
    make_table1(input_root, table_root)
    make_table2(input_root, table_root, models, ph_tests)
    make_table3(table_root, models, ph_tests)
    return [table_root / name for name in [
        "Table1.csv", "Table1A_PTEN_reconstruction.csv", "Table1B_survival_cohort.csv",
        "Table1.html", "Table2.csv", "Table2.html", "Table3.csv", "Table3.html"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="Analysis output root containing data_processed/ and results/")
    parser.add_argument("--out", type=Path,
                        help="Output root; tables/ is created here (default: --input)")
    args = parser.parse_args()
    paths = export_tables(args.input, args.out or args.input)
    print(json.dumps({"tables": [str(path) for path in paths]}, indent=2))


if __name__ == "__main__":
    main()
