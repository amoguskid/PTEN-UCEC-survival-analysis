# Data Dictionary

## Phase 4 raw inputs

`Supplementary_Table_S2_PTEN_discrete_CNA_table.tsv` uses `SAMPLE_ID` for the
sample identifier and `PTEN` for the discrete CNA value. Ten `NP` calls remain
missing and are not recoded. The supplied
`Supplementary_Table_S3_PTEN_mutation_table.tsv` uses `Gene` and `Sample ID`.
The finalization script also accepts the documented cBioPortal alternatives
`sampleId` or `Tumor_Sample_Barcode` for sample identifiers and `Hugo_Symbol`
or `HUGO_SYMBOL` for the gene name.

## `data_processed/subtype_analysis_dataset.csv`

| Column | Meaning |
|---|---|
| `PATIENT_ID` | TCGA patient identifier |
| `SAMPLE_ID` | Exact PanCancer Atlas selected primary sample identifier |
| `SUBTYPE` | Original PanCancer Atlas subtype label |
| `SUBTYPE_READABLE` | Readable subtype label used in figures |
| `IN_PANCANPATHWAYS_FREEZE` | Whether the patient was included in the PanCancer pathway-analysis freeze |
| `Altered` | PTEN status coded 1 for altered and 0 for unaltered |
| `PTEN` | Original PTEN alteration annotation from the supplied matrix |
| `PTEN_GROUP` | Readable altered/unaltered group |

## `results/subtype_statistics.csv`

One row per molecular subtype, including total tumors, altered and unaltered counts, percentage altered, and Wilson 95% confidence limits.

## `results/subtype_chi_square_test.csv`

The overall 4 × 2 Pearson chi-square result, degrees of freedom, p-value, Cramér's V, and expected-count check.

## `results/subtype_pairwise_fisher_holm.csv`

Exploratory pairwise subtype comparisons using Fisher's exact test with Holm correction for six comparisons.

## `data_processed/subtype_excluded_samples.csv`

Documents:
- Firehose PTEN samples not selected by the PanCancer Atlas clinical sample file; and
- selected PanCancer Atlas samples excluded because subtype was missing.

## `data_processed/adjusted_cox_complete_case_dataset.csv`

The fixed 459-patient cohort used for all three primary Cox models. It contains overall-survival time and event coding, PTEN status, centered age per 10 years, early/advanced stage, and molecular subtype.

## `results/PTEN_adjustment_comparison.csv`

One row for the PTEN coefficient from each primary model:

- PTEN only;
- PTEN + age + stage; and
- PTEN + age + stage + molecular subtype.

Columns include cohort size, deaths, concordance index, hazard ratio, 95% confidence limits, and p-value.

## `results/all_adjusted_cox_PH_tests.csv`

Rank-transform proportional-hazards diagnostics for every covariate in the three primary Cox models. `violates_PH_at_0.05` marks diagnostic p-values below 0.05.

## `results/stage_subtype_stratified_sensitivity_coefficients.csv`

Sensitivity model that stratifies baseline hazards by stage and subtype while estimating PTEN and age effects.

## `data_processed/PTEN_status_reconstruction.csv`

One row per sample in the supplied PTEN matrix. Main fields:

| Column | Meaning |
|---|---|
| `SAMPLE_ID` | TCGA sample identifier |
| `PORTAL_ALTERED` | Supplied PTEN altered/unaltered label |
| `PTEN_CNA` | Discrete PTEN CNA call |
| `MUTATION_PRESENT` | At least one PTEN mutation appears in the supplied mutation table |
| `HIGH_LEVEL_CNA_PRESENT` | PTEN CNA equals -2 or +2 |
| `RECONSTRUCTED_ALTERED` | Mutation or high-level CNA rule |
| `ALTERATION_SOURCE` | Mutation only, CNA only, both, or none |
| `MATCHES_PORTAL_LABEL` | Sample-level agreement check |

## `results/PTEN_status_reconstruction_summary.csv`

Counts for mutations, high-level CNAs, overlap, reconstructed altered samples, portal-labeled altered samples, exact matches, and mismatches.

## `results/PTEN_status_reconstruction_validation.csv`

Expected and observed values for the six reconstruction checks, with a Boolean
`matches_expected` field. Phase 4 stops before figure generation if any check
fails.

## `results/final_key_findings.csv`

Four-row summary used for the project narrative and poster: direct reconstruction, independent survival rerun, subtype association, and fully adjusted Cox result.

## `results/clean_run_final_checksum_comparison.csv`

SHA-256 comparison for all generated CSV outputs and the three final poster PNGs across two fresh executions.
