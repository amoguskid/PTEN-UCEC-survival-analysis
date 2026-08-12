# PTEN/UCEC Reproducible Survival Analysis
## Purpose

This package contains four linked analyses:

1. **Phase 1 — Independent survival rerun:** Rebuilds the PTEN altered-versus-unaltered overall-survival analysis from the supplied Firehose Legacy exports.
2. **Phase 2 — New molecular-subtype analysis:** Tests whether PTEN alteration frequency differs across the four TCGA-UCEC molecular subtypes using supplied PanCancer Atlas clinical files.
3. **Phase 3 — Adjusted survival analysis:** Fits PTEN-only, age/stage-adjusted, and age/stage/subtype-adjusted Cox models on the same complete-case cohort, with proportional-hazards diagnostics and a stratified sensitivity model.
4. **Phase 4 — Direct PTEN reconstruction and final figures:** Reconstructs PTEN status from the supplied mutation and discrete CNA tables and creates three final figures.

The package is designed so another person can keep only the raw files and master notebook, run all cells, and recreate the processed datasets, exclusion records, statistics, and figures.

## Research questions

### Survival rerun

Is overall survival different between PTEN-altered and PTEN-unaltered tumors in the supplied exported patient-level data?

### Molecular-subtype extension

Is PTEN alteration frequency associated with TCGA molecular subtype in uterine corpus endometrial carcinoma?

### Adjusted survival extension

Does the observed association between PTEN alteration status and overall survival remain after accounting for age, clinical stage, and TCGA molecular subtype?

### Direct PTEN reconstruction

Can the supplied portal-defined PTEN altered/unaltered label be reproduced exactly from any PTEN mutation or a high-level discrete CNA (GISTIC -2 or +2)?

## Data sources

### Firehose Legacy exports

- PTEN altered/unaltered sample matrix
- Patient clinical and survival table
- Portal survival summary
- Supporting mutation and copy-number tables
- Manuscript access date: January 19, 2026

### PanCancer Atlas clinical files

- `data_clinical_patient.txt`
- `data_clinical_sample.txt`
- Study: Uterine Corpus Endometrial Carcinoma (TCGA, PanCancer Atlas)
- Study ID: `ucec_tcga_pan_can_atlas_2018`
- Files supplied for this project on July 14, 2026

The raw TCGA-derived data files are not redistributed in this public repository. The source studies are publicly available through cBioPortal: [Uterine Corpus Endometrial Carcinoma, Firehose Legacy](https://www.cbioportal.org/study/summary?id=ucec_tcga) (study ID: `ucec_tcga`) and [Uterine Corpus Endometrial Carcinoma, PanCancer Atlas](https://www.cbioportal.org/study/summary?id=ucec_tcga_pan_can_atlas_2018) (study ID: `ucec_tcga_pan_can_atlas_2018`).

## Exact input-file acquisition and construction

Create `data_raw/` and put the seven files below in it. cBioPortal's **Download** tab
provides the study archive; its **Query** results provide the table-download menus.
Use all samples in the indicated study and query the single gene `PTEN`. Do not use
the PanCancer Atlas release for S1--S5: release mixing is limited to the explicitly
documented subtype join.

| Required filename | Study ID | Direct or derived | Exact construction and required columns |
|---|---|---|---|
| `Supplementary_Table_S1_PTEN_altered_unaltered_sample_matrix.tsv` | `ucec_tcga` | Derived from the PTEN query's downloaded sample-by-gene/OncoPrint data | Retain one row for each of the 549 queried samples. `studyID:sampleId` must be `ucec_tcga:<SAMPLE_ID>`; `PTEN` is the downloaded PTEN alteration annotation; set integer `Altered` to 1 when the query reports any PTEN mutation or discrete CNA of -2 or +2, otherwise 0. Required columns: `studyID:sampleId`, `Altered`, `PTEN`. |
| `Supplementary_Table_S2_PTEN_discrete_CNA_table.tsv` | `ucec_tcga` | Derived by transposing/filtering the study archive's discrete GISTIC CNA file to PTEN | Select the `PTEN` gene row, transpose it to one row per queried sample, and prepend `ucec_tcga:` to each sample ID. Keep GISTIC values as integers (-2, -1, 0, 1, 2). Required columns: `studyID:sampleId`, `PTEN`; include the same 549 samples and use the required filename. |
| `Supplementary_Table_S3_PTEN_mutation_table.tsv` | `ucec_tcga` | Direct PTEN-only mutation-table export from Query Results → Mutations (a filtered copy of `data_mutations_extended.txt` is equivalent) | Export/filter records where `Hugo_Symbol` is `PTEN`; do not deduplicate distinct variants. Required identifier is `studyID:sampleId`, `SAMPLE_ID`, or `Tumor_Sample_Barcode`; if building from the archive, rename/copy `Tumor_Sample_Barcode` to `studyID:sampleId` and prefix it with `ucec_tcga:`. Retain `Hugo_Symbol`; an optional `PTEN` annotation column is accepted. |
| `Supplementary_Table_S4_cBioPortal_clinical_survival_table.tsv` | `ucec_tcga` | Derived patient-level table from the Query Results clinical-data export | Select one representative tumor sample per patient and retain the portal patient/sample IDs and clinical fields. Required columns: `patientId`, `sampleId`, `OS_MONTHS`, `OS_STATUS`; the notebook also carries through `AGE`, `CLINICAL_STAGE`, `GRADE`, `HISTOLOGICAL_DIAGNOSIS`, and `RACE`, so retain those columns (blank values are allowed). Preserve the 500 exported patient rows and the exact required filename. |
| `Supplementary_Table_S5_cBioPortal_PTEN_survival_summary_table.tsv` | `ucec_tcga` | Derived by downloading/transcribing the PTEN query's Comparison/Survival summary, not patient-level data | Make one row per displayed survival endpoint, including the `overall` row. Required columns (including spelling and punctuation): `Survival Type`, `Number of Patients`, `# in Altered group`, `# in Unaltered group`, `p-Value`. Values must be the portal results for that saved all-sample PTEN query. |
| `data_clinical_patient.txt` | `ucec_tcga_pan_can_atlas_2018` | Direct, unedited file from the PanCancer Atlas study archive | Required columns: `PATIENT_ID`, `SUBTYPE`, `CANCER_TYPE_ACRONYM`, `IN_PANCANPATHWAYS_FREEZE`. Keep cBioPortal's leading `#` metadata lines and filename. |
| `data_clinical_sample.txt` | `ucec_tcga_pan_can_atlas_2018` | Direct, unedited file from the PanCancer Atlas study archive | Required columns: `PATIENT_ID`, `SAMPLE_ID`, `SAMPLE_TYPE`. Keep cBioPortal's leading `#` metadata lines and filename. The supplied selection contains 529 unique primary samples. |

Before running, check that tab delimiters and headers have not been changed by a
spreadsheet program. S1, S2, S4, and S5 are deliberately named derived
supplementary tables: they are not files with those names in a cBioPortal archive.
The two `data_clinical_*.txt` files are direct archive files. Phase 4 validates the
expected reconstruction counts and stops rather than silently continuing if any
input was built from a different query or release.


## Why the subtype analysis uses 507 tumors rather than the 498-patient survival cohort

The subtype-frequency question does not require survival data. Restricting it to patients with complete survival information would unnecessarily discard eligible tumors and could introduce selection bias.

The Phase 2 cohort is therefore built from the exact 529 primary samples selected in the PanCancer Atlas sample file:

- 529 selected primary samples
- 529 exact matches to the Firehose PTEN matrix
- 22 missing molecular-subtype labels
- 507 tumors included in the subtype association analysis

The Firehose PTEN matrix contains 549 sample rows. Twenty sample rows were not selected in the PanCancer Atlas sample file and are documented in the exclusion audit.

## Key reproducible results

### Phase 1

- 498 analyzable patients
- 146 PTEN-altered
- 352 PTEN-unaltered
- 85 deaths
- Log-rank p ≈ 0.0398
- Unadjusted hazard ratio ≈ 0.593
- 95% CI ≈ 0.359–0.981

This supports the same general survival association as the portal analysis, but it is not an exact replication of the portal's 543-patient cohort.

### Phase 2

- 507 tumors with PTEN status and molecular subtype
- Pearson chi-square ≈ 53.27
- 3 degrees of freedom
- p ≈ 1.60 × 10⁻¹¹
- Cramér's V ≈ 0.324

Observed PTEN-altered frequencies:

| Molecular subtype | Altered / total | Percent |
|---|---:|---:|
| POLE-ultramutated | 16 / 49 | 32.7% |
| MSI-hypermutated | 59 / 148 | 39.9% |
| Copy-number low | 76 / 147 | 51.7% |
| Copy-number high | 22 / 163 | 13.5% |

These results show an unadjusted association between PTEN alteration frequency and molecular subtype. They do not establish causation or independent prognostic value.

### Phase 3

All three primary Cox models used the same 459-patient complete-case cohort with 77 deaths:

| Model | PTEN HR | 95% CI | p-value |
|---|---:|---:|---:|
| PTEN only | 0.595 | 0.354–1.001 | 0.0504 |
| PTEN + age + stage | 0.826 | 0.483–1.411 | 0.4840 |
| PTEN + age + stage + subtype | 1.071 | 0.589–1.948 | 0.8214 |

The unadjusted association attenuated after adjustment. A stage- and subtype-stratified sensitivity model gave a similar null PTEN result (HR 1.139, 95% CI 0.623–2.080; p=0.6732). This supports the cautious conclusion that PTEN alteration was not independently associated with overall survival in this exported cohort.

### Phase 4

- 549 samples were reconstructed from the mutation and CNA tables.
- 161 samples had at least one PTEN mutation.
- 26 samples had a high-level PTEN CNA (-2 or +2); six also had a mutation.
- The reconstructed label matched the supplied portal label for 549/549 samples.

The three recommended figures are:

- `figures/final_01_kaplan_meier.png`
- `figures/final_02_subtype_frequency.png`
- `figures/final_03_adjusted_cox_forest.png`

Vector PDF versions are saved beside the PNG files.

## Folder structure

PTEN-UCEC-survival-analysis/
├── README.md
├── requirements.txt
├── PTEN_UCEC_OCSEF_master_notebook.ipynb
├── PTEN_amplification_excluded_sensitivity.py
├── VALIDATION_REPORT.md
├── ANALYSIS_NOTES.md
└── DATA_DICTIONARY.md

## How to run

1. Install Python 3.11 or newer.
2. Download or clone this repository and open a terminal in the project folder.
3. Create a folder named `data_raw/`.
4. Obtain the public TCGA-derived source files from the cBioPortal studies identified above and place the required files in `data_raw/` using these filenames:

   - `Supplementary_Table_S1_PTEN_altered_unaltered_sample_matrix.tsv`
   - `Supplementary_Table_S2_PTEN_discrete_CNA_table.tsv`
   - `Supplementary_Table_S3_PTEN_mutation_table.tsv`
   - `Supplementary_Table_S4_cBioPortal_clinical_survival_table.tsv`
   - `Supplementary_Table_S5_cBioPortal_PTEN_survival_summary_table.tsv`
   - `data_clinical_patient.txt`
   - `data_clinical_sample.txt`

5. Install the required packages:

   ```bash
   pip install -r requirements.txt
   ```

6. Open `PTEN_UCEC_OCSEF_master_notebook.ipynb`.
7. Select **Restart Kernel and Run All Cells**.
8. Confirm that all 20 code cells finish without errors and that the `data_processed/`, `results/`, and `figures/` folders are created.
9. After the notebook finishes, run the amplification-excluded sensitivity analysis:

   ```bash
   python PTEN_amplification_excluded_sensitivity.py
   ```

10. Review `results/PTEN_status_reconstruction_validation.csv`, `results/PTEN_adjustment_comparison.csv`, `results/final_key_findings.csv`, and the amplification-excluded sensitivity output files.
## Interpretation of findings

> I exactly reconstructed the supplied PTEN label from mutation and high-level copy-number data. My independent survival rerun supported the portal's general unadjusted association but did not reproduce its exact cohort. PTEN alteration frequency differed substantially across molecular subtypes, and the survival association disappeared after adjustment for age, stage, and subtype. Therefore, this dataset does not support PTEN alteration as an independent survival predictor, although PTEN remains biologically important in UCEC.

## Important limitations

- Firehose Legacy PTEN status and PanCancer Atlas subtype labels come from different TCGA study releases.
- The datasets were harmonized by exact TCGA sample and patient identifiers, but release-specific processing differences may remain.
- Twenty-two PanCancer Atlas patients lacked subtype labels.
- The chi-square analysis is unadjusted.
- The subtype analysis tests alteration frequency, not subtype-specific survival or treatment response.
- Portal-defined PTEN status combines mutation and discrete copy-number alteration information.
- The adjusted Cox analysis used 459 complete cases and 77 deaths, limiting precision.
- Stage and the MSI indicator showed evidence of non-proportional hazards in the fully adjusted model; the stage/subtype-stratified sensitivity analysis preserved the null PTEN result.
- The adjusted model did not include every possible clinical factor, such as grade, histology, treatment, or comorbidity.
- No external cohort was used for validation.
