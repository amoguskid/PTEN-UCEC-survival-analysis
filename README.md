# PTEN/UCEC Reproducible Survival Analysis

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/amoguskid/PTEN-UCEC-survival-analysis/blob/main/PTEN_UCEC_OCSEF_master_notebook.ipynb)

## Purpose

This package contains four original linked analyses and the integrated JEI manuscript revision:

1. **Phase 1 — Independent survival rerun:** Rebuilds the PTEN altered-versus-unaltered overall-survival analysis from the supplied Firehose Legacy exports.
2. **Phase 2 — New molecular-subtype analysis:** Tests whether PTEN alteration frequency differs across the four TCGA-UCEC molecular subtypes using supplied PanCancer Atlas clinical files.
3. **Phase 3 — Adjusted survival analysis:** Fits PTEN-only, age/stage-adjusted, and age/stage/subtype-adjusted Cox models on the same complete-case cohort, with proportional-hazards diagnostics and a stratified sensitivity model.
4. **Phase 4 — Direct PTEN reconstruction and final figures:** Reconstructs PTEN status from the supplied mutation and discrete CNA tables and creates three final figures.

**JEI revision — notebook sections 20–25:** Adds the reviewer-requested secondary PFS and PFI analyses, the revised manuscript Figure 3, new Figure 4, and Tables 1–3. Overall survival remains the primary endpoint. The notebook contains 54 cells, including 27 code cells, and includes saved figure and table displays.

The repository code and original seven source exports are sufficient to regenerate the processed datasets, exclusion records, statistics, and figures. The setup cell obtains the code in Colab, installs the pinned dependencies, accepts the source-file upload, and downloads the TCGA-CDR workbook described below.

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

### TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR)

The JEI revision uses the TCGA-CDR workbook for the secondary progression endpoints. The setup cell automatically downloads `TCGA-CDR-SupplementalTableS1.xlsx` from the official NCI/GDC source linked on the [TCGA Pan-Cancer Clinical Data Resource publication page](https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018): [workbook download](https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81). It saves the workbook under `reference_inputs/` and checks this SHA-256 digest before analysis:

```text
ea594c0fbb6731477c7ac511fab449ca9c38b0d42d269591ed9f5c4090e75a5a
```

Source: Liu et al. (2018), DOI [10.1016/j.cell.2018.02.052](https://doi.org/10.1016/j.cell.2018.02.052). The analysis also records input checksums and the actual software versions used in each run.

The raw TCGA-derived data files are not redistributed in this public repository. The source studies are publicly available through cBioPortal: [Uterine Corpus Endometrial Carcinoma, Firehose Legacy](https://www.cbioportal.org/study/summary?id=ucec_tcga) (study ID: `ucec_tcga`) and [Uterine Corpus Endometrial Carcinoma, PanCancer Atlas](https://www.cbioportal.org/study/summary?id=ucec_tcga_pan_can_atlas_2018) (study ID: `ucec_tcga_pan_can_atlas_2018`).

### Required input files

Create `data_raw/` and place the seven files below in it. The S1–S5 filenames
refer to project-specific supplementary tables, not files distributed under those
names in a cBioPortal study archive. Use the Firehose Legacy release for S1–S5;
the PanCancer Atlas release is used only for the documented molecular-subtype
join.

| Required filename | Study | Source and required columns |
|---|---|---|
| `Supplementary_Table_S1_PTEN_altered_unaltered_sample_matrix.tsv` | `ucec_tcga` | Derived from the all-sample PTEN query. One row per sample (549 rows) with `studyID:sampleId`, the downloaded `PTEN` annotation, and binary `Altered` (mutation or GISTIC -2/+2). |
| `Supplementary_Table_S2_PTEN_discrete_CNA_table.tsv` | `ucec_tcga` | PTEN row from the discrete GISTIC CNA data, transposed to one row per sample. Required columns: `SAMPLE_ID`, `PTEN`. Preserve `NP` where a CNA call is unavailable. |
| `Supplementary_Table_S3_PTEN_mutation_table.tsv` | `ucec_tcga` | PTEN mutation records from the mutation export. Required columns in the supplied file: `Gene`, `Sample ID`. Retain separate variant records. |
| `Supplementary_Table_S4_cBioPortal_clinical_survival_table.tsv` | `ucec_tcga` | Patient-level clinical export with one representative tumor sample per patient. Required columns: `patientId`, `sampleId`, `OS_MONTHS`, `OS_STATUS`; retain `AGE`, `CLINICAL_STAGE`, `GRADE`, `HISTOLOGICAL_DIAGNOSIS`, and `RACE` when available. The supplied table has 500 rows. |
| `Supplementary_Table_S5_cBioPortal_PTEN_survival_summary_table.tsv` | `ucec_tcga` | Survival summary from the same saved all-sample PTEN query. Required columns: `Survival Type`, `Number of Patients`, `# in Altered group`, `# in Unaltered group`, `p-Value`. |
| `data_clinical_patient.txt` | `ucec_tcga_pan_can_atlas_2018` | Unedited study-archive file. Required fields: `PATIENT_ID`, `SUBTYPE`, `CANCER_TYPE_ACRONYM`, `IN_PANCANPATHWAYS_FREEZE`. Keep the leading cBioPortal metadata lines. |
| `data_clinical_sample.txt` | `ucec_tcga_pan_can_atlas_2018` | Unedited study-archive file. Required fields: `PATIENT_ID`, `SAMPLE_ID`, `SAMPLE_TYPE`. Keep the leading cBioPortal metadata lines. The supplied selection contains 529 unique primary samples. |

Keep the files tab-delimited and do not allow spreadsheet software to change the
headers or TCGA identifiers. The notebook validates the expected cohort counts,
and Phase 4 independently validates all six PTEN reconstruction totals before it
writes the final outputs.


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

The three original Phase 4 figures are:

- `figures/final_01_kaplan_meier.png`
- `figures/final_02_subtype_frequency.png`
- `figures/final_03_adjusted_cox_forest.png`

Vector PDF versions are saved beside the PNG files.

### JEI revision outputs

Sections 20–25 run the revision scripts and display these manuscript outputs:

| Section | Output | Contents |
|---|---|---|
| 20 | Revision analyses | OS, PFS, and PFI models, diagnostics, cohort and descriptive audits, and the mutation-or-deletion sensitivity definition. |
| 21 | `figures/Figure3_Cox_models.png` and `.pdf` | Revised manuscript Figure 3: the complete-case PTEN-only estimate and all six fully adjusted OS coefficients, with numerical HRs, confidence intervals, and p values. |
| 22 | `figures/Figure4_PFS_PFI_Kaplan_Meier.png` and `.pdf` | PFS and PFI Kaplan–Meier panels with confidence bands, censoring marks, and numbers at risk. |
| 23 | `tables/Table1.csv`, `.html`; `tables/Table2.csv`, `.html` | PTEN reconstruction and survival-cohort characteristics, including age and observed Kaplan–Meier medians; all six fully adjusted OS coefficients. |
| 24 | `tables/Table3.csv` and `.html` | PTEN estimates from five models for each secondary progression endpoint. |
| 25 | `JEI_figures_and_tables.zip` | Downloadable figures, tables, and model/audit results after an output audit. |

The relative figure and table paths above are under `jei_revision_outputs/` after a rerun. Saved manuscript references are under `jei_revision/reference_outputs/{figures,tables,results}`. These references contain figure, table, and result outputs; raw exports and patient-level datasets are not redistributed. Patient-level files generated by a rerun remain in `jei_revision_outputs/data_processed/`.

The original `figures/final_03_adjusted_cox_forest.png` shows only the PTEN estimate across the three primary models. It remains an original Phase 4 output; the revised manuscript Figure 3 is `Figure3_Cox_models.png`.

PFS and PFI use the original 498 eligible patients and the same 459 complete cases. PFS counts death from any cause; PFI censors deaths without tumor. Both endpoints include new primary tumors. The progression analyses are secondary, use nominal p values, and are not independent validation. The complete-case PTEN-only and age/stage progression models violated the proportional-hazards assumption for PTEN, so their HRs are summary estimates.

The JEI analysis script also runs the mutation-or-deletion sensitivity definition. It reclassifies three amplification-only cases as unaltered while retaining those patients. No additional standalone sensitivity-script step is required to produce the JEI outputs. Kaplan–Meier curves and medians use observed follow-up without extrapolation.

## Folder structure

| Path | Purpose |
|---|---|
| `PTEN_UCEC_OCSEF_master_notebook.ipynb` | Original Phases 1–4 plus JEI sections 20–25; setup and saved displays. |
| `ocsef_finalization.py` | Original Phase 4 reconstruction and figures. |
| `PTEN_amplification_excluded_sensitivity.py` | Standalone script from the original workflow. |
| `jei_revision/` | Integrated revision analysis, figure and table export scripts. |
| `jei_revision/requirements.txt` | Pinned dependency closure used for the complete notebook, including all original phases. |
| `jei_revision/reference_outputs/{figures,tables,results}` | Saved manuscript reference outputs without raw or patient-level datasets. |
| `data_raw/` | Original seven source exports supplied by the person running the notebook. |
| `reference_inputs/` | TCGA-CDR workbook downloaded and verified by setup. |
| `data_processed/`, `results/`, `figures/` | Locally generated original-phase outputs. |
| `jei_revision_outputs/` | Locally generated JEI revision outputs, including patient-level `data_processed/`. |
| `VALIDATION_REPORT.md`, `ANALYSIS_NOTES.md`, `DATA_DICTIONARY.md` | Execution evidence, analysis details, and input definitions. |
| `README.md`, `LICENSE`, `requirements.txt` | Project documentation, MIT license, and original-phase dependency list. |

## How to run

### Google Colab

1. Use the **Open in Colab** badge above to open the GitHub master notebook. The saved outputs can be viewed before rerunning it.
2. Start a fresh session and run the section 0 setup cell. It clones the repository when needed and installs `jei_revision/requirements.txt` for all phases. If Colab requests a restart after installation, restart the session and run all cells again.
3. Upload the original seven source exports together, or upload a ZIP containing them, when the setup cell prompts. The TCGA-CDR workbook is downloaded automatically and checked against its recorded checksum.
4. Select **Runtime → Run all**. Run all 27 code cells in order, including the original analyses and JEI sections 20–25. Section 25 audits the revision outputs and downloads `JEI_figures_and_tables.zip`.

### Local Python

1. Use Python 3.12 for the complete notebook and the dependencies in `jei_revision/requirements.txt`. The original Phase 1–4 execution used Python 3.13.5; see `VALIDATION_REPORT.md` for recorded execution and environment details.
2. Download or clone this repository and open a terminal in the project folder.
3. Create a folder named `data_raw/`.
4. Place the original source exports described above in `data_raw/` using these filenames:

   - `Supplementary_Table_S1_PTEN_altered_unaltered_sample_matrix.tsv`
   - `Supplementary_Table_S2_PTEN_discrete_CNA_table.tsv`
   - `Supplementary_Table_S3_PTEN_mutation_table.tsv`
   - `Supplementary_Table_S4_cBioPortal_clinical_survival_table.tsv`
   - `Supplementary_Table_S5_cBioPortal_PTEN_survival_summary_table.tsv`
   - `data_clinical_patient.txt`
   - `data_clinical_sample.txt`

5. Install the required packages:

   ```bash
   python -m pip install -r jei_revision/requirements.txt
   ```

6. Open `PTEN_UCEC_OCSEF_master_notebook.ipynb`.
7. Select **Restart Kernel and Run All Cells**. Setup downloads and verifies the TCGA-CDR workbook. The notebook runs Phases 1–3 in order, invokes `ocsef_finalization.py` for Phase 4, and then runs the JEI revision scripts and displays in sections 20–25. Later steps depend on the earlier result files.
8. Confirm that all 27 code cells finish without errors and review the section 25 output audit. The complete notebook contains 54 cells in total.
9. Review the original reconstruction and adjustment results in `results/`, and Figures 3–4, Tables 1–3, model diagnostics, checksums, and software versions in `jei_revision_outputs/`. Section 25 creates `JEI_figures_and_tables.zip`; patient-level datasets remain outside this archive.

## Interpretation of findings

> The supplied PTEN label was exactly reconstructed from mutation and high-level copy-number data. The independent survival rerun supported the portal's general unadjusted association but did not reproduce its exact cohort. PTEN alteration frequency differed substantially across molecular subtypes, and the survival association disappeared after adjustment for age, stage, and subtype. Therefore, this dataset does not support PTEN alteration as an independent survival predictor, although PTEN remains biologically important in UCEC.

## How to cite

Related publication:

> Liu N. A multi-omics secondary-data analysis of PTEN in uterine corpus
> endometrial carcinoma (UCEC) using public bioinformatics resources.
> *National High School Journal of Science*. Published August 2, 2026.
> https://nhsjs.com/2026/a-multi-omics-secondary-data-analysis-of-pten-in-uterine-corpus-endometrial-carcinoma-ucec-using-public-bioinformatics-resources/

The repository contains a distinct reproducibility and molecular-subtype
extension. Cite it by its title, author (Nathan Liu), year (2026), and repository
URL. If the AMIA 2026 High School Scholars abstract is accepted, its citation
will be added here.

The analysis code is available under the MIT License; see `LICENSE`.

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
