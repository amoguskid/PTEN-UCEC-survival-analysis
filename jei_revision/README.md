# JEI revision analyses

This extension regenerates manuscript-style Figures 1–4 and manuscript
Tables 1–3 from the original TCGA UCEC cohorts. Overall survival remains the
primary endpoint. PFS and PFI are reviewer-requested secondary endpoints.

The master notebook contains 56 cells, including 28 code cells, with the revision
analysis and display cells in sections 20–25. Section 20 generates and displays
Figures 1–2; sections 21–22 display Figures 3–4. Section 25 checks 19 required
output files and archives all four manuscript figures with the tables and results.
Its saved outputs show the revised figures and tables. Running all cells requires
the original seven source exports listed in the repository README. These exports
are not redistributed in the public repository. The setup cell downloads the
TCGA-CDR workbook from the NCI publication page and verifies its checksum.

## Run from a terminal

Install `jei_revision/requirements.txt` from the repository root:

```bash
python -m pip install -r jei_revision/requirements.txt
```

Place the original seven exports in `data_raw/`. Before running these commands,
run the notebook's section 0 setup cell to download and verify
`reference_inputs/TCGA-CDR-SupplementalTableS1.xlsx`, or download it from the
[official GDC source](https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81)
and verify the SHA-256 digest recorded in the repository README. Run the original
notebook phases first so that `results/subtype_statistics.csv` is available for
Figure 2. Then run from the repository root:

```bash
python jei_revision/analyze_survival_revision.py --raw data_raw --mutations data_raw/Supplementary_Table_S3_PTEN_mutation_table.tsv --cdr reference_inputs/TCGA-CDR-SupplementalTableS1.xlsx --out jei_revision_outputs
python jei_revision/describe_os_cohort.py --clinical data_raw/Supplementary_Table_S4_cBioPortal_clinical_survival_table.tsv --os jei_revision_outputs/data_processed/OS_eligible_cohort.csv --audit jei_revision_outputs/data_processed/all_498_patients_endpoint_and_covariate_audit.csv --out jei_revision_outputs/results
python jei_revision/make_figure3.py --coefficients jei_revision_outputs/results/all_cox_coefficients.csv --out jei_revision_outputs/figures
python jei_revision/supplemental_descriptive_audits.py --labels data_raw/Supplementary_Table_S1_PTEN_altered_unaltered_sample_matrix.tsv --mutations data_raw/Supplementary_Table_S3_PTEN_mutation_table.tsv --os jei_revision_outputs/data_processed/OS_eligible_cohort.csv --out jei_revision_outputs/results
python jei_revision/export_manuscript_tables.py --input jei_revision_outputs --out jei_revision_outputs
python jei_revision/make_figures1_2.py --os jei_revision_outputs/data_processed/OS_eligible_cohort.csv --subtype-statistics results/subtype_statistics.csv --out jei_revision_outputs/figures
```

Regenerated manuscript figures are saved under `jei_revision_outputs/figures/`:
`Figure1_Kaplan_Meier`, `Figure2_Subtype_frequency`, `Figure3_Cox_models`, and
`Figure4_PFS_PFI_Kaplan_Meier`, each as PNG and vector PDF. Saved references for
all four figures are under `jei_revision/reference_outputs/figures/`; table and
model references are in the adjacent `tables/` and `results/` folders. The original
Phase 4 `final_01`, `final_02`, and `final_03` figures remain in the repository-root
`figures/` output folder. Regenerated patient-level datasets are saved locally in
`jei_revision_outputs/data_processed`, separate from the references and archive.

## Scope and interpretation

- Both progression endpoints use the original 498 patients and 459 complete cases.
- TCGA-CDR PFS counts death from any cause. PFI censors deaths without tumor.
  Both endpoints include new primary tumors. Neither is restricted to confirmed
  UCEC-specific events.
- Intermediate progression HRs are summary estimates because PTEN violated
  proportional hazards in those models. Primary adjusted and stratified results
  remain nonsignificant. Secondary p values are nominal.
- The mutation-or-deletion sensitivity definition reclassifies three
  amplification-only cases while retaining those patients.
- Kaplan–Meier curves and medians use observed follow-up without extrapolation.
- No subtype-specific survival or null-only mutation analyses are added.

The recorded revision environment used Python 3.12.14 and the pinned dependencies. The
notebook records the actual environment used whenever it is rerun in Colab.
See `VALIDATION_REPORT.md` at the repository root for execution evidence.

Source: [TCGA Pan-Cancer Clinical Data Resource](https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018),
Liu et al. (2018), DOI 10.1016/j.cell.2018.02.052.
