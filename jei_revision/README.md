# JEI revision analyses

This extension regenerates the revised Figure 3, new Figure 4, and manuscript
Tables 1–3 from the original TCGA UCEC cohorts. Overall survival remains the
primary endpoint. PFS and PFI are reviewer-requested secondary endpoints.

The master notebook contains the setup and all display cells in sections 20–25.
Its saved outputs show the revised figures and tables. Running all cells requires
the original seven source exports listed in the repository README. These exports
are not redistributed in the public repository. The setup cell downloads the
TCGA-CDR workbook from the NCI publication page and verifies its checksum.

## Run from a terminal

Install `requirements.txt`, then run from the repository root:

```bash
python jei_revision/analyze_survival_revision.py --raw data_raw --mutations data_raw/Supplementary_Table_S3_PTEN_mutation_table.tsv --cdr reference_inputs/TCGA-CDR-SupplementalTableS1.xlsx --out jei_revision_outputs
python jei_revision/describe_os_cohort.py --clinical data_raw/Supplementary_Table_S4_cBioPortal_clinical_survival_table.tsv --os jei_revision_outputs/data_processed/OS_eligible_cohort.csv --audit jei_revision_outputs/data_processed/all_498_patients_endpoint_and_covariate_audit.csv --out jei_revision_outputs/results
python jei_revision/make_figure3.py --coefficients jei_revision_outputs/results/all_cox_coefficients.csv --out jei_revision_outputs/figures
python jei_revision/supplemental_descriptive_audits.py --labels data_raw/Supplementary_Table_S1_PTEN_altered_unaltered_sample_matrix.tsv --mutations data_raw/Supplementary_Table_S3_PTEN_mutation_table.tsv --os jei_revision_outputs/data_processed/OS_eligible_cohort.csv --out jei_revision_outputs/results
python jei_revision/export_manuscript_tables.py --input jei_revision_outputs --out jei_revision_outputs
```

The `reference_outputs` folder contains the validated figure, table, and model
outputs for the revised manuscript. Regenerated patient-level datasets are saved
locally in `jei_revision_outputs/data_processed`, separate from those references.

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

The validated revision used Python 3.12.14 and the pinned dependencies. The
notebook records the actual environment used whenever it is rerun in Colab.

Source: [TCGA Pan-Cancer Clinical Data Resource](https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018),
Liu et al. (2018), DOI 10.1016/j.cell.2018.02.052.
