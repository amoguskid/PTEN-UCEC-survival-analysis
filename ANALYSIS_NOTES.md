# Analysis Notes

## Cohort logic

### Phase 1

The survival rerun begins with the 549-row PTEN matrix, matches 500 records to the supplied clinical export, and excludes two additional patients with unusable survival time. The final survival cohort contains 498 patients.

### Phase 2

The molecular-subtype analysis deliberately does not begin with the survival-complete cohort.

The PanCancer Atlas sample file contains one selected primary sample for each of 529 patients. Every selected sample exactly matches a sample in the Firehose PTEN matrix. The patient file provides subtype labels for 507 of these 529 patients. The 22 missing labels are all marked as outside the PanCancer pathway-analysis freeze in the supplied patient file.

This sample-level merge also resolves the only duplicated Firehose patient: one patient has two Firehose samples with discordant PTEN status, while the PanCancer Atlas sample file selects the primary sample ending in `-01`.

### Phase 3

The adjusted Cox analysis starts from the 498-patient survival cohort. It uses one fixed complete-case cohort for all three primary models so that movement in the PTEN hazard ratio is caused by covariate adjustment rather than a changing sample set.

- 498 patients available for survival analysis
- 2 missing age values
- 37 missing or unrecognized molecular-subtype labels
- 459 complete cases
- 77 deaths
- 139 PTEN-altered and 320 PTEN-unaltered patients

Age is centered at the cohort median and scaled per 10 years. Stage I–II is coded as early stage and Stage III–IV as advanced stage. Copy-number high is the subtype reference category.

### Phase 4

PTEN status is independently reconstructed from the supplied files using the same rule represented by the portal label:

- any PTEN mutation in the mutation table; or
- high-level discrete CNA equal to -2 or +2.

Low-level gains (+1) and shallow losses (-1) are not classified as altered. The reconstruction matched all 549 supplied labels exactly.

### Integrated JEI revision: sections 20–25

The master notebook now contains 54 cells, including 27 code cells. Section 0 sets up Colab, installs the pinned dependency closure in `jei_revision/requirements.txt` for all phases, and checks the inputs. A fresh rerun requires the original seven exports listed in `README.md`, supplied as individual files or a ZIP. These raw exports are not redistributed.

The setup cell automatically downloads the TCGA-CDR workbook from the official NCI/GDC source linked on the [TCGA Pan-Cancer Clinical Data Resource publication page](https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018). It verifies SHA-256 `ea594c0fbb6731477c7ac511fab449ca9c38b0d42d269591ed9f5c4090e75a5a` before using the workbook. Revision outputs record input checksums and the actual runtime package versions. The local revision environment uses Python 3.12.14; execution evidence belongs in `VALIDATION_REPORT.md`.

Overall survival remains the primary endpoint. Reviewer-requested PFS and PFI analyses use the original 498-patient eligible cohort and the same 459 complete cases. Patient identifiers are matched exactly to TCGA-CDR. PFS counts death from any cause; PFI censors deaths without tumor. Both endpoints include new primary tumors, so neither is limited to confirmed UCEC-specific events. The secondary p values are nominal, and these analyses do not provide an independent validation cohort.

`jei_revision/analyze_survival_revision.py` reconstructs both the original PTEN definition and a mutation-or-deletion sensitivity definition. The latter reclassifies the three amplification-only cases as unaltered while retaining them in the cohort. The other amplification cases also have PTEN mutations and remain altered. The notebook runs this sensitivity definition within the revision analysis; it requires no separate invocation of the original standalone script.

| Notebook section | Saved display and generated artifact |
|---|---|
| 20 | Revision models, diagnostics, descriptive audits, and sensitivity analysis. |
| 21 | Revised manuscript Figure 3, including the complete-case PTEN-only estimate and all six fully adjusted OS coefficients. |
| 22 | Figure 4: PFS and PFI Kaplan–Meier curves with confidence bands, censoring marks, and numbers at risk. |
| 23 | Tables 1–2: reconstruction and survival-cohort characteristics including age; all six fully adjusted OS coefficients. |
| 24 | Table 3: PTEN estimates from five models for each progression endpoint. |
| 25 | Output audit and `JEI_figures_and_tables.zip` download. |

Reruns write revision artifacts to `jei_revision_outputs/`: figures as PNG and vector PDF, and manuscript tables as CSV and HTML. The original Phase 4 `final_03_adjusted_cox_forest` plots PTEN estimates across three models; the revised manuscript figure is `Figure3_Cox_models`, which also shows every fully adjusted coefficient. Figure 4 is `Figure4_PFS_PFI_Kaplan_Meier`.

Saved reference outputs are under `jei_revision/reference_outputs/{figures,tables,results}` without raw exports or patient-level datasets. Patient-level files generated by a rerun remain under `jei_revision_outputs/data_processed/`. The download archive includes only figures, tables, and model/audit result files.

Table 1 reports age availability and missingness, age summaries, and Kaplan–Meier medians from observed data. The unaltered median is 110.6 months using half-up rounding of 110.55; the altered curve does not reach 50% survival, so its median is not reached. Figure 4 displays a 120-month window, while the models use all recorded follow-up. Curves and medians are not extrapolated.

## Statistical plan

- Primary table: four molecular subtypes × two PTEN-status groups
- Primary test: Pearson chi-square
- Assumption check: all expected counts must be at least five
- Effect size: Cramér's V
- Subtype estimates: PTEN-altered percentage with Wilson 95% confidence interval
- Cell interpretation: adjusted standardized residuals
- Exploratory post-hoc tests: pairwise Fisher exact tests with Holm correction

## Prespecified subtype mapping

- `UCEC_POLE` → POLE-ultramutated
- `UCEC_MSI` → MSI-hypermutated
- `UCEC_CN_LOW` → Copy-number low
- `UCEC_CN_HIGH` → Copy-number high

The raw labels remain in all processed files.

## Cross-release limitation

The PTEN status matrix comes from the Firehose Legacy analysis, while the subtype labels come from the PanCancer Atlas clinical data. Exact TCGA identifiers permit deterministic matching, but the source releases are not identical. This must be stated on the poster and during judging.

## Adjusted-model diagnostics

The proportional-hazards test flagged advanced stage and the MSI subtype indicator in the fully adjusted model. A sensitivity model stratified by stage and subtype allowed different baseline hazards across those groups. Its PTEN estimate remained null (HR 1.139, 95% CI 0.623–2.080; p=0.6732), supporting the primary conclusion without claiming that every model assumption was perfect.

For the secondary progression endpoints, PTEN violated proportional hazards in the complete-case PTEN-only and age/stage models. Their HRs are summary estimates. The fully adjusted and stratified progression models did not detect a significant PTEN association. Table 3 reports the diagnostics and qualifications alongside those estimates.

## Scope of the extension

The published manuscript relied mainly on portal-generated descriptive outputs. The current continuation work adds direct PTEN-label reconstruction, patient-level cohort cleaning and exclusion audits, an independent survival rerun, molecular-subtype association testing, multivariable and stratified Cox models, model diagnostics, reproducible code, and final figures. These new analyses are distinct from the previously published portal outputs and define the scope of this extension.

The integrated JEI revision adds the requested progression endpoints and revised manuscript figures and tables within that scope. No subtype-specific survival or null-only mutation analyses are added.
