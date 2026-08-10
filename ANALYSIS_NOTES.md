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

## What is new for OCSEF

The published manuscript relied mainly on portal-generated descriptive outputs. The current continuation work adds direct PTEN-label reconstruction, patient-level cohort cleaning and exclusion audits, an independent survival rerun, molecular-subtype association testing, multivariable and stratified Cox models, model diagnostics, reproducible code, and poster-ready figures. These new analyses—not the previously published portal outputs—should dominate the abstract, display, and interview.
