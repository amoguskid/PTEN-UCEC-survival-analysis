# Validation Report

## Validation dates

- JEI revision notebook validation: **October 7, 2026**
- Clean-run analysis validation: **July 15, 2026**
- Public repository completion update: **August 13, 2026**

The July date records the two clean runs described below. The August date records the later update that added the complete reproducibility files to the public repository; it is not a second analysis-validation date.

## JEI revision validation — October 7, 2026

The revised master notebook contains 54 cells, including 27 code cells. All 27
code cells were executed from top to bottom in one fresh IPython session with
Python 3.12.14 and the analysis versions pinned in
`jei_revision/requirements.txt`. Execution and display capture used IPython
directly because the validation workspace did not permit Jupyter kernel
sockets. No code cells failed. This was a local validation; the Colab setup and
upload/download paths are provided for the same notebook.

The run recreated the original four phases and the JEI revision outputs in
sections 20–25. It saved the actual execution outputs in the notebook, including
the revised Figure 3, the PFS/PFI Kaplan–Meier Figure 4, and Tables 1–3. The two
figures were visually checked, including the separation between the Figure 3
axis label and footnote. All 15 required revision output files were present.

The executed outputs were compared with the revised manuscript and the
validated revision analysis package:

- All 165 formatted cells in Table 1A, Table 1B, Table 2, and Table 3 matched the manuscript exactly, including headers.
- All 30 model rows, 78 coefficient rows, and 78 proportional-hazards diagnostic rows matched the reference CSV values exactly.
- All cohort counts, six age summaries, six Kaplan–Meier summaries, six log-rank tests, and 24 Figure 4 risk counts matched exactly.
- The seven endpoint/cohort and PTEN reconstruction datasets matched the reference CSV rows exactly.

The primary cohorts remained 498 endpoint-eligible patients and 459 complete
cases. The new endpoints contained 140 PFS events and 117 PFI events in the
498-patient cohort, and 125 PFS events and 107 PFI events among complete cases.
The fully adjusted PTEN estimates were HR 1.18 (95% CI 0.76–1.84; p=0.469) for
PFS and HR 1.14 (95% CI 0.69–1.86; p=0.610) for PFI. PFS includes deaths from
any cause; PFI censors deaths without tumor.

The OS median audit confirmed an observed Kaplan–Meier crossing at 110.55
months in the PTEN-unaltered group, reported as 110.6 months. Ten patients
remained at risk immediately before that time. The PTEN-altered curve did not
reach 0.5. Neither the curves nor medians were extrapolated.

Aggregate reference outputs are stored in `jei_revision/reference_outputs/`.
Raw exports and patient-level processed datasets are not redistributed. The
notebook can be opened in Colab from GitHub to view its saved outputs; a fresh
rerun requires the original seven source exports, while the TCGA-CDR workbook
is downloaded from the official NCI endpoint and checked against its SHA-256.
No new subtype-specific survival or null-only mutation analyses were added.

## Original validation record

## Clean-run procedure

1. Preserved the raw files, project documents, master notebook, and finalization script.
2. Deleted generated files from `data_processed/`, `results/`, and `figures/`.
3. Executed all 20 code cells from top to bottom with the pinned package versions in `requirements.txt`.
4. Confirmed that every required processed dataset, statistical table, diagnostic, and final figure was recreated.
5. Repeated the clean-run procedure in a fresh Python process.
6. Compared SHA-256 checksums for all generated CSV outputs and the three final poster PNGs.

## Result

- Code cells executed per run: **20**
- Execution errors: **0**
- Required outputs missing: **0**
- Phase 1 count checks passed: **6/6**
- Phase 2 count checks passed: **9/9**
- Direct PTEN reconstruction checks passed: **3/3**
- CSV and final-poster PNG outputs compared: **55**
- Identical across two clean runs: **55/55**

The final comparison is saved as `results/clean_run_final_checksum_comparison.csv`.

## Phase 1 — independent survival rerun

- Analyzable patients: **498**
- PTEN-altered: **146**
- PTEN-unaltered: **352**
- Deaths: **85**
- Log-rank p: **0.03977181**
- Unadjusted Cox HR: **0.59336**
- 95% CI: **0.35879–0.98130**

The direction was consistent with the portal output, but the 498-patient rerun was not an exact replication of the portal's 543-patient comparison.

## Phase 2 — molecular-subtype association

- PanCancer Atlas selected primary samples: **529**
- Exact PTEN-matrix matches: **529**
- Nonmissing subtype labels: **507**
- Missing subtype labels: **22**
- Chi-square: **53.272180**
- Degrees of freedom: **3**
- p-value: **1.60404940047e-11**
- Cramér's V: **0.324150**
- Minimum expected cell count: **16.7199**

Subtype-specific PTEN alteration frequencies:

- POLE-ultramutated: **16/49 (32.7%)**
- MSI-hypermutated: **59/148 (39.9%)**
- Copy-number low: **76/147 (51.7%)**
- Copy-number high: **22/163 (13.5%)**

## Phase 3 — adjusted Cox analysis

The same complete-case cohort was used in all three primary models:

- Patients: **459**
- Deaths: **77**
- PTEN-altered: **139**
- PTEN-unaltered: **320**

PTEN estimates:

- PTEN only: **HR 0.595, 95% CI 0.354–1.001, p=0.0504**
- PTEN + age + stage: **HR 0.826, 95% CI 0.483–1.411, p=0.4840**
- PTEN + age + stage + subtype: **HR 1.071, 95% CI 0.589–1.948, p=0.8214**

Advanced stage and the MSI indicator showed evidence of non-proportional hazards in the fully adjusted model. A stage/subtype-stratified sensitivity model gave **HR 1.139, 95% CI 0.623–2.080, p=0.6732** for PTEN, preserving the null conclusion.

## Phase 4 — direct PTEN reconstruction

- Samples reconstructed: **549**
- Samples with at least one PTEN mutation: **161**
- Samples with a high-level PTEN CNA (-2 or +2): **26**
- Samples with both: **6**
- Reconstructed PTEN-altered samples: **181**
- Portal-labeled PTEN-altered samples: **181**
- Exact sample-level matches: **549/549**
- Mismatches: **0**

## Final scientific interpretation

The supplied PTEN label was exactly reproducible from mutation and high-level CNA data. The independent survival rerun supported the portal's general unadjusted association but used a smaller cohort. PTEN alteration frequency differed substantially across molecular subtypes, and the PTEN survival estimate moved toward the null after adjustment for age, stage, and subtype. The data therefore do not support PTEN alteration as an independent survival predictor in this cohort. This observational secondary-data analysis does not establish causation or treatment response.
