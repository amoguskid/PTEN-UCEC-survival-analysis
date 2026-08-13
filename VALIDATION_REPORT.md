# Validation Report

## Validation dates

- Clean-run analysis validation: **July 15, 2026**
- Public repository completion update: **August 13, 2026**

The July date records the two clean runs described below. The August date records the later update that added the complete reproducibility files to the public repository; it is not a second analysis-validation date.

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
