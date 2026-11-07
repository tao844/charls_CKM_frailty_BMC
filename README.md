# Frailty and mortality across cardiovascular-kidney-metabolic syndrome stages: a joint-exposure analysis of CHARLS 2011–2020

Analysis code and aggregate result tables for the manuscript submitted to the **Yonsei Medical Journal** (Original Article).
**This repository contains no individual-level data.** CHARLS data are protected by a data-use agreement and must be obtained separately from the CHARLS team (http://charls.pku.edu.cn/).

## Cohort and primary numbers (rebuilt chain)

- Source: official CHARLS 2011 raw data; analytic cohort flow 17,705 → 13,256 (five-step exclusion chain in `results/Rebuild_flow_official.csv`)
- Main analysis set (mortality): 12,754 with derivable FI and CKM stage; complete-case model sample **n = 12,182, 1,438 deaths**
- Risk set (cognition): 9,140 free of baseline cognitive impairment; **1,908 events, 731 competing deaths**
- Exposures: CKM stage 0–2 vs 3–4 (AHA framework, published staging rules) × frailty status (41-item deficit-accumulation FI; robust <0.10, pre-frail 0.10–0.21, frail >0.21)
- Primary mortality model: age-as-time-scale Cox with delayed entry, 10 pre-specified covariates + age, household-clustered robust variance
- Cognition: education-specific cutoffs (each education stratum wave-1 mean − 1 SD), Cox with death as censoring + Fine–Gray subdistribution hazards

## Repository structure

| Path | Content |
|---|---|
| `R/00_setup.R` – `R/03_alt_fi_nockm.R` | Cohort derivation, death-module rebuild, exposure construction (stage-2 pipeline) |
| `R/04_mortality_final.R`, `R/05_cognition_final.R` | Primary outcome models (stage-2 pipeline) |
| `R/06_mi_m20.R`, `R/06b_mi_pooled_compare.R` | Multiple imputation (m = 20, seed 20261104) and pooled comparison |
| `R/07_cognition_lmm.R`, `R/07b_cog_lmm_sens.R`, `R/07c_cog_lmm_mechanism.R` | Cognitive trajectory mixed models |
| `R/08_discrimination.R` | Incremental discrimination (ΔC) |
| `R/09_figures.R` | Figure 1 (cohort flow), Figure 2 (KM curves), Figure 3 (Fine–Gray forest plot) |
| `R/_ms_spec.R` | Final primary models used for the manuscript tables (mortality six-cell table, within-stratum contrasts, trends, RERI; cognition Cox and Fine–Gray panels) |
| `R/10_timefixed_supplement.R` | Time-fixed (T1) specification of the mortality model, full cohort and excluding the 2011–2013 interval deaths; reproduces Supplementary Table S4 |
| `results/` | Aggregate result tables only (29 CSV + `_ms_spec.json`): main HR/SHR tables, rates, stepwise adjustment, RERI, trend tests, PH tests, time-varying HR appendix, MI pooled estimates, LMM slopes, flow and validation tables. `Supplementary_YMJ.docx` (manuscript supplement, Tables S1–S5) is generated from these aggregate tables |
| `sessionInfo.txt` | R session information for reproducibility |

## Reproducibility notes

1. Obtain CHARLS data (all waves incl. exit interviews / mortality follow-up) through the CHARLS application process; individual-level inputs are **not** part of this repository.
2. R ≥ 4.4 with packages: `survival` (required), `mice`, `cmprsk`, `nlme`, `haven`, `dplyr`, `ggplot2`, `jsonlite`.
3. Run pipeline in order: `00_setup.R` → `01_alt_cohort_flow.R` → `02_final_death_module.R` → `03_alt_exposures.R` → `03_alt_fi_nockm.R` → `_ms_spec.R`; secondary analyses `04`–`08` as applicable; figures via `09_figures.R`.
4. Frailty item list and CKM staging rules are fixed and documented in the manuscript supplement (`results/S_CKM_staging_map.csv`, `results/CKM_noblood_staging_rule.csv`).

## Citation / correspondence

Wanqing Tao, Department of Geriatrics, Southern Central Hospital of Yunnan Province (The First People's Hospital of Honghe State), Gejiu 661000, Yunnan, China. Email: 779562081@qq.com; ORCID: 0009-0000-8464-8565.
