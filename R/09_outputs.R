## ============================================================
## 09_outputs.R —— 汇总清单（跑完 01–08 后运行，核对回传完整性）
## ============================================================
files <- list.files("output/tables", pattern = "\\.csv$", full.names = TRUE)
cat("共", length(files), "个结果文件：\n")
print(basename(files))

## ---- 回传前自查（缺一不可）----
checklist <- c(
  "Flow_steps_check.csv", "Gap_11982_vs_10556.csv",
  "Unexplained_exclusions_profile.csv",
  "Included_vs_excluded.csv",              # Major 6
  "FI_noCKM.csv",                          # Major 7b（个体级，仅本地留档，不进公开包）
  "Table1_riskset.csv",                    # Minor 风险集基线表
  "CKM_stage_distribution.csv", "S_CKM_staging_map.csv",   # Major 7a
  "Death_main_6groups.csv", "Death_events_rates.csv", "Death_stepwise.csv",
  "Death_within_stratum_contrasts.csv", "Death_trend_tests.csv",
  "Death_RERI.csv", "Death_PH_tests.csv",
  "CI_main_6groups.csv", "CI_events_rates.csv", "CI_stepwise.csv",
  "CI_within_stratum_contrasts.csv", "CI_trend_tests.csv",
  "CI_education_cutoff_values.csv", "CI_confirmed_2waves.csv",
  "CI_FineGray_secondary.csv",
  "CI_IPC_weighted.csv", "CI_RERI.csv", "CI_PH_tests.csv",
  "Death_MI20_pooled_full.csv", "CI_MI20_pooled.csv",    # Major 11
  "LMM_joint_time.csv", "LMM_group_slopes.csv",          # Major 8c
  "Discrimination_death.csv",                            # Major 10
  "MI20_rvi_fmi_summary.csv",
  "Death_CC_vs_MI20.csv", "CI_CC_vs_MI20.csv"
)
cat("\n== 回传自查清单 ==\n")
miss <- character(0)
for (f in checklist) {
  ok <- file.exists(file.path("output/tables", f))
  cat(ifelse(ok, "[v]", "[x]"), f, "\n")
  if (!ok) miss <- c(miss, f)
}
cat("\n缺失", length(miss), "个文件"); if (length(miss)) cat("：", paste(miss, collapse = ", "))
cat("\n")

## 其他必备产物
for (f in c("output/sessionInfo.txt", "output/deaths.csv")) {
  cat(ifelse(file.exists(f), "[v]", "[x]"), f, "\n")
}
cat("09 完成——请回传 output 文件夹\n")
