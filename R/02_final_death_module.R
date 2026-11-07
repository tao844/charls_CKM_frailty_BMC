## ============================================================
## 02_final_death_module.R —— 死亡结局构建（Major 5）
## 产出：output/deaths.csv（ID, died, dtime_yr）
##   died：随访期内死亡 =1
##   dtime_yr：自基线（2011 首次有效访视月）至"首次记录死亡/末次在访波次"的月数/12
## 数据源：CHARLS_exit_csv/.../analysis_db_mortality.xlsx
##   （自 Harmonized CHARLS Version D 推导：died + indv_l_died 定结局，
##     iyear/imonth 定时间；论文方法学仅写"自 Harmonized CHARLS Version D 推导"）
## 校验关卡：全队列 10,556 死亡 = 871；风险集 8,916 死亡 = 684
## ============================================================
source("00_setup.R")

mort_path <- file.path("CHARLS_exit_csv",
  "CHARLS — 死亡（Mortality） (dbcodebook.cn)", "analysis_db_mortality.xlsx")
mort <- readxl::read_excel(mort_path)
cat("死亡库行数：", nrow(mort), "；其中死亡：", sum(mort$Die_status == 1, na.rm = TRUE), "\n")

## ---- 匹配到分析样本 ----
deaths <- data.frame(
  ID       = bl$ID,
  died     = as.integer(mort$Die_status[match(bl$ID, mort$id)] == 1),
  dtime_yr = mort$Die_time_months[match(bl$ID, mort$id)] / 12,
  wave_yr  = mort$year[match(bl$ID, mort$id)]          # 首次有效访视波年份（全队列均为 2011）
)
## Major 5 敏感性：死亡时间前移 12 个月（下限 1/12 年，避免非正时间）
deaths$dtime_yr_sens <- ifelse(deaths$died == 1,
                               pmax(deaths$dtime_yr - 1, 1/12),
                               deaths$dtime_yr)
stopifnot(!any(is.na(deaths$died)), !any(is.na(deaths$dtime_yr)))

## ---- 校验（代码包三关 + 本项目两关）----
stopifnot(all(deaths$ID %in% bl$ID))
stopifnot(all(deaths$dtime_yr > 0 | deaths$died == 0))
n_dead_all <- sum(deaths$died)
n_dead_rs  <- sum(deaths$died[bl$cog_impair_bl == 0])
cat("死亡人数：", n_dead_all, "（全队列；预期 871）\n")
cat("风险集死亡人数：", n_dead_rs, "（预期 684）\n")
cat("基线认知障碍者中死亡：", n_dead_all - n_dead_rs,
    "（应恰为 871-684=187）\n")
if (n_dead_all != 871) stop("全队列死亡数 ", n_dead_all, " != 871，停下报告")
if (n_dead_rs  != 684) stop("风险集死亡数 ", n_dead_rs, " != 684，停下报告")

## ---- 死亡时间分布（年）----
cat("死亡时间分布（年）：\n")
print(round(quantile(deaths$dtime_yr[deaths$died == 1],
      c(.05, .25, .5, .75, .95)), 2))
yrcut <- cut(deaths$dtime_yr[deaths$died == 1],
             breaks = c(0, 2, 4, 5, 7, 9, Inf),
             labels = c("0-2", "2-4", "4-5", "5-7", "7-9", ">9"))
cat("各时间段死亡人数：\n"); print(table(yrcut))

## ---- 精度构成表（回复 Major 5：死亡时间如何确定）----
dd <- deaths[deaths$died == 1, ]
## 判定：死亡月是否落在各确认波访视窗口（24/48/84/108 月 ±3）内
wave_m <- c(24, 48, 84, 108)
dd$timing_src <- vapply(dd$dtime_yr * 12, function(mo) {
  hit <- which(abs(wave_m - mo) <= 3)
  if (length(hit)) paste0("确认波访视时点(", wave_m[hit], "月)") else "波次间隔内（区间）"
}, character(1))
prec <- as.data.frame(table(首次有效访视年份 = dd$wave_yr, 时间精度 = dd$timing_src))
write.csv(prec, "output/tables/Death_timing_precision.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("死亡时间精度构成：\n"); print(prec)

write.csv(deaths, "output/deaths.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("02 完成\n")
