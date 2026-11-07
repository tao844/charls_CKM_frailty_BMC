## ============================================================
## 09c_death_counts_two_cohorts.R —— 六组死亡数双口径对照表（修正表）
## 背景：审核打回"六组死亡数口径混用"——主文（全队列 10,556/871）与
##       风险集（8,916/684）两套死亡数必须分列、禁止交叉引用。
## 输出：output/tables/Death_events_both_cohorts.csv（唯一权威对照表）
## ============================================================
## 注意（2026-11-05 复核结论）：本表两口径必须用同一分组变量 joint6
## （= CKM_Stage_bl + frailty_bl 派生，与组分完全自洽，已交叉验证）。
## FI_noCKM 敏感性重分组 joint6_nc（frail3_nc 派生，与 frailty_bl 一致率 87.3%，
## 组大小不同，如 3~4 衰前为 885）属于 Major 7b 专用，禁止混入本表。
source("00_setup.R")
deaths <- read.csv("output/deaths.csv", fileEncoding = "UTF-8")

d_all <- merge(bl, deaths, by = "ID", all.x = TRUE)
stopifnot(!any(is.na(d_all$died)))
rs2  <- merge(rs, deaths, by = "ID", all.x = TRUE)
stopifnot(!any(is.na(rs2$died)))

agg <- function(dd) {
  ev  <- aggregate(died ~ joint6, data = dd, FUN = sum)
  py  <- aggregate(dtime_yr ~ joint6, data = dd, FUN = sum)
  n   <- aggregate(died ~ joint6, data = dd, FUN = length)
  m <- merge(merge(n, ev, by = "joint6"), py, by = "joint6")
  names(m) <- c("group", "n", "deaths", "person_years")
  m$rate100py <- round(m$deaths / m$person_years * 100, 2)
  m$person_years <- round(m$person_years, 1)
  m[match(jlev, m$group), ]
}

f <- agg(d_all)
r <- agg(rs2)
stopifnot(sum(f$deaths) == 871, sum(r$deaths) == 684,
          sum(f$n) == 10556, sum(r$n) == 8916)

out <- data.frame(
  group = f$group,
  n_full = f$n, deaths_full = f$deaths, py_full = f$person_years, rate_full = f$rate100py,
  n_riskset = r$n, deaths_riskset = r$deaths, py_riskset = r$person_years, rate_riskset = r$rate100py)
out$total_row <- ""
out <- rbind(out, data.frame(group = "合计", n_full = sum(f$n), deaths_full = sum(f$deaths),
  py_full = round(sum(f$person_years), 1), rate_full = round(sum(f$deaths) / sum(f$person_years) * 100, 2),
  n_riskset = sum(r$n), deaths_riskset = sum(r$deaths),
  py_riskset = round(sum(r$person_years), 1),
  rate_riskset = round(sum(r$deaths) / sum(r$person_years) * 100, 2), total_row = ""))
write.csv(out, "output/tables/Death_events_both_cohorts.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
print(out)
cat("\n校验：全队列", sum(f$deaths), "例死亡 / 风险集", sum(r$deaths), "例死亡 —— 与 871/684 关卡一致\n")
cat("09c 完成\n")
