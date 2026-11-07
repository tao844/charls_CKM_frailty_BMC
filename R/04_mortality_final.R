## ============================================================
## 04_mortality_final.R —— 死亡主分析（Major 2/3/9/10/12）
## 【口径 2026-11-04 PI 定死】主分析集 = 全队列 10,556（871 例死亡）——返修核心改动；
##   风险集 8,916（684 例）作为敏感性分析重跑主模型。
## 主分析 = 标准全因死亡 Cox 模型，无竞争事件（Major 3a）
## 统一 13 协变量 + 家庭聚类稳健方差
## 输出顺序：主模型 HR 表 → 事件数/人年 → 阶梯调整 → 组内对比 → 趋势 → RERI → cox.zph
##          → 敏感性（9a-9f 族 / 风险集 / 死亡时间前移 12 个月）
## ============================================================
source("00_setup.R")
deaths <- read.csv("output/deaths.csv", fileEncoding = "UTF-8")

## ---- 全队列分析帧：对 bl 施加与 00 中 rs 相同的协变量预处理 ----
bl$CKM34  <- as.integer(bl$CKM_Stage_bl >= 3)
bl$frail2 <- as.integer(bl$frailty_bl == "衰弱")
bl$frailord <- as.integer(bl$frailty_bl == "衰弱前期") +
  2 * as.integer(bl$frailty_bl == "衰弱")
bl$ln_hhcperc <- log(bl$hhcperc_bl)
bl$exercise_bl <- factor(ifelse(is.na(bl$exercise_bl), "missing",
  ifelse(bl$exercise_bl == 1, "yes", "no")), levels = c("no", "yes", "missing"))
bl$pension_bl <- factor(ifelse(is.na(bl$pension_bl), "missing",
  ifelse(bl$pension_bl == 1, "yes", "no")), levels = c("no", "yes", "missing"))
bl$nation_bl <- factor(ifelse(is.na(bl$nation_bl), "missing",
  ifelse(bl$nation_bl == 1, "minority", "han")), levels = c("han", "minority", "missing"))

d_all <- merge(bl, deaths, by = "ID", all.x = TRUE)
stopifnot(!any(is.na(d_all$died)))
cat("校验：全队列死亡数 =", sum(d_all$died), "（应为 871）\n")
if (sum(d_all$died) != 871) stop("全队列死亡数 != 871，停下报告")

## 风险集帧（00 已预处理 rs；补死亡信息）
rs <- merge(rs, deaths, by = "ID", all.x = TRUE)
stopifnot(!any(is.na(rs$died)))
if (sum(rs$died) != 684) stop("风险集死亡数 != 684，停下报告")

## ---- 1) 主模型：全队列，六联合暴露组，统一 13 协变量，家庭聚类稳健方差 ----
f_main <- as.formula(paste("Surv(dtime_yr, died) ~ joint6 +",
                           paste(cov13, collapse = " + ")))
fit_main <- coxph(f_main, data = d_all, cluster = householdID,
                  x = TRUE, model = TRUE)
s_main <- summary(fit_main)
print(s_main$coefficients[grep("joint6", rownames(s_main$coefficients)), ])
write.csv(data.frame(term = rownames(s_main$coefficients), s_main$coefficients,
                     s_main$conf.int[, c("lower .95", "upper .95"), drop = FALSE]),
          "output/tables/Death_main_6groups.csv", row.names = FALSE, fileEncoding = "UTF-8")

## 每格事件数与人年（Table 2/3 用）
ev <- aggregate(died ~ joint6, data = d_all, FUN = sum)
py <- aggregate(dtime_yr ~ joint6, data = d_all, FUN = sum)
rate <- merge(ev, py, by = "joint6")
rate$rate100py <- round(rate$died / rate$dtime_yr * 100, 2)
write.csv(rate, "output/tables/Death_events_rates.csv", row.names = FALSE, fileEncoding = "UTF-8")
print(rate)

## ---- 2) 阶梯调整（Major 2：unadjusted → +age/sex → 全模型）----
fit0 <- coxph(Surv(dtime_yr, died) ~ joint6, data = d_all, cluster = householdID)
fit1 <- coxph(Surv(dtime_yr, died) ~ joint6 + age_bl + gender, data = d_all,
              cluster = householdID)
jn <- paste0("joint6", jlev[-1])
step_coef <- data.frame(
  group = jlev[-1],
  HR_unadj  = round(exp(coef(fit0)[jn]), 3),
  HR_agesex = round(exp(coef(fit1)[jn]), 3),
  HR_full   = round(exp(coef(fit_main)[jn]), 3))
print(step_coef)
write.csv(step_coef, "output/tables/Death_stepwise.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
blocks <- list(demo = c("age_bl", "gender", "edu_bl", "rural_bl", "marry_bl", "nation_bl"),
               behav = c("smoken_bl", "drinkev_bl", "exercise_bl", "sleep_bl"),
               econ = c("ln_hhcperc", "ins_bl", "pension_bl"))
dr <- sapply(names(blocks), function(nm) {
  f <- as.formula(paste("Surv(dtime_yr, died) ~ joint6 +",
                        paste(unlist(blocks[1:which(names(blocks) == nm)]), collapse = "+")))
  round(exp(coef(coxph(f, data = d_all, cluster = householdID))[jn]), 3)
})
write.csv(as.data.frame(dr), "output/tables/Death_blockwise_drivers.csv",
          fileEncoding = "UTF-8")

## ---- 3) 组内直接对比（Major 2，含 95%CI）----
cn <- names(coef(fit_main))
mkL <- function(num, den) { L <- numeric(length(cn)); L[cn == num] <- 1; L[cn == den] <- -1; L }
contr <- rbind(
  contrast_cox(fit_main, mkL("joint6CKM 3~4期 衰弱", "joint6CKM 3~4期 无衰弱"), "CKM3-4: 衰弱 vs 健壮"),
  contrast_cox(fit_main, mkL("joint6CKM 3~4期 衰弱前期", "joint6CKM 3~4期 无衰弱"), "CKM3-4: 衰弱前期 vs 健壮"),
  contrast_cox(fit_main, mkL("joint6CKM 0~2期 衰弱", "joint6CKM 0~2期 无衰弱"), "CKM0-2: 衰弱 vs 健壮"))
print(contr)
write.csv(contr, "output/tables/Death_within_stratum_contrasts.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

## ---- 4) 趋势检验（Major 2：各 CKM 层内 frailty 等级趋势）----
tr <- do.call(rbind, lapply(c(0, 1), function(s) {
  sub <- d_all[d_all$CKM34 == s, ]
  ft <- coxph(as.formula(paste(
    "Surv(dtime_yr, died) ~ frailord +", paste(cov13, collapse = "+"))),
    data = sub, cluster = householdID)
  data.frame(CKM = ifelse(s == 1, "3-4期", "0-2期"),
             HR_per_level = round(exp(coef(ft)["frailord"]), 3),
             p_trend = signif(summary(ft)$coefficients["frailord", "Pr(>|z|)"], 3))
}))
print(tr)
write.csv(tr, "output/tables/Death_trend_tests.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

## ---- 5) 乘加两尺度交互（Major 10）----
fit_ix <- coxph(Surv(dtime_yr, died) ~ CKM34 * frail2 + age_bl + gender +
                  edu_bl + rural_bl + marry_bl + nation_bl + smoken_bl +
                  drinkev_bl + exercise_bl + sleep_bl + ln_hhcperc +
                  ins_bl + pension_bl,
                data = d_all, cluster = householdID)
si <- summary(fit_ix)$coefficients
print(si[c("CKM34", "frail2", "CKM34:frail2"), ])
reri <- reri_cox(fit_ix, "CKM34", "frail2", "CKM34:frail2")
print(reri)
write.csv(reri, "output/tables/Death_RERI.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(round(exp(cbind(est = coef(fit_ix), summary(fit_ix)$conf.int)), 3),
          "output/tables/Death_interaction_components.csv",
          fileEncoding = "UTF-8")

## ---- 6) PH 诊断（Major 12：联合暴露项专项）----
zph <- cox.zph(fit_main)
print(zph)
write.csv(round(zph$table, 4), "output/tables/Death_PH_tests.csv",
          fileEncoding = "UTF-8")

## ---- 7) 敏感性分析族 ----
cat("\n-- 7a +CES-D --\n")
fit_cesd <- coxph(update(f_main, . ~ . + cesd10_bl), data = d_all, cluster = householdID)
s_cesd <- summary(fit_cesd)
gj <- grep("joint6", names(coef(fit_cesd)))
print(round(cbind(HR = round(exp(coef(fit_cesd)[gj]), 3),
                  lo = round(s_cesd$conf.int[gj, "lower .95"], 3),
                  hi = round(s_cesd$conf.int[gj, "upper .95"], 3)), 3))
write.csv(data.frame(term = names(coef(fit_cesd))[gj],
                     HR = round(exp(coef(fit_cesd)[gj]), 3),
                     lo = round(s_cesd$conf.int[gj, "lower .95"], 3),
                     hi = round(s_cesd$conf.int[gj, "upper .95"], 3),
                     p = signif(s_cesd$coefficients[gj, "Pr(>|z|)"], 3)),
          "output/tables/Death_sens_CESD.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("\n-- 7b +BMI（Major 7c）--\n")
fit_bmi <- coxph(update(f_main, . ~ . + bmi_bl), data = d_all, cluster = householdID)
s_bmi <- summary(fit_bmi)
gj <- grep("joint6", names(coef(fit_bmi)))
print(round(cbind(HR = round(exp(coef(fit_bmi)[gj]), 3),
                  lo = round(s_bmi$conf.int[gj, "lower .95"], 3),
                  hi = round(s_bmi$conf.int[gj, "upper .95"], 3)), 3))
write.csv(data.frame(term = names(coef(fit_bmi))[gj],
                     HR = round(exp(coef(fit_bmi)[gj]), 3),
                     lo = round(s_bmi$conf.int[gj, "lower .95"], 3),
                     hi = round(s_bmi$conf.int[gj, "upper .95"], 3),
                     p = signif(s_bmi$coefficients[gj, "Pr(>|z|)"], 3)),
          "output/tables/Death_sens_BMI.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("\n-- 7c 年龄为时间尺度（Major 9）--\n")
f_age <- as.formula(paste("Surv(age_bl, age_bl + dtime_yr, died) ~ joint6 +",
                          paste(setdiff(cov13, "age_bl"), collapse = " + ")))
fit_age <- coxph(f_age, data = d_all, cluster = householdID)
print(round(exp(coef(fit_age))[grep("joint6", names(coef(fit_age)))], 3))
s_age <- summary(fit_age)
gj <- grep("joint6", names(coef(fit_age)))
write.csv(data.frame(term = names(coef(fit_age))[gj],
                     HR = round(exp(coef(fit_age)[gj]), 3),
                     lo = round(s_age$conf.int[gj, "lower .95"], 3),
                     hi = round(s_age$conf.int[gj, "upper .95"], 3),
                     p = signif(s_age$coefficients[gj, "Pr(>|z|)"], 3)),
          "output/tables/Death_sens_agetimescale.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("\n-- 7d 有血样子样本（Major 6）--\n")
sub_b <- d_all[!is.na(d_all$bl_crea_bl), ]
fit_b <- coxph(f_main, data = sub_b, cluster = householdID)
print(round(exp(coef(fit_b))[grep("joint6", names(coef(fit_b)))], 3))
s_b <- summary(fit_b)
gj <- grep("joint6", names(coef(fit_b)))
write.csv(data.frame(term = names(coef(fit_b))[gj],
                     HR = round(exp(coef(fit_b)[gj]), 3),
                     lo = round(s_b$conf.int[gj, "lower .95"], 3),
                     hi = round(s_b$conf.int[gj, "upper .95"], 3),
                     p = signif(s_b$coefficients[gj, "Pr(>|z|)"], 3)),
          "output/tables/Death_sens_bloodsample.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("\n-- 7e 五分类 CKM 分期（Major 7a）--\n")
d_all$CKM5 <- factor(d_all$CKM_Stage_bl, levels = 0:4)
fit_5stg <- coxph(update(f_main, . ~ CKM5 + . - joint6), data = d_all,
                  cluster = householdID)
print(round(exp(coef(fit_5stg)), 3))
s5 <- summary(fit_5stg)
g5 <- grep("CKM5", rownames(s5$coefficients))
write.csv(data.frame(term = rownames(s5$coefficients)[g5],
                     HR = round(exp(coef(fit_5stg)[g5]), 3),
                     lo = round(s5$conf.int[g5, "lower .95"], 3),
                     hi = round(s5$conf.int[g5, "upper .95"], 3),
                     p = signif(s5$coefficients[g5, "Pr(>|z|)"], 3)),
          "output/tables/Death_sens_CKM5stage.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("\n-- 7f FI_noCKM 重分组重跑（Major 7b）--\n")
finc <- read.csv("output/tables/FI_noCKM.csv", fileEncoding = "UTF-8")
d_all$joint6_nc <- factor(finc$joint6_nc[match(d_all$ID, finc$ID)],
                          levels = c("CKM 0~2期 无衰弱", "CKM 0~2期 衰弱前期", "CKM 0~2期 衰弱",
                                     "CKM 3~4期 无衰弱", "CKM 3~4期 衰弱前期", "CKM 3~4期 衰弱"))
fit_nc <- coxph(update(f_main, . ~ . - joint6 + joint6_nc), data = d_all,
                cluster = householdID, x = TRUE)
s_nc <- summary(fit_nc)
print(s_nc$coefficients[grep("joint6_nc", rownames(s_nc$coefficients)), ])
write.csv(data.frame(term = rownames(s_nc$coefficients), s_nc$coefficients,
                     s_nc$conf.int[, c("lower .95", "upper .95"), drop = FALSE]),
          "output/tables/Death_main_FInoCKM.csv", row.names = FALSE, fileEncoding = "UTF-8")

## ---- 8) 敏感性：风险集 8,916（684 例死亡）重跑主模型 ----
cat("\n-- 8) 风险集敏感性（主模型原样重跑，n=8,916 / 684 事件）--\n")
fit_rs <- coxph(f_main, data = rs, cluster = householdID, x = TRUE, model = TRUE)
s_rs <- summary(fit_rs)
print(s_rs$coefficients[grep("joint6", rownames(s_rs$coefficients)), ])
write.csv(data.frame(term = rownames(s_rs$coefficients), s_rs$coefficients,
                     s_rs$conf.int[, c("lower .95", "upper .95"), drop = FALSE]),
          "output/tables/Death_main_riskset.csv", row.names = FALSE, fileEncoding = "UTF-8")
ev2 <- aggregate(died ~ joint6, data = rs, FUN = sum)
py2 <- aggregate(dtime_yr ~ joint6, data = rs, FUN = sum)
rate2 <- merge(ev2, py2, by = "joint6"); rate2$rate100py <- round(rate2$died / rate2$dtime_yr * 100, 2)
write.csv(rate2, "output/tables/Death_events_rates_riskset.csv", row.names = FALSE, fileEncoding = "UTF-8")
print(rate2)

## ---- 9) 敏感性：死亡时间前移 12 个月（Major 5）----
cat("\n-- 9) 死亡时间前移 12 个月（dtime_yr_sens 重跑主模型）--\n")
fit_ts <- coxph(update(f_main, Surv(dtime_yr_sens, died) ~ .),
                data = d_all, cluster = householdID)
s_ts <- summary(fit_ts)
ts_out <- data.frame(
  term = grep("joint6", names(coef(fit_ts)), value = TRUE),
  HR = round(exp(coef(fit_ts)[grep("joint6", names(coef(fit_ts)))]), 3),
  lo = round(s_ts$conf.int[grep("joint6", names(coef(fit_ts))), "lower .95"], 3),
  hi = round(s_ts$conf.int[grep("joint6", names(coef(fit_ts))), "upper .95"], 3),
  p  = signif(s_ts$coefficients[grep("joint6", rownames(s_ts$coefficients)), "Pr(>|z|)"], 3))
print(ts_out)
write.csv(ts_out, "output/tables/Death_main_time_sens.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

## ---- 10) PH 附录（已定口径：tt(joint6) 报 3/6/9 年 HR(t)；edu_bl/nation_bl 以 strata 吸收）----
## 因子型 joint6 需先转数值哑变量才能进 tt()
jd <- model.matrix(~ joint6, data = d_all)[, -1]
colnames(jd) <- paste0("jg", seq_len(ncol(jd)))
d_all <- cbind(d_all, jd)
f_ph <- as.formula(paste(
  "Surv(dtime_yr, died) ~", paste(colnames(jd), collapse = " + "),
  "+ strata(edu_bl) + strata(nation_bl) +",
  paste(setdiff(cov13, c("edu_bl", "nation_bl")), collapse = " + "),
  "+", paste(paste0("tt(", colnames(jd), ")"), collapse = " + ")))
fit_tt <- coxph(f_ph, data = d_all, cluster = householdID,
                tt = function(x, t, ...) x * t)
co <- coef(fit_tt)
hr_t <- do.call(rbind, lapply(seq_along(jlev[-1]), function(i) {
  bg <- unname(co[paste0("jg", i)]); bt <- unname(co[paste0("tt(jg", i, ")")])
  data.frame(group = jlev[-1][i],
             HR_3y = round(exp(bg + bt * 3), 3),
             HR_6y = round(exp(bg + bt * 6), 3),
             HR_9y = round(exp(bg + bt * 9), 3))
}))
print(hr_t)
write.csv(hr_t, "output/tables/Death_PH_timevarying_HR.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cat("04 完成\n")
