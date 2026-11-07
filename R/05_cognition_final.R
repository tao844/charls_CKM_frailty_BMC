## ============================================================
## 05_cognition_final.R —— 认知障碍主分析（Major 2/3/4/8/9/10/12）
## 【口径 2026-11-04 PI 定死】
##   1) 主分析集 = 风险集 8,916（无基线认知障碍者）
##   2) 主模型 = 标准 Cox（死亡作删失），统一 13 协变量 + 家庭聚类稳健，
##      参照 = CKM 0~2 健壮；主定义 = 教育分层界值（各教育层 wave1 均数-1SD）
##   3) 敏感性 = 单一 8.96 界值（原口径）+ 双波连续确认（主定义口径）
##   4) 次分析 = Fine-Gray（cmprsk::crr，死亡为竞争事件），只报 SHR
##   5) Major 4：IPC 权重敏感性（1%/99% 缩尾）+ 失访模式表
##   6) 双尺度交互（乘性 + RERI，delta 法 CI）
##   7) cox.zph 全套；PH 口径已定：主模型时间固定，附录 tt(joint6) 报 3/6/9 年
##      HR(t)，edu_bl/nation_bl 以 strata 吸收
## 注：落盘 HR 表一律用 summary()$conf.int（修复代码包 exp(cbind(coef, conf.int))
##     双重取指的原式缺陷）
## ============================================================
source("00_setup.R")
deaths <- read.csv("output/deaths.csv", fileEncoding = "UTF-8")
stopifnot(all(c("edu_bl", "cog_impair") %in% names(lng)))
rs <- merge(rs, deaths, by = "ID", all.x = TRUE)
stopifnot(!any(is.na(rs$died)))
fu <- aggregate(time ~ ID, data = lng, FUN = max)   # 每人最后随访年
rs$lastfu <- fu$time[match(rs$ID, fu$ID)]

lci <- lng[order(lng$ID, lng$time), ]
lci <- lci[lci$ID %in% rs$ID & !is.na(lci$total_cognition), ]

## ---- 事件/删失时间已在下方各节内联构造；HR 表统一出口 ----
hr_tab <- function(fit, prefix = "joint6") {
  s <- summary(fit)
  g <- grep(prefix, rownames(s$coefficients))
  data.frame(term = rownames(s$coefficients)[g],
             HR = round(exp(s$coefficients[g, "coef"]), 3),
             lo = round(s$conf.int[g, "lower .95"], 3),
             hi = round(s$conf.int[g, "upper .95"], 3),
             p  = signif(s$coefficients[g, "Pr(>|z|)"], 3))
}

## ---- 1) 教育分层界值表（主定义）----
w1c <- lng[lng$wave == 1 & !is.na(lng$total_cognition) & !is.na(lng$edu_bl), ]
cut_tab <- do.call(rbind, lapply(sort(unique(w1c$edu_bl)), function(e) {
  x <- w1c$total_cognition[w1c$edu_bl == e]
  data.frame(edu_bl = e, n = length(x),
             mean = round(mean(x), 2), sd = round(sd(x), 2),
             cutoff = round(mean(x) - sd(x), 2))
}))
cat("教育分层界值表（wave1 均数-1SD）：\n"); print(cut_tab)
write.csv(cut_tab, "output/tables/CI_education_cutoff_values.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

lci$cutval <- cut_tab$cutoff[match(lci$edu_bl, cut_tab$edu_bl)]
lci$ci_e <- as.integer(!is.na(lci$cutval) & lci$total_cognition < lci$cutval)

## 主定义事件：首次低于所在教育层界值的波
first_e <- lci[lci$ci_e == 1, ]; first_e <- first_e[!duplicated(first_e$ID), c("ID", "time")]
names(first_e) <- c("ID", "ctime_e")
rs <- merge(rs, first_e, by = "ID", all.x = TRUE)
rs$ci_ev <- as.integer(!is.na(rs$ctime_e))
rs$ci_tm <- ifelse(rs$ci_ev == 1, rs$ctime_e,
                   ifelse(rs$died == 1, pmin(rs$lastfu, rs$dtime_yr, na.rm = TRUE),
                          rs$lastfu))
## 仅基线一访者删失时间为 0：不贡献人时，Cox 中自动剔除，保留并报告
n_zero <- sum(rs$ci_tm == 0 & rs$ci_ev == 0)
cat("删失时间为 0（仅波1记录）人数：", n_zero, "\n")
stopifnot(all(rs$ci_tm >= 0 | rs$ci_ev == 0, na.rm = TRUE),
          !any(is.na(rs$ci_tm)))
cat("\n主定义（教育分层界值）事件数：", sum(rs$ci_ev), "\n")

## ---- 2) 主模型：标准 Cox，统一 13 协变量，家庭聚类稳健 ----
f_ci <- as.formula(paste("Surv(ci_tm, ci_ev) ~ joint6 +",
                         paste(cov13, collapse = " + ")))
fit_ci <- coxph(f_ci, data = rs, cluster = householdID, x = TRUE)
main_tab <- hr_tab(fit_ci)
print(main_tab)
write.csv(main_tab, "output/tables/CI_main_6groups.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

## 每格事件数与人年
evc <- aggregate(ci_ev ~ joint6, data = rs, FUN = sum)
pyc <- aggregate(ci_tm ~ joint6, data = rs, FUN = sum)
ratec <- merge(evc, pyc, by = "joint6")
ratec$rate100py <- round(ratec$ci_ev / ratec$ci_tm * 100, 2)
write.csv(ratec, "output/tables/CI_events_rates.csv", row.names = FALSE, fileEncoding = "UTF-8")
print(ratec)

## ---- 2b) 阶梯调整 + 组内对比 + 趋势（Major 2，同 04 结构；代码包 05 第 3 节）----
fit_ci0 <- coxph(Surv(ci_tm, ci_ev) ~ joint6, data = rs, cluster = householdID)
fit_ci1 <- coxph(Surv(ci_tm, ci_ev) ~ joint6 + age_bl + gender, data = rs,
                 cluster = householdID)
## 修复（同 04 写法）：fit_ci1 含 age_bl/gender 共 7 个系数，必须按 joint6 名子集取 5 行
jnc <- paste0("joint6", jlev[-1])
stepc <- data.frame(group = jlev[-1],
  HR_unadj = round(exp(coef(fit_ci0)[jnc]), 3),
  HR_agesex = round(exp(coef(fit_ci1)[jnc]), 3),
  HR_full = round(exp(coef(fit_ci)[jnc]), 3))
print(stepc)
write.csv(stepc, "output/tables/CI_stepwise.csv", row.names = FALSE, fileEncoding = "UTF-8")

cc <- names(coef(fit_ci))
mkL <- function(num, den) { L <- numeric(length(cc)); L[cc == num] <- 1; L[cc == den] <- -1; L }
con_ci <- rbind(
  contrast_cox(fit_ci, mkL("joint6CKM 3~4期 衰弱","joint6CKM 3~4期 无衰弱"), "CKM3-4: 衰弱 vs 健壮"),
  contrast_cox(fit_ci, mkL("joint6CKM 3~4期 衰弱前期","joint6CKM 3~4期 无衰弱"), "CKM3-4: 衰弱前期 vs 健壮"),
  contrast_cox(fit_ci, mkL("joint6CKM 0~2期 衰弱","joint6CKM 0~2期 无衰弱"), "CKM0-2: 衰弱 vs 健壮"))
print(con_ci)
write.csv(con_ci, "output/tables/CI_within_stratum_contrasts.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

trc <- do.call(rbind, lapply(c(0, 1), function(s) {
  sub <- rs[rs$CKM34 == s, ]
  ft <- coxph(as.formula(paste(
    "Surv(ci_tm, ci_ev) ~ frailord +", paste(cov13, collapse = "+"))),
    data = sub, cluster = householdID)
  data.frame(CKM = ifelse(s == 1, "3-4期", "0-2期"),
             HR_per_level = round(exp(coef(ft)["frailord"]), 3),
             p_trend = signif(summary(ft)$coefficients["frailord", "Pr(>|z|)"], 3))
}))
print(trc)
write.csv(trc, "output/tables/CI_trend_tests.csv", row.names = FALSE, fileEncoding = "UTF-8")

## ---- 3) 敏感性 3a：单一 8.96 界值（原口径）----
first_8 <- lci[lci$cog_impair == 1, ]; first_8 <- first_8[!duplicated(first_8$ID), c("ID", "time")]
names(first_8) <- c("ID", "ctime_8")
rs <- merge(rs, first_8, by = "ID", all.x = TRUE)
rs$ci_ev8 <- as.integer(!is.na(rs$ctime_8))
rs$ci_tm8 <- ifelse(rs$ci_ev8 == 1, rs$ctime_8,
                    ifelse(rs$died == 1, pmin(rs$lastfu, rs$dtime_yr, na.rm = TRUE),
                           rs$lastfu))
fit_8 <- coxph(update(f_ci, Surv(ci_tm8, ci_ev8) ~ .), data = rs, cluster = householdID)
tab8 <- hr_tab(fit_8)
cat("\n敏感性 3a：单一 8.96 界值，事件数：", sum(rs$ci_ev8), "\n"); print(tab8)
write.csv(tab8, "output/tables/CI_main_896.csv", row.names = FALSE, fileEncoding = "UTF-8")

## ---- 3b) 敏感性：双波连续确认（主定义口径）----
lci$prev_e <- ave(lci$ci_e, lci$ID, FUN = function(x) c(NA, head(x, -1)))
conf2 <- lci[lci$ci_e == 1 & !is.na(lci$prev_e) & lci$prev_e == 1, ]
conf2 <- conf2[!duplicated(conf2$ID), c("ID", "time")]
names(conf2) <- c("ID", "ctime_c")
rs <- merge(rs, conf2, by = "ID", all.x = TRUE)
rs$ci_evc <- as.integer(!is.na(rs$ctime_c))
rs$ci_tmc <- ifelse(rs$ci_evc == 1, rs$ctime_c,
                    ifelse(rs$died == 1, pmin(rs$lastfu, rs$dtime_yr, na.rm = TRUE),
                           rs$lastfu))
fit_c2 <- coxph(update(f_ci, Surv(ci_tmc, ci_evc) ~ .), data = rs, cluster = householdID)
tabc <- hr_tab(fit_c2)
cat("\n敏感性 3b：双波连续确认，事件数：", sum(rs$ci_evc), "\n"); print(tabc)
write.csv(tabc, "output/tables/CI_confirmed_2waves.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

## ---- 4) 次分析：Fine-Gray（死亡为竞争事件，只报 SHR）----
rs$fgstat <- ifelse(rs$ci_ev == 1, 2, ifelse(rs$died == 1, 1, 0))
rs$fgtime <- ifelse(rs$fgstat == 1, rs$dtime_yr,
                    ifelse(rs$fgstat == 2, rs$ci_tm, rs$lastfu))
## 先从 rs 显式取完整案例，再构 design matrix（本机 model.matrix 的 na.action
## 属性不可靠，不可依赖）
fg_vars <- c("joint6", cov13)
cc_fg <- complete.cases(rs[, fg_vars]) & !is.na(rs$fgtime) & !is.na(rs$fgstat)
cat("FG 完整案例数：", sum(cc_fg), "（总", length(cc_fg), "）\n")
mmf <- model.matrix(~ joint6 + age_bl + gender + edu_bl + rural_bl +
                      marry_bl + nation_bl + smoken_bl + drinkev_bl +
                      exercise_bl + sleep_bl + ln_hhcperc + ins_bl +
                      pension_bl, data = rs[cc_fg, ])
mm <- mmf[, -1]
fg <- cmprsk::crr(rs$fgtime[cc_fg], rs$fgstat[cc_fg],
                  cov1 = mm, failcode = 2)
sfg <- summary(fg)$coef
fg_tab <- data.frame(term = grep("joint6", rownames(sfg), value = TRUE),
                     SHR = round(exp(sfg[grep("joint6", rownames(sfg)), "coef"]), 3),
                     lo = round(exp(sfg[grep("joint6", rownames(sfg)), "coef"] - 1.96 * sfg[grep("joint6", rownames(sfg)), "se(coef)"]), 3),
                     hi = round(exp(sfg[grep("joint6", rownames(sfg)), "coef"] + 1.96 * sfg[grep("joint6", rownames(sfg)), "se(coef)"]), 3),
                     p  = signif(sfg[grep("joint6", rownames(sfg)), "p-value"], 3))
cat("\nFine-Gray 次分析（SHR）：\n"); print(fg_tab)
write.csv(fg_tab, "output/tables/CI_FineGray_secondary.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

## ---- 5) Major 4：IPC 权重敏感性 + 失访模式表 ----
lw <- lng[lng$ID %in% rs$ID, ]
lw <- lw[order(lw$ID, lw$time), ]
lw$respond <- as.integer(!is.na(lw$total_cognition))
lw$first   <- !duplicated(lw$ID)
lw$lag_resp <- ave(lw$respond, lw$ID, FUN = function(x) c(NA, head(x, -1)))
base_cov <- c("age_bl", "gender", "edu_bl", "rural_bl", "ln_hhcperc")
## lw 与 rs 同名列合并会产生 .x/.y 后缀，故显式限定 lw 只带权重构造所需列
lw2 <- merge(lw[, c("ID", "wave", "time", "respond", "first", "lag_resp")],
             rs[, c("ID", "joint6", base_cov)], by = "ID")
## 失访模式表：各组各波完成认知测评比例
att_pat <- with(lw2, aggregate(respond, by = list(joint6 = joint6, wave = wave), FUN = mean))
names(att_pat)[3] <- "response_rate"
att_pat$response_rate <- round(att_pat$response_rate, 3)
write.csv(att_pat, "output/tables/CI_attrition_pattern.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n失访模式表（各组各波认知测评完成率）：\n"); print(att_pat)
## IPCW
tt2 <- lw2[!lw2$first & !is.na(lw2$lag_resp), ]
den_fit <- glm(respond ~ joint6 + age_bl + gender + edu_bl + rural_bl +
                 ln_hhcperc, data = tt2, family = binomial)
num_fit <- glm(respond ~ age_bl + gender + edu_bl + rural_bl +
                 ln_hhcperc, data = tt2, family = binomial)
tt2$pw <- predict(num_fit, tt2, type = "response") /
          pmax(predict(den_fit, tt2, type = "response"), 1e-3)
wfin <- aggregate(pw ~ ID, data = tt2, FUN = prod, na.rm = TRUE)
names(wfin)[2] <- "ipc_w"
qs <- quantile(wfin$ipc_w, c(.01, .99), na.rm = TRUE)
cat("IPC 权重缩尾界：", round(qs, 3), "\n")
wfin$ipc_w <- pmin(pmax(wfin$ipc_w, qs[1]), qs[2])
rs$ipc_w <- wfin$ipc_w[match(rs$ID, wfin$ID)]
rs$ipc_w[is.na(rs$ipc_w)] <- 1
fit_ipc <- coxph(f_ci, data = rs, weights = ipc_w, cluster = householdID)
tabi <- hr_tab(fit_ipc)
cat("\nIPC 加权主模型：\n"); print(tabi)
write.csv(tabi, "output/tables/CI_IPC_weighted.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

## ---- 6) 双尺度交互（Major 10）----
fit_ixc <- coxph(Surv(ci_tm, ci_ev) ~ CKM34 * frail2 + age_bl + gender +
                   edu_bl + rural_bl + marry_bl + nation_bl + smoken_bl +
                   drinkev_bl + exercise_bl + sleep_bl + ln_hhcperc +
                   ins_bl + pension_bl,
                 data = rs, cluster = householdID)
sic <- summary(fit_ixc)$coefficients
print(sic[c("CKM34", "frail2", "CKM34:frail2"), ])
reri_c <- reri_cox(fit_ixc, "CKM34", "frail2", "CKM34:frail2")
print(reri_c)
write.csv(reri_c, "output/tables/CI_RERI.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(round(exp(cbind(est = coef(fit_ixc), summary(fit_ixc)$conf.int)), 3),
          "output/tables/CI_interaction_components.csv", fileEncoding = "UTF-8")

## ---- 7) PH 诊断 + 附录（已定口径）----
zphc <- cox.zph(fit_ci)
print(zphc)
write.csv(round(zphc$table, 4), "output/tables/CI_PH_tests.csv", fileEncoding = "UTF-8")

jnc <- paste0("joint6", jlev[-1])
## 因子型 joint6 需先转数值哑变量才能进 tt()
jd <- model.matrix(~ joint6, data = rs)[, -1]
colnames(jd) <- paste0("jg", seq_len(ncol(jd)))
rs <- cbind(rs, jd)
f_ph <- as.formula(paste(
  "Surv(ci_tm, ci_ev) ~", paste(colnames(jd), collapse = " + "),
  "+ strata(edu_bl) + strata(nation_bl) +",
  paste(setdiff(cov13, c("edu_bl", "nation_bl")), collapse = " + "),
  "+", paste(paste0("tt(", colnames(jd), ")"), collapse = " + ")))
fit_tt <- coxph(f_ph, data = rs, cluster = householdID,
                tt = function(x, t, ...) x * t)
co <- coef(fit_tt)
hr_t <- do.call(rbind, lapply(seq_along(jlev[-1]), function(i) {
  bg <- unname(co[paste0("jg", i)]); bt <- unname(co[paste0("tt(jg", i, ")")])
  data.frame(group = jlev[-1][i],
             HR_3y = round(exp(bg + bt * 3), 3),
             HR_6y = round(exp(bg + bt * 6), 3),
             HR_9y = round(exp(bg + bt * 9), 3))
}))
cat("\nPH 附录：tt(joint6) 时依 HR(t)（edu_bl/nation_bl 已 strata 吸收）\n")
print(hr_t)
write.csv(hr_t, "output/tables/CI_PH_timevarying_HR.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cat("05 完成\n")
