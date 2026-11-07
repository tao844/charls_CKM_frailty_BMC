## ============================================================
## 09b_checks.R —— 阶段 3 杂项核查（审稿意见→证据映射表取证）
##   1) Major 6：无血样者的 CKM 分期规则核查（汇总级，无个体行输出）
##   2) Minor 9：连续 FI（原始 FI_bl 与 FI_noCKM）的死亡敏感性分析
## ============================================================
source("00_setup.R")
deaths <- read.csv("output/deaths.csv", fileEncoding = "UTF-8")

## ---- 1) Major 6：无血样者分期规则核查 ----
## 血样（2011 采集）可用性代理：肌酐 bl_crea_bl 非缺失
have_blood <- !is.na(bl$bl_crea_bl)
cat("有血样：", sum(have_blood), "  无血样：", sum(!have_blood), "\n")

## 各分期所需关键变量的缺失结构（按有无血样分层，仅汇总计数）
key_vars <- c("bl_glu_bl", "bl_hbalc_bl", "bl_crea_bl", "bl_cysc_bl",
              "bl_cho_bl", "bl_hdl_bl", "bl_tg_bl", "bl_hdl_bl")
miss_tab <- do.call(rbind, lapply(key_vars, function(v) {
  data.frame(var = v,
             miss_noblood = sum(is.na(bl[[v]])[!have_blood]),
             miss_blood   = sum(is.na(bl[[v]])[ have_blood]))
}))
print(miss_tab)

## 无血样者的分期分布 vs 有血样者
stg_by_blood <- as.data.frame.matrix(
  table(factor(bl$CKM_Stage_bl, levels = 0:4), blood = have_blood))
stg_by_blood$stage <- 0:4
stg_by_blood <- stg_by_blood[, c("stage", "FALSE", "TRUE")]
names(stg_by_blood)[2:3] <- c("no_blood_n", "blood_n")
cat("\n分期 × 血样可用性：\n"); print(stg_by_blood)

## 无血样者进入 3~4 期的组分来源核查：
## 3 期可由 eGFR<60 或 fmh 高危派生；无血样者 eGFR/fmh 应缺失——
## 若仍有 3~4 期，只能来自自报/体检派生组分（hearte/kidneye/血压/腰围/BMI）
nb <- bl[!have_blood, ]
nb34 <- nb[nb$CKM_Stage_bl >= 3, ]
comp <- data.frame(
  n_noblood_3  = sum(nb$CKM_Stage_bl == 3),
  n_noblood_4  = sum(nb$CKM_Stage_bl == 4),
  hearte_self  = sum(nb34$hearte == 1, na.rm = TRUE),          # 4 期自报心脏病
  kidneye_self = sum(nb34$kidneye_bl == 1, na.rm = TRUE),      # 自报肾病（3 期等价组分）
  hibpe_tx     = sum(nb34$hibpe_1_bl == 1, na.rm = TRUE),      # 高血压自报或用药
  diabe_tx     = sum(nb34$diabe_1_bl == 1, na.rm = TRUE),      # 糖尿病自报或用药
  dyslipe_self = sum(nb34$dyslipe_bl == 1, na.rm = TRUE),      # 血脂异常自报
  eGFR_avail   = sum(!is.na(nb34$eGFR_bl)),
  fmh_avail    = sum(!is.na(nb34$fmh_score_bl)))
cat("\n无血样且 3~4 期者的组分来源（汇总）：\n"); print(comp)
write.csv(comp, "output/tables/CKM_noblood_components.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

## 规则结论表（逐条）
rules <- data.frame(
  item = c("无血样者总人数", "无血样者中 CKM 3~4 期人数",
           "3~4 期中 eGFR 可得人数", "3~4 期中 fmh_score 可得人数",
           "3~4 期可解释人数（自报心脏病/肾病任一=1）",
           "派生规则结论"),
  value = c(sum(!have_blood),
            nrow(nb34),
            comp$eGFR_avail,
            comp$fmh_avail,
            sum(nb34$hearte == 1 | nb34$kidneye_bl == 1, na.rm = TRUE),
            "无血样者（n=2,978）无法进入需血样的组分（eGFR/HbA1c/FBG/血脂/fmh）；其分期仅由自报与体检测量派生：4 期=自报心脏病 hearte（n=314）；3 期仅 2 人，来自自报肾病/高危代谢病组合（eGFR 与 fmh 均不可得）；0~2 期由 BMI/腰围/血压/自报代谢病派生。初稿'暴露无缺失'表述作废，按此规则如实报告"),
  stringsAsFactors = FALSE)
write.csv(rules, "output/tables/CKM_noblood_staging_rule.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

## ---- 2) Minor 9：连续 FI 的死亡敏感性（补齐认知端已有分析）----
d_all <- merge(bl, deaths, by = "ID", all.x = TRUE)
stopifnot(!any(is.na(d_all$died)))
finc <- read.csv("output/tables/FI_noCKM.csv", fileEncoding = "UTF-8")
d_all$FI_noCKM <- finc$FI_noCKM[match(d_all$ID, finc$ID)]
## 与 04 相同的协变量预处理（bl 侧为原始三分类因子版本）
d_all$exercise_bl <- factor(ifelse(is.na(d_all$exercise_bl), "missing",
  ifelse(d_all$exercise_bl == 1, "yes", "no")), levels = c("no", "yes", "missing"))
d_all$pension_bl <- factor(ifelse(is.na(d_all$pension_bl), "missing",
  ifelse(d_all$pension_bl == 1, "yes", "no")), levels = c("no", "yes", "missing"))
d_all$nation_bl <- factor(ifelse(is.na(d_all$nation_bl), "missing",
  ifelse(d_all$nation_bl == 1, "minority", "han")), levels = c("han", "minority", "missing"))
d_all$CKM34 <- as.integer(d_all$CKM_Stage_bl >= 3)
d_all$ln_hhcperc <- log(d_all$hhcperc_bl)

cc13 <- paste(cov13, collapse = " + ")
res <- list()
for (v in c("FI_bl", "FI_noCKM")) {
  z <- scale(d_all[[v]])  # 每 1 SD
  d_all$FI_sd <- as.numeric(z)
  f1 <- coxph(as.formula(paste("Surv(dtime_yr, died) ~ CKM34 + FI_sd +", cc13)),
              data = d_all, cluster = householdID)
  f2 <- coxph(as.formula(paste("Surv(dtime_yr, died) ~ CKM34 * FI_sd +", cc13)),
              data = d_all, cluster = householdID)
  s1 <- summary(f1); s2 <- summary(f2)
  res[[v]] <- rbind(
    data.frame(FI = v, term = "CKM 3~4期", HR = round(exp(coef(f1)["CKM34"]), 3),
               lo = round(s1$conf.int["CKM34", "lower .95"], 3),
               hi = round(s1$conf.int["CKM34", "upper .95"], 3),
               p = signif(s1$coefficients["CKM34", "Pr(>|z|)"], 3)),
    data.frame(FI = v, term = "FI 每 1 SD", HR = round(exp(coef(f1)["FI_sd"]), 3),
               lo = round(s1$conf.int["FI_sd", "lower .95"], 3),
               hi = round(s1$conf.int["FI_sd", "upper .95"], 3),
               p = signif(s1$coefficients["FI_sd", "Pr(>|z|)"], 3)),
    data.frame(FI = v, term = "交互 CKM34×FI", HR = round(exp(coef(f2)["CKM34:FI_sd"]), 3),
               lo = round(s2$conf.int["CKM34:FI_sd", "lower .95"], 3),
               hi = round(s2$conf.int["CKM34:FI_sd", "upper .95"], 3),
               p = signif(s2$coefficients["CKM34:FI_sd", "Pr(>|z|)"], 3)))
}
out <- do.call(rbind, res)
print(out)
write.csv(out, "output/tables/Death_FI_continuous.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("09b 完成\n")
