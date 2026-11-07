# 10_timefixed_supplement.R
# T1（time-fixed 尺度）joint6 六格：全队列与剔除 2011-13 段两版
# 产物：投稿包 Supplementary_YMJ.docx S4 表；附 1,074 人排除子组死亡率核验
# 复核记录：与 PI 独立复算 v2 表逐位一致（2026-11-07）
suppressMessages({library(dplyr); library(survival)})
proj <- "D:/00 陶万清的工作/07 working stata/01 CKM_frailty_revision"
setwd(proj)
s2 <- readRDS("output/rebuild/bl_cov_step2.rds")
s3 <- readRDS("output/rebuild/exposures_step3.rds")
dd <- s2 %>% inner_join(s3, by = "ID")
covs <- c("age","gender","edu_bl","rural_bl","marry_bl","smoken_bl","drinkev_bl",
          "sleep_bl","ln_hhcperc","ins_bl","pension_bl")
fiok <- !is.na(dd$FI41) & !is.na(dd$CKM2)
cc <- complete.cases(dd[, covs])
d_model <- dd[fiok & cc, ]
cat("model n =", nrow(d_model), " deaths =", sum(d_model$died), "\n")
d_model$joint6 <- factor(d_model$joint6)
d_model$joint6 <- relevel(d_model$joint6, ref = "0-2_robust")

f <- as.formula(paste("Surv(dtime_yr, died) ~ joint6 +", paste(covs, collapse = " + ")))
m1 <- coxph(f, data = d_model, cluster = householdID)
cf <- summary(m1)$conf.int
cat("=== T1 joint6 (含2011-13) ===\n")
print(round(cf[rownames(cf) %in% grep("joint6", rownames(cf), value = TRUE), c("exp(coef)","lower .95","upper .95")], 3))

# v2: 剔除 2011-13 段死亡（dtime_yr <= 2）
d2 <- d_model[!(d_model$died == 1 & d_model$dtime_yr <= 2), ]
cat("\nv2 n =", nrow(d2), " deaths =", sum(d2$died), "\n")
m2 <- coxph(f, data = d2, cluster = householdID)
cf2 <- summary(m2)$conf.int
cat("=== T1 joint6 (剔除2011-13) ===\n")
print(round(cf2[rownames(cf2) %in% grep("joint6", rownames(cf2), value = TRUE), c("exp(coef)","lower .95","upper .95")], 3))

# A2 子组核验：572 协变量不完整 vs 502 暴露不可判（在 13,256 分析队列内）
analytic <- dd  # s2 已是分析队列？
cat("\ns2 n =", nrow(s2), "\n")
g572 <- dd[fiok & !cc, ]; g502 <- dd[!fiok, ]
cat("572组: n=", nrow(g572), " 死亡率=", round(mean(g572$died, na.rm=TRUE)*100,1), "%\n")
cat("502组: n=", nrow(g502), " 死亡率=", round(mean(g502$died, na.rm=TRUE)*100,1), "%\n")
