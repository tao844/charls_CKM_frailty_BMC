## ============================================================
## 08_discrimination.R —— 增量预测价值（Major 10）
## 基线模型 = CKM 分期(0~2/3~4) + 年龄 + 性别；
## 扩展模型 = 基线 + 衰弱（frailty_bl 三分类）
## 指标：Harrell C 变化（bootstrap CI）；timeROC（可选）；NRI/IDI（可选）
## 若指标不显著：按审稿人建议改写措辞（"associated with"而非"stratifies/predicts"）
## 【2026-11-04 口径】免装包方案：Harrell C ΔC + 聚类 bootstrap CI 为主；
##   timeROC/survIDINRI 未安装，自动跳过；NRI/IDI 走附录/降级路线由 PI 定夺
## ============================================================
source("00_setup.R")
deaths <- read.csv("output/deaths.csv", fileEncoding = "UTF-8")
rs <- merge(rs, deaths, by = "ID", all.x = TRUE)
rs$frailty3 <- factor(rs$frailty_bl, levels = c("无衰弱","衰弱前期","衰弱"))

fA <- Surv(dtime_yr, died) ~ CKM34 + age_bl + gender
fB <- Surv(dtime_yr, died) ~ CKM34 + age_bl + gender + frailty3
fitA <- coxph(fA, data = rs); fitB <- coxph(fB, data = rs)
CA <- concordance(fitA)$concordance
CB <- concordance(fitB)$concordance
cat("C(CKM+年龄+性别) =", round(CA, 3), "  C(+衰弱) =", round(CB, 3),
    "  ΔC =", round(CB - CA, 4), "\n")

## ΔC 的 cluster bootstrap CI（按 householdID 抽）
set.seed(1)
hh <- unique(rs$householdID)
dC <- replicate(200, {
  s <- sample(hh, replace = TRUE)
  idx <- unlist(lapply(s, function(h) which(rs$householdID == h)))
  d <- rs[idx, ]
  tryCatch(concordance(coxph(fB, data = d))$concordance -
           concordance(coxph(fA, data = d))$concordance,
           error = function(e) NA)
})
ci <- quantile(dC, c(.025,.975), na.rm = TRUE)
out <- data.frame(C_base = round(CA,3), C_ext = round(CB,3),
                  dC = round(CB-CA,4), dC_lo = round(ci[1],4),
                  dC_hi = round(ci[2],4))
write.csv(out, "output/tables/Discrimination_death.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
print(out)

## timeROC（5 年，可选）
if (requireNamespace("timeROC", quietly = TRUE)) {
  ## 线性预测子
  rs$lpA <- predict(fitA, type = "lp"); rs$lpB <- predict(fitB, type = "lp")
  rocA <- timeROC::timeROC(rs$dtime_yr, rs$died, rs$lpA, cause = 1,
                           times = 5, iid = TRUE)
  rocB <- timeROC::timeROC(rs$dtime_yr, rs$died, rs$lpB, cause = 1,
                           times = 5, iid = TRUE)
  ct <- timeROC::compare(rocA, rocB)
  print(ct); saveRDS(list(rocA = rocA, rocB = rocB, ct = ct),
                     "output/timetROC_death.rds")
}

## NRI/IDI（可选，survIDINRI）
if (requireNamespace("survIDINRI", quietly = TRUE)) {
  t0 <- 9
  z0 <- model.matrix(fA, rs)[, -1]; z1 <- model.matrix(fB, rs)[, -1]
  nri <- survIDINRI::IDI.INF(rs$dtime_yr, rs$died, z0, z1, t0 = t0, npert = 300)
  print(survIDINRI::IDI.INF.OUT(nri))
}
## 认知结局的增量价值同法（把结局换成 ci_ev/ci_tm）——可选。
cat("08 完成\n")
