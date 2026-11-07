## ============================================================
## 07_lmm.R —— 连续认知轨迹线性混合模型（Major 8c）
## 与离散认知障碍结局互补；若有需要可升级为联合模型（JM/JMbayes）
## ============================================================
source("00_setup.R")
deaths <- read.csv("output/deaths.csv", fileEncoding = "UTF-8")

lk <- lng[lng$ID %in% rs$ID & !is.na(lng$total_cognition), ]
## 修复 2026-11-04：lng 自带 gender/edu_bl/age_bl，直接 merge 会产生 .x/.y
## 后缀导致模型变量找不到；先剔除 lng 侧同名列，统一使用 rs（分析风险集）口径
lk <- lk[, !names(lk) %in% c("gender", "edu_bl", "age_bl")]
lk <- merge(lk, rs[, c("ID","joint6","gender","edu_bl","age_bl")], by = "ID")
lk$jt <- as.numeric(factor(lk$joint6, levels = jlev))

if (requireNamespace("nlme", quietly = TRUE)) {
  library(nlme)
  m_lmm <- lme(total_cognition ~ joint6 * time + gender + edu_bl,
               random = ~ time | ID, data = lk, na.action = na.omit)
  s <- summary(m_lmm)$tTable
  write.csv(round(s, 4), "output/tables/LMM_joint_time.csv", fileEncoding = "UTF-8")
  print(s[grep("joint6|time", rownames(s)), ])
  ## 各组年均下降速度 = 固定 time 主效应 + 交互项
  gd <- data.frame(group = jlev,
                   slope = round(fixef(m_lmm)["time"] +
                     c(0, fixef(m_lmm)[grep(":time", names(fixef(m_lmm)))]), 3))
  write.csv(gd, "output/tables/LMM_group_slopes.csv",
            row.names = FALSE, fileEncoding = "UTF-8")
  print(gd)
}
## 报告要点（写入返修稿）：组×时间交互 = 各联合组认知下降速度差异；
## 与"界值定义"结论互相印证；如需 joint model 再议。
cat("07 完成\n")
