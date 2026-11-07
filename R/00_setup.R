## ============================================================
## 00_setup.R —— 包、路径、统一协变量定义
## ============================================================
rm(list = ls())

need <- c("survival", "mice")
for (p in need) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p)
library(survival)

opt <- c("cmprsk", "nlme", "lme4", "timeROC", "survIDINRI", "mstate")
for (p in opt) cat(p, ":", requireNamespace(p, quietly = TRUE), "\n")

dir.create("output", showWarnings = FALSE)
dir.create("output/tables", showWarnings = FALSE)

## ---------- 数据 ----------
bl  <- readRDS("CKMx衰弱_基线数据.rds")     # 10,556 人基线
lng <- readRDS("CKMx衰弱_分析长数据.rds")   # 波次纵向
ful <- readRDS("charls_全量清洗数据.rds")   # 全量清洗（96,628 行，含各波）

## 风险集：无基线认知障碍（8,916 人）—— Table 1 主体将改用此集
rs <- bl[bl$cog_impair_bl == 0, ]
cat("全样本", nrow(bl), " 风险集", nrow(rs), "\n")

## ---------- 联合暴露组（六组，参照 = CKM 0~2 健壮）----------
jlev <- c("CKM 0~2期 无衰弱","CKM 0~2期 衰弱前期","CKM 0~2期 衰弱",
          "CKM 3~4期 无衰弱","CKM 3~4期 衰弱前期","CKM 3~4期 衰弱")
rs$joint6 <- factor(rs$joint, levels = jlev)
bl$joint6 <- factor(bl$joint, levels = jlev)

## 二分类暴露（交互/加性分析用）
rs$CKM34   <- as.integer(rs$CKM_Stage_bl >= 3)
rs$frail2  <- as.integer(rs$frailty_bl == "衰弱")            # 衰弱 vs 其余
rs$frailord <- as.integer(rs$frailty_bl == "衰弱前期") +
  2*as.integer(rs$frailty_bl == "衰弱")          # 0/1/2 趋势用

## ---------- 统一协变量集（Major 9）----------
## 主要模型（13 个，死亡与认知两个结局完全一致）：
##   社会人口学：性别、教育、城乡、婚姻、民族
##   行为/经济：吸烟、饮酒、锻炼、睡眠时长、家庭人均消费(log)、医保、养老保险
##   年龄：以随访时间为时间尺度时作为协变量；以年龄为时间尺度的敏感性分析中移入时间轴
##  DAG 定位：以上均为基线混杂因素。
##  抑郁（CES-D 10）与 BMI 可能是"暴露→结局"通路上的中介（Major 7c/9），
##  故不进主要模型，仅在敏感性分析中加入（见 04/05 第 9 节）。
cov13 <- c("gender","edu_bl","rural_bl","marry_bl","nation_bl",
           "smoken_bl","drinkev_bl","exercise_bl","sleep_bl",
           "hhcperc_bl","ins_bl","pension_bl","age_bl")
rs$ln_hhcperc <- log(rs$hhcperc_bl)
cov13[cov13 == "hhcperc_bl"] <- "ln_hhcperc"
cov13 <- c(cov13, "age_bl")  # 确保 age_bl 在列（去重）
cov13 <- unique(cov13)

## ---------- 高缺失协变量的处理（2026-11-04 数据体检发现）----------
## 缺失率：exercise_bl 57%、pension_bl 18%、nation_bl 14%、sleep_bl 0.4%。
## 若整列删除，完整案例仅 3,128/10,556（死亡 75 例）——不可接受。
## 处理：三个高缺失变量各设 "missing" 水平（missing-indicator 法），
## 主模型保留全样本；06 号多重插补作为敏感性分析互证。
## 【核对编码】exercise_bl：假设 1=锻炼、0=不锻炼；pension_bl：1=领取；nation_bl：1=少数民族。
## 如编码不同，改下方 ifelse 条件即可。
rs$exercise_bl <- factor(ifelse(is.na(rs$exercise_bl), "missing",
                                ifelse(rs$exercise_bl == 1, "yes", "no")),
                         levels = c("no", "yes", "missing"))
rs$pension_bl  <- factor(ifelse(is.na(rs$pension_bl), "missing",
                                ifelse(rs$pension_bl == 1, "yes", "no")),
                         levels = c("no", "yes", "missing"))
rs$nation_bl   <- factor(ifelse(is.na(rs$nation_bl), "missing",
                                ifelse(rs$nation_bl == 1, "minority", "han")),
                         levels = c("han", "minority", "missing"))

cat("协变量缺失情况（风险集）：\n")
print(sapply(rs[, cov13], function(x) sum(is.na(x))))

## ---------- 通用小工具 ----------
## 对比（线性组合）估计：HR、95%CI、P 值 —— 用于组内对比（Major 2）
contrast_cox <- function(fit, L, label = "") {
  b <- coef(fit); V <- vcov(fit)
  est <- sum(L * b); se <- sqrt(drop(t(L) %*% V %*% L))
  z <- est / se
  data.frame(contrast = label,
             HR = round(exp(est), 3),
             lo = round(exp(est - 1.96*se), 3),
             hi = round(exp(est + 1.96*se), 3),
             p = signif(2*pnorm(-abs(z)), 3))
}

## RERI 加性交互（delta 法 CI）—— Major 10
reri_cox <- function(fit, nmA, nmB, nmAB) {
  b <- coef(fit); V <- vcov(fit)
  a <- b[nmA]; bb <- b[nmB]; c <- b[nmAB]
  g <- c(-exp(a), -exp(bb), exp(c))
  V3 <- V[c(nmA, nmB, nmAB), c(nmA, nmB, nmAB)]
  se <- sqrt(drop(t(g) %*% V3 %*% g))
  R <- exp(c) - exp(a) - exp(bb) + 1
  data.frame(RERI = round(R, 3),
             lo = round(R - 1.96*se, 3),
             hi = round(R + 1.96*se, 3))
}

sink("output/sessionInfo.txt")
print(sessionInfo())
sink()
cat("00 完成\n")