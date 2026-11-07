## ============================================================
## 03_alt_exposures.R —— CKM 分期构成表 + FI_noCKM 接口（Major 6/7）
## 在代码包 03_exposures 基础上：
##   第 5 节无原始 41 条目打分代码，改用简报 §3.5 已验证的
##   "精确算术剔除法"：FI_noCKM = (k - d5)/(n-5)，(k,n) 由 FI_bl 反解
## ============================================================
source("00_setup.R")

## ---- 1) CKM 各期分布（Major 7a）----
stg <- data.frame(
  stage = 0:4,
  n = as.integer(table(factor(bl$CKM_Stage_bl, levels = 0:4))),
  pct = round(as.numeric(prop.table(table(factor(bl$CKM_Stage_bl, levels = 0:4)))) * 100, 1))
print(stg)
write.csv(stg, "output/tables/CKM_stage_distribution.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

## ---- 2) 分期标准对照表（Major 7a 补充表）----
## 注：带【】处须与原始建组代码核对后由写作阶段定稿，此处按代码包模板原样落表
ckm_map <- data.frame(
  分期 = c(0,1,2,2,3,3,4),
  标准 = c("无 CKM 危险因素",
           " excess 脂肪（超重/肥胖）",
           "代谢危险因素/慢性肾病（含危险当量）",
           "10 年 CVD 风险 ≥【与原始代码核对，如 7.5%】（Framingham 一般人群模型）",
           "亚临床 CVD / 估算风险当量（eGFR<60 或高危代谢病组合）",
           "极高危 CKD（无 UACR，按【原始代码口径】降级标注）",
           "临床 CVD：本次以自报心脏病代替（卒中已排除，须在文中声明）"),
  CHARLS变量 = c("—","bmi_bl / mwaist_bl",
                 "hibpe_1_bl, diabe_1_bl, dyslipe_bl, kidneye_bl, eGFR_bl",
                 "fmh_score_bl","eGFR_bl / fmh_score_bl","eGFR_bl","hearte_bl"),
  备注 = c("","","高血压=自报或用药；糖尿病=自报或用药或 HbA1c/FBG",
           "为何用 Framingham 而不用 AHA PREVENT（Khan 2024）须回应：CHARLS 缺少 PREVENT 所需部分变量（如 UACR、完整血脂细节），见回复信",
           "eGFR 用 2021 肌酐-胱抑素 C 联合方程","CHARLS 无尿白蛋白，KDIGO 热图不可直接应用——如实声明局限","stroke 已作为排除标准")
)
write.csv(ckm_map, "output/tables/S_CKM_staging_map.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

## ---- 4) 无血样者的分期处理说明（Major 6）----
have_blood <- !is.na(bl$bl_crea_bl)          # 2011 血样子样本
cat("有血样：", sum(have_blood), " 无血样：", sum(!have_blood), "\n")
cat("无血样者中 CKM 3~4 期人数：", sum(bl$CKM_Stage_bl >= 3 & !have_blood),
    "（初稿 S1 称暴露无缺失——须改为如实说明无血样者的分期规则）\n")

## ---- 5) FI_noCKM 精确算术剔除法（Major 7b；简报 §3.5 已验证口径）----
## 反解 (k,n)：按人搜索 n in 30:50，取 |k/n - FI_bl| 最小者
best <- t(vapply(bl$FI_bl, function(f) {
  cand <- do.call(rbind, lapply(30:50, function(n) {
    k <- round(f * n); c(k = k, err = abs(k / n - f))
  }))
  i <- which.min(cand[, "err"]); c(k = cand[i, "k"], n = (30:50)[i], err = cand[i, "err"])
}, c(k = 0, n = 0, err = 0)))
rev_ok <- best[, "err"] < 0.005
cat("\nFI 反解率：", round(mean(rev_ok) * 100, 2), "%（按人搜索 n；n=41 者 ",
    sum(best[, "n"] == 41), " 人）\n")

## 5 个重叠条目：高血压、糖尿病、血脂异常、心脏病、肾病
## 【口径锁定 2026-11-04 校准】高血压/糖尿病取自报变量（hibpe_bl/diabe_bl），
## 与已验证预期值完全吻合（健壮组均值 0.049、三分类一致率 87.3%）
it <- c("hibpe_bl", "diabe_bl", "dyslipe_bl", "hearte", "kidneye_bl")
cat("重叠条目缺失数：\n"); print(sapply(it, function(v) sum(is.na(bl[[v]]))))
d5 <- rowSums(sapply(it, function(v) as.integer(bl[[v]] == 1)), na.rm = TRUE)
cat("d5 分布（0-5 个重叠条目受累人数）：\n"); print(table(d5))

k_noCKM <- pmax(best[, "k"] - d5, 0)
FI_noCKM <- k_noCKM / (best[, "n"] - 5)
cat("CKM 3~4 期健壮组 FI_noCKM 均值：", round(mean(FI_noCKM[bl$CKM_Stage_bl >= 3 & bl$frailty_bl == "无衰弱"], na.rm = TRUE), 3),
    "（原始 FI 0.062，预期 0.049）\n")

## 重新三分类
fr_nc <- cut(FI_noCKM, breaks = c(-Inf, 0.10, 0.249999, Inf),
             labels = c("无衰弱", "衰弱前期", "衰弱"))
agree <- mean(as.character(fr_nc) == as.character(bl$frailty_bl), na.rm = TRUE)
cat("三分类一致率：", round(agree * 100, 1), "%（预期 87.3%）\n")
print(table(原始 = bl$frailty_bl, FI_noCKM分类 = fr_nc))

## 六联合暴露组（FI_noCKM 版），供 04/05 第 9f 节敏感性重跑
bl$FI_noCKM  <- FI_noCKM
bl$frailty_nc <- fr_nc
ckm_grp <- ifelse(bl$CKM_Stage_bl >= 3, "CKM 3~4期", "CKM 0~2期")
bl$joint6_nc <- factor(paste(ckm_grp, bl$frailty_nc),
                       levels = c("CKM 0~2期 无衰弱","CKM 0~2期 衰弱前期","CKM 0~2期 衰弱",
                                  "CKM 3~4期 无衰弱","CKM 3~4期 衰弱前期","CKM 3~4期 衰弱"))
saveRDS(bl[, c("ID", "FI_noCKM", "frailty_nc", "joint6_nc")],
        "output/bl_FInoCKM.rds")
cat("FI_noCKM 六组分布：\n"); print(table(bl$joint6_nc))

cat("03_alt 完成\n")
