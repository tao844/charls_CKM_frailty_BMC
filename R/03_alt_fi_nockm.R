## ============================================================
## 03_alt_fi_nockm.R —— FI_noCKM 精确算术剔除法（Major 7b）
## 方法（交接简报 §3.5 已验证口径）：
##   1) 由 FI_bl 反解 (k,n)：按人搜索 n in 30:50，取 |k/n - FI_bl| 最小者
##   2) 5 个重叠条目（高血压、糖尿病、血脂异常、心脏病、肾病）受累数 d5
##   3) FI_noCKM = (k - d5)/(n - 5)
##   4) 按 <0.10 / 0.10-0.24 / >=0.25 重切三分类
## 【口径锁定 2026-11-04 校准】高血压/糖尿病取自报变量（hibpe_bl/diabe_bl），
##   与已验证预期值完全吻合（健壮组均值 0.049、一致率 87.3%）
## 过关标准：
##   1 反解率 >= 99.9%；2 三分类一致率 85%-90%；
##   3 CKM3-4 健壮组 FI_noCKM 均值约 0.049；4 写出 FI_noCKM.csv
## ============================================================
source("00_setup.R")

## ---- 1) 反解 (k,n) ----
best <- t(vapply(bl$FI_bl, function(f) {
  cand <- do.call(rbind, lapply(30:50, function(n) {
    k <- round(f * n); c(k = k, err = abs(k / n - f))
  }))
  i <- which.min(cand[, "err"]); c(k = cand[i, "k"], n = (30:50)[i], err = cand[i, "err"])
}, c(k = 0, n = 0, err = 0)))
rev_rate <- mean(best[, "err"] < 0.005)
cat("FI (k,n) 反解成功率：", round(rev_rate * 100, 2), "%（n=41 者 ",
    sum(best[, "n"] == 41), " 人）\n")
if (rev_rate < 0.999) stop("反解率 ", round(rev_rate * 100, 2), "% < 99.9%，停下报告")

## ---- 2) 5 个重叠条目受累数 d5 ----
it <- c("hibpe_bl", "diabe_bl", "dyslipe_bl", "hearte", "kidneye_bl")
cat("重叠条目缺失数：\n"); print(sapply(it, function(v) sum(is.na(bl[[v]]))))
d5 <- rowSums(sapply(it, function(v) as.integer(bl[[v]] == 1)), na.rm = TRUE)
cat("d5 分布（0-5 个重叠条目受累人数）：\n"); print(table(d5))

## ---- 3) FI_noCKM 与重分类 ----
FI_noCKM <- pmax(best[, "k"] - d5, 0) / (best[, "n"] - 5)
m34rob <- mean(FI_noCKM[bl$CKM_Stage_bl >= 3 & bl$frailty_bl == "无衰弱"], na.rm = TRUE)
cat("CKM 3~4 期健壮组：原 FI 均值 ",
    round(mean(bl$FI_bl[bl$CKM_Stage_bl >= 3 & bl$frailty_bl == "无衰弱"], na.rm = TRUE), 3),
    " -> FI_noCKM 均值 ", round(m34rob, 3), "（预期约 0.049）\n")
if (abs(m34rob - 0.049) > 0.005) stop("CKM3-4 健壮组 FI_noCKM 均值 ", round(m34rob, 3),
                                      " 偏离 0.049 超容差，停下报告")

frail3_nc <- cut(FI_noCKM, breaks = c(-Inf, 0.10, 0.249999, Inf),
                 labels = c("无衰弱", "衰弱前期", "衰弱"))
agree <- mean(as.character(frail3_nc) == as.character(bl$frailty_bl), na.rm = TRUE)
cat("FI 三分类一致率：", round(agree * 100, 1), "%（应在 85%-90%）\n")
print(table(原始 = bl$frailty_bl, FI_noCKM分类 = frail3_nc))
if (agree < 0.85 || agree > 0.90) stop("三分类一致率 ", round(agree * 100, 1),
                                       "% 不在 85%-90%，停下报告")

## ---- 4) 六联合暴露组（FI_noCKM 版）与产物 ----
ckm_grp <- ifelse(bl$CKM_Stage_bl >= 3, "CKM 3~4期", "CKM 0~2期")
joint6_nc <- factor(paste(ckm_grp, frail3_nc),
                    levels = c("CKM 0~2期 无衰弱", "CKM 0~2期 衰弱前期", "CKM 0~2期 衰弱",
                               "CKM 3~4期 无衰弱", "CKM 3~4期 衰弱前期", "CKM 3~4期 衰弱"))
cat("FI_noCKM 六组分布：\n"); print(table(joint6_nc))

out <- data.frame(ID = bl$ID, FI_noCKM = round(FI_noCKM, 4),
                  frail3_nc = as.character(frail3_nc),
                  joint6_nc = as.character(joint6_nc))
write.csv(out, "output/tables/FI_noCKM.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("FI_noCKM.csv 写出：", nrow(out), "行\n")
cat("03_alt_fi_nockm 完成\n")
