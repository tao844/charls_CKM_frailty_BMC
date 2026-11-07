## ============================================================
## 06b_mi_pooled_compare.R —— MI pooled 与完整案例主模型并排对照
## 重点：exercise_bl 三分类在插补后是否保持方向与量级
## 完整案例定义：missing-indicator 三分类（missing 作独立水平，主模型口径）
## MI 口径：pmm 插补后 exercise/pension/nation 为二分类（无 missing 水平）
## ============================================================
source("00_setup.R")
deaths <- read.csv("output/deaths.csv", fileEncoding = "UTF-8")

## ---- 完整案例模型 1：死亡（全队列；同 04）----
bl$CKM34 <- as.integer(bl$CKM_Stage_bl >= 3)
bl$ln_hhcperc <- log(bl$hhcperc_bl)
bl$exercise_bl <- factor(ifelse(is.na(bl$exercise_bl), "missing",
  ifelse(bl$exercise_bl == 1, "yes", "no")), levels = c("no", "yes", "missing"))
bl$pension_bl <- factor(ifelse(is.na(bl$pension_bl), "missing",
  ifelse(bl$pension_bl == 1, "yes", "no")), levels = c("no", "yes", "missing"))
bl$nation_bl <- factor(ifelse(is.na(bl$nation_bl), "missing",
  ifelse(bl$nation_bl == 1, "minority", "han")), levels = c("han", "minority", "missing"))
d_all <- merge(bl, deaths, by = "ID", all.x = TRUE)
f_cc <- as.formula(paste("Surv(dtime_yr, died) ~ joint6 +",
                         paste(cov13, collapse = " + ")))
cc_d <- coxph(f_cc, data = d_all, cluster = householdID)
s_d <- summary(cc_d)

## ---- 完整案例模型 2：认知（风险集，教育分层界值主定义；同 05）----
rs <- merge(rs, deaths, by = "ID", all.x = TRUE)
fu <- aggregate(time ~ ID, data = lng, FUN = max)
rs$lastfu <- fu$time[match(rs$ID, fu$ID)]
lci <- lng[order(lng$ID, lng$time), ]
lci <- lci[lci$ID %in% rs$ID & !is.na(lci$total_cognition), ]
w1c <- lng[lng$wave == 1 & !is.na(lng$total_cognition) & !is.na(lng$edu_bl), ]
cut_tab <- do.call(rbind, lapply(sort(unique(w1c$edu_bl)), function(e) {
  x <- w1c$total_cognition[w1c$edu_bl == e]
  data.frame(edu_bl = e, cutoff = mean(x) - sd(x))
}))
lci$cutval <- cut_tab$cutoff[match(lci$edu_bl, cut_tab$edu_bl)]
lci$ci_e <- as.integer(!is.na(lci$cutval) & lci$total_cognition < lci$cutval)
first_e <- lci[lci$ci_e == 1, ]; first_e <- first_e[!duplicated(first_e$ID), c("ID", "time")]
names(first_e) <- c("ID", "ctime_e")
rs <- merge(rs, first_e, by = "ID", all.x = TRUE)
rs$ci_ev <- as.integer(!is.na(rs$ctime_e))
rs$ci_tm <- ifelse(rs$ci_ev == 1, rs$ctime_e,
                   ifelse(rs$died == 1, pmin(rs$lastfu, rs$dtime_yr, na.rm = TRUE),
                          rs$lastfu))
cc_c <- coxph(update(f_cc, Surv(ci_tm, ci_ev) ~ .), data = rs, cluster = householdID)
s_c <- summary(cc_c)

cc_row <- function(s, pat) {
  g <- grep(pat, rownames(s$coefficients), value = TRUE)
  data.frame(term = g,
             HR_CC = round(exp(s$coefficients[g, "coef"]), 3),
             lo_CC = round(s$conf.int[g, "lower .95"], 3),
             hi_CC = round(s$conf.int[g, "upper .95"], 3),
             p_CC = signif(s$coefficients[g, "Pr(>|z|)"], 3))
}

## ---- 并排对照 ----
mi_d <- read.csv("output/tables/Death_MI20_pooled_full.csv", fileEncoding = "UTF-8")
mi_c <- read.csv("output/tables/CI_MI20_pooled.csv", fileEncoding = "UTF-8")

cmp_d <- merge(cc_row(s_d, "^joint6"), mi_d[mi_d$term %in% cc_row(s_d, "^joint6")$term,
              c("term", "HR", "lo", "hi", "p.value")],
              by = "term", all = TRUE)
names(cmp_d)[6:9] <- c("HR_MI", "lo_MI", "hi_MI", "p_MI")
## exercise 对照：CC 三分类（yes/missing）vs MI 二分类（yes）——手工对齐，避免笛卡尔积
ex_compare <- function(cc_ex, mi_df) {
  mi_ex <- mi_df[grep("exercise_f", mi_df$term), c("HR", "lo", "hi", "p.value")]
  pad1 <- function(v) c(v, NA)[seq_len(nrow(cc_ex))]  # MI 仅 1 行（yes），缺失水平补 NA
  data.frame(term   = cc_ex$term,
             HR_CC = cc_ex$HR_CC, lo_CC = cc_ex$lo_CC, hi_CC = cc_ex$hi_CC, p_CC = cc_ex$p_CC,
             HR_MI = pad1(mi_ex$HR), lo_MI = pad1(mi_ex$lo),
             hi_MI = pad1(mi_ex$hi),  p_MI = pad1(mi_ex$p.value))
}
cmp_d_ex <- ex_compare(cc_row(s_d, "exercise_bl"), mi_d)
cmp_c <- merge(cc_row(s_c, "^joint6"), mi_c[mi_c$term %in% cc_row(s_c, "^joint6")$term,
              c("term", "HR", "lo", "hi", "p.value")],
              by = "term", all = TRUE)
names(cmp_c)[6:9] <- c("HR_MI", "lo_MI", "hi_MI", "p_MI")
cmp_c_ex <- ex_compare(cc_row(s_c, "exercise_bl"), mi_c)

cat("\n== 死亡（全队列）：CC vs MI pooled ==\n"); print(cmp_d)
cat("\n-- exercise_bl（CC 三分类 vs MI 二分类）--\n"); print(cmp_d_ex)
cat("\n== 认知（风险集）：CC vs MI pooled ==\n"); print(cmp_c)
cat("\n-- exercise_bl（CC 三分类 vs MI 二分类）--\n"); print(cmp_c_ex)

write.csv(rbind(cmp_d, cmp_d_ex), "output/tables/Death_CC_vs_MI20.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(rbind(cmp_c, cmp_c_ex), "output/tables/CI_CC_vs_MI20.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("06b 完成\n")
