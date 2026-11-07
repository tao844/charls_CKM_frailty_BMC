## ============================================================
## 07_cognition_lmm.R —— 认知 LMM 最小检验（PI 裁定口径）
## 风险集：基线无认知障碍（教育分层界值 = 各教育层 wave1 认知均数 −1SD，界值表输出）
## 时点：同一套真访谈年月（r*iwy/iwm），time = 访视时点 − 基线访视时点（年）
## 暴露：FI 三分类（frail3，重建 41 项口径）；调整 age gender edu_bl；随机截距+随机斜率
## 结局：cog_h = (imrc+dlrc)/2 + ser7 + orient + draw（wave1-4 同构），按基线 SD 标准化
## 判据：衰弱组斜率差（frail vs robust）< 0.02 SD/年 → 关闭该线
## ============================================================
suppressMessages(library(haven)); suppressMessages(library(dplyr)); suppressMessages(library(lme4))
proj <- "D:/00 陶万清的工作/07 working stata/01 CKM_frailty_revision"
core <- readRDS(file.path(proj, "output/rebuild/bl_cov_step2.rds"))
ex   <- readRDS(file.path(proj, "output/rebuild/exposures_step3.rds"))

hh <- read_dta(file.path(proj, "CHARLS_official/H_CHARLS_D_Data (1)/H_CHARLS_D_Data.dta"))
hh$ID <- sprintf("%012.0f", as.numeric(as.character(hh$ID)))

mk_wave <- function(i) {
  data.frame(ID = hh$ID,
             cog = (as.numeric(hh[[paste0("r", i, "imrc")]]) + as.numeric(hh[[paste0("r", i, "dlrc")]])) / 2 +
                    as.numeric(hh[[paste0("r", i, "ser7")]]) + as.numeric(hh[[paste0("r", i, "orient")]]) +
                    as.numeric(hh[[paste0("r", i, "draw")]]),
             t_visit = as.numeric(hh[[paste0("r", i, "iwy")]]) + (as.numeric(hh[[paste0("r", i, "iwm")]]) - 0.5) / 12,
             wave = i)
}
long <- rbind(mk_wave(1), mk_wave(2), mk_wave(3), mk_wave(4)) %>%
  filter(!is.na(cog), !is.na(t_visit))

d <- core %>% select(ID, age, gender, edu_bl, bl_time) %>%
  left_join(ex %>% select(ID, frail3), by = "ID") %>%
  left_join(long, by = "ID") %>%
  mutate(frail3 = factor(frail3, levels = c("robust", "prefrail", "frail"))) %>%
  filter(!is.na(frail3))

## 基线认知与界值（wave1，真访谈时点口径）
bl <- d %>% filter(wave == 1) %>% select(ID, cog_bl = cog, edu_bl)
thr <- bl %>% group_by(edu_bl) %>%
  summarise(n = n(), mean = mean(cog_bl), sd = sd(cog_bl), threshold = mean - sd, .groups = "drop")
cat("== 教育分层界值表（wave1 cog_h 均数 −1SD）==\n")
print(round(thr, 3))
write.csv(round(thr, 3), file.path(proj, "output/rebuild/cog_thresholds_by_edu.csv"), row.names = FALSE, fileEncoding = "UTF-8")

bl <- bl %>% left_join(thr %>% select(edu_bl, threshold), by = "edu_bl") %>%
  mutate(impaired = cog_bl <= threshold)
cat("\n基线认知障碍（教育分层界值）:", sum(bl$impaired, na.rm = TRUE), "/", nrow(bl), "\n")

risk_ids <- bl$ID[!bl$impaired & !is.na(bl$impaired)]
cat("风险集:", length(risk_ids), "\n")

## 标准化（wave1 风险集 SD）
sd0 <- sd(bl$cog_bl[bl$ID %in% risk_ids])
d <- d %>% filter(ID %in% risk_ids) %>%
  left_join(bl %>% select(ID, cog_bl), by = "ID") %>%
  mutate(cog_z = cog / sd0,
         time_yr = t_visit - bl_time) %>%
  filter(!is.na(time_yr), time_yr > -0.5, time_yr < 10)

cat("\n风险集 LMM 数据: 人", length(unique(d$ID)), " 观测", nrow(d), "\n")
print(d %>% group_by(wave, frail3) %>% summarise(n = n(), mean_z = round(mean(cog_z), 3), .groups = "drop"))

## LMM：随机截距+随机斜率
m <- lmer(cog_z ~ frail3 * time_yr + age + gender + factor(edu_bl) + (1 + time_yr | ID), data = d)
cf <- summary(m)$coefficients
cat("\n== LMM 固定效应（cog_z / 年）==\n")
print(round(cf[rownames(cf) %in% c("frail3prefrail", "frail3frail", "time_yr",
      "frail3prefrail:time_yr", "frail3frail:time_yr"), ], 4))

## 各组斜率与判据
b_r  <- cf["time_yr", "Estimate"]
b_p  <- b_r + cf["frail3prefrail:time_yr", "Estimate"]
b_f  <- b_r + cf["frail3frail:time_yr", "Estimate"]
se_f <- sqrt(vcov(m)["time_yr", "time_yr"] + vcov(m)["frail3frail:time_yr", "frail3frail:time_yr"] +
             2 * vcov(m)["time_yr", "frail3frail:time_yr"])
diff <- cf["frail3frail:time_yr", "Estimate"]
cat("\n== 判据 ==\n")
cat("斜率（SD/年）: 健壮", round(b_r, 4), " 衰前", round(b_p, 4), " 衰弱", round(b_f, 4), "\n")
cat("衰弱 vs 健壮 斜率差:", round(diff, 4), " (SE", round(se_f, 4), ")  95%CI:",
    round(diff - 1.96 * se_f, 4), "至", round(diff + 1.96 * se_f, 4), "\n")
cat("判据 |差| < 0.02 SD/年 :", abs(diff) < 0.02, "\n")

saveRDS(m, file.path(proj, "output/rebuild/_cog_lmm.rds"))
write.csv(round(cf, 4), file.path(proj, "output/rebuild/cog_lmm_fixed_effects.csv"), fileEncoding = "UTF-8")
cat("\n07 认知 LMM 完成\n")
