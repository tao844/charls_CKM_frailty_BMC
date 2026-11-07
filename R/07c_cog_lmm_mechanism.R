## ============================================================
## 07c_cog_lmm_mechanism.R —— PI 意见四：LMM 机制检验（回归均数假设）
## 规格：a) 全认知样本（含基线障碍者） b) 高基线子组（cog_bl > 本教育层均数）
## 判读：若限制在高基线者后 衰弱-健壮 斜率差收窄 → 支持回归均数；不收窄 → 另寻解释
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
long <- rbind(mk_wave(1), mk_wave(2), mk_wave(3), mk_wave(4)) %>% filter(!is.na(cog), !is.na(t_visit))
d <- core %>% select(ID, age, gender, edu_bl, bl_time) %>%
  left_join(ex %>% select(ID, frail3), by = "ID") %>% left_join(long, by = "ID") %>%
  mutate(frail3 = factor(frail3, levels = c("robust", "prefrail", "frail"))) %>%
  filter(!is.na(frail3))
bl <- d %>% filter(wave == 1) %>% select(ID, cog_bl = cog, edu_bl)
thr <- bl %>% group_by(edu_bl) %>% summarise(threshold = mean(cog_bl) - sd(cog_bl),
                                             mean_c = mean(cog_bl), .groups = "drop")
bl <- bl %>% left_join(thr, by = "edu_bl") %>% mutate(impaired = cog_bl <= threshold)
sd0 <- sd(bl$cog_bl)  # 全样本基线 SD 作统一标准化
d <- d %>% left_join(bl %>% select(ID, cog_bl, impaired, mean_c), by = "ID") %>%
  mutate(cog_z = cog / sd0, time_yr = t_visit - bl_time) %>%
  filter(!is.na(cog_bl), !is.na(time_yr), time_yr > -0.5, time_yr < 10)

run_lmm <- function(dd, tag) {
  m <- lmer(cog_z ~ frail3 * time_yr + age + gender + factor(edu_bl) + (1 + time_yr | ID), data = dd)
  cf <- summary(m)$coefficients
  V <- vcov(m)
  diff <- cf["frail3frail:time_yr", "Estimate"]
  se <- sqrt(V["time_yr", "time_yr"] + V["frail3frail:time_yr", "frail3frail:time_yr"] + 2 * V["time_yr", "frail3frail:time_yr"])
  cat(sprintf("%-38s 人=%d 观测=%d | 衰弱-健壮斜率差=%.4f (95%%CI %.4f 至 %.4f)\n",
              tag, length(unique(dd$ID)), nrow(dd), diff, diff - 1.96 * se, diff + 1.96 * se))
}

cat("== 机制检验（统一 sd0 =", round(sd0, 2), "）==\n")
run_lmm(d, "a) 全认知样本（含基线障碍）")
run_lmm(d %>% filter(!impaired), "b) 风险集（基线>界值，=07 主分析口径）")
run_lmm(d %>% filter(cog_bl > mean_c), "c) 高基线子组（>本层均数）")
cat("\n07c 完成\n")
