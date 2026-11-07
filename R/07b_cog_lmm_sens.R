## ============================================================
## 07b_cog_lmm_sens.R —— 认知 LMM 敏感性：排除 2018 前死亡者（PI 意见四）
## 目的：检验"衰弱组下降更慢"是否由死亡竞争下的幸存者选择驱动
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

d <- core %>% select(ID, age, gender, edu_bl, bl_time, died, dtime_yr) %>%
  left_join(ex %>% select(ID, frail3), by = "ID") %>%
  left_join(long, by = "ID") %>%
  mutate(frail3 = factor(frail3, levels = c("robust", "prefrail", "frail"))) %>%
  filter(!is.na(frail3))

bl <- d %>% filter(wave == 1) %>% select(ID, cog_bl = cog, edu_bl)
thr <- bl %>% group_by(edu_bl) %>% summarise(threshold = mean(cog_bl) - sd(cog_bl), .groups = "drop")
bl <- bl %>% left_join(thr, by = "edu_bl") %>% mutate(impaired = cog_bl <= threshold)
risk_ids <- bl$ID[!bl$impaired & !is.na(bl$impaired)]
sd0 <- sd(bl$cog_bl[bl$ID %in% risk_ids])

d0 <- d %>% filter(ID %in% risk_ids) %>%
  left_join(bl %>% select(ID, cog_bl), by = "ID") %>%
  mutate(cog_z = cog / sd0, time_yr = t_visit - bl_time) %>%
  filter(!is.na(time_yr), time_yr > -0.5, time_yr < 10)

run_lmm <- function(dd, tag) {
  m <- lmer(cog_z ~ frail3 * time_yr + age + gender + factor(edu_bl) + (1 + time_yr | ID), data = dd)
  cf <- summary(m)$coefficients
  b_r <- cf["time_yr", "Estimate"]; b_f <- b_r + cf["frail3frail:time_yr", "Estimate"]
  se <- sqrt(vcov(m)["time_yr", "time_yr"] + vcov(m)["frail3frail:time_yr", "frail3frail:time_yr"] +
             2 * vcov(m)["time_yr", "frail3frail:time_yr"])
  diff <- cf["frail3frail:time_yr", "Estimate"]
  cat("\n== ", tag, " ==\n", sep = "")
  cat("人:", length(unique(dd$ID)), " 观测:", nrow(dd), "\n")
  cat("斜率: 健壮", round(b_r, 4), " 衰弱", round(b_f, 4),
      " | 衰弱-健壮 =", round(diff, 4), "(95%CI", round(diff - 1.96 * se, 4), "至", round(diff + 1.96 * se, 4), ")\n")
  invisible(m)
}

run_lmm(d0, "主分析（全风险集）")

## 敏感性：排除 2018 前死亡者（dtime<7 且 died）+ 要求至少 2 次观测
alive18 <- core$ID[!(core$died == 1 & core$dtime_yr < 7)]
d1 <- d0[d0$ID %in% alive18, ]
cnt <- d1 %>% group_by(ID) %>% summarise(k = n(), .groups = "drop") %>% filter(k >= 2)
d1 <- d1[d1$ID %in% cnt$ID, ]
run_lmm(d1, "敏感性：排除 2018 前死亡者且≥2 次观测")

## 敏感性2：排除全部死亡者（纯存活轨迹）
alive_all <- core$ID[core$died == 0]
d2 <- d0[d0$ID %in% alive_all, ]
cnt2 <- d2 %>% group_by(ID) %>% summarise(k = n(), .groups = "drop") %>% filter(k >= 2)
d2 <- d2[d2$ID %in% cnt2$ID, ]
run_lmm(d2, "敏感性：仅存活者（排除全部死亡）")
cat("\n07b 完成\n")
