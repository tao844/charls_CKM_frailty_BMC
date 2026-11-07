## ============================================================
## 06_mi_m20.R —— 多重插补 m=20（Major 11）
## 【口径 2026-11-04 PI 定死】
##   mice m=20，seed = 20261104；
##   插补变量集 = 统一 13 协变量（原始未缺失指示化版本）+ died + cog_event
##                + 各自 Nelson-Aalen + householdID（meth 置 ""，仅用于聚类）；
##   非插补变量 meth 一律 ""；方法默认 pmm；
##   合并：死亡 = 全队列 Cox；认知 = 风险集 Cox（教育分层界值主定义）
## 对照输出：pool 后 HR vs 完整案例主模型（重点 exercise_bl 方向与量级）
##          + 每模型 rvi/fmi 摘要
## ============================================================
source("00_setup.R")
library(mice)
deaths <- read.csv("output/deaths.csv", fileEncoding = "UTF-8")

## ---- 结局帧构建 ----
fr <- bl  # 原始协变量（未缺失指示化）
fr <- merge(fr, deaths, by = "ID", all.x = TRUE)
stopifnot(!any(is.na(fr$died)))

## 认知事件（教育分层界值主定义，同 05）
fu <- aggregate(time ~ ID, data = lng, FUN = max)
fr$lastfu <- fu$time[match(fr$ID, fu$ID)]
lci <- lng[order(lng$ID, lng$time), ]
lci <- lci[!is.na(lci$total_cognition), ]
w1c <- lng[lng$wave == 1 & !is.na(lng$total_cognition) & !is.na(lng$edu_bl), ]
cut_tab <- do.call(rbind, lapply(sort(unique(w1c$edu_bl)), function(e) {
  x <- w1c$total_cognition[w1c$edu_bl == e]
  data.frame(edu_bl = e, cutoff = mean(x) - sd(x))
}))
lci$cutval <- cut_tab$cutoff[match(lci$edu_bl, cut_tab$edu_bl)]
lci$ci_e <- as.integer(!is.na(lci$cutval) & lci$total_cognition < lci$cutval)
first_e <- lci[lci$ci_e == 1, ]; first_e <- first_e[!duplicated(first_e$ID), c("ID", "time")]
names(first_e) <- c("ID", "ctime_e")
fr <- merge(fr, first_e, by = "ID", all.x = TRUE)
fr$cog_event <- as.integer(!is.na(fr$ctime_e))
fr$ci_tm <- ifelse(fr$cog_event == 1, fr$ctime_e,
                   ifelse(fr$died == 1, pmin(fr$lastfu, fr$dtime_yr, na.rm = TRUE),
                          fr$lastfu))
## 非风险集者认知结局不定义
fr$cog_event[fr$cog_impair_bl != 0] <- NA
fr$ci_tm[fr$cog_impair_bl != 0] <- NA

## Nelson-Aalen 累积 hazard
bh_d <- basehaz(coxph(Surv(dtime_yr, died) ~ 1, data = fr))
fr$nelson <- bh_d$hazard[findInterval(fr$dtime_yr, bh_d$time)]
tmp_c <- fr[fr$cog_impair_bl == 0, ]
bh_c <- basehaz(coxph(Surv(ci_tm, cog_event) ~ 1, data = tmp_c))
fr$nelson_ci <- NA
fr$nelson_ci[fr$cog_impair_bl == 0] <-
  bh_c$hazard[findInterval(tmp_c$ci_tm, bh_c$time)]

## ---- 插补框架 ----
## CKM34/frail2/frailord 为派生确定性列（joint6/CKM_Stage_bl 的函数），
## 留在框内会引发近共线（mice ridge 警告 1,429 条），剔除；pool 模型亦不依赖它们
fr$CKM34  <- as.integer(fr$CKM_Stage_bl >= 3)
imp_vars <- c("joint6", "CKM_Stage_bl",
              "cog_impair_bl", "householdID",
              "age_bl", "gender", "edu_bl", "rural_bl", "marry_bl", "nation_bl",
              "smoken_bl", "drinkev_bl", "exercise_bl", "sleep_bl",
              "hhcperc_bl", "ins_bl", "pension_bl",
              "died", "dtime_yr", "nelson",
              "cog_event", "ci_tm", "nelson_ci")
imp0 <- fr[, imp_vars]
meth <- make.method(imp0)
no_imp <- c("joint6", "CKM_Stage_bl",
            "cog_impair_bl", "householdID",
            "age_bl", "gender", "rural_bl", "marry_bl",
            "died", "dtime_yr", "nelson", "ci_tm")
meth[no_imp] <- ""
## Nelson-Aalen 与对应时间/事件变量函数相依或带结构性缺失块、householdID 仅为
## 聚类用：保留在插补变量集内（作为插补目标），但不作为其他变量的预测因子
predM <- make.predictorMatrix(imp0)
predM[, c("householdID", "nelson", "nelson_ci", "ci_tm", "cog_event")] <- 0
cat("插补方法表（非 \"\" 者）：\n")
print(meth[meth != ""])

## ---- mice m=20（已存在则复用，避免重复耗时）----
if (file.exists("output/mi_m20.rds")) {
  imp <- readRDS("output/mi_m20.rds")
  cat("复用已存在的 output/mi_m20.rds\n")
} else {
  imp <- mice(imp0, m = 20, method = meth, predictorMatrix = predM,
              seed = 20261104, printFlag = FALSE)
  saveRDS(imp, "output/mi_m20.rds")
  cat("mice 完成并保存 output/mi_m20.rds\n")
}

## ---- 每次插补内重构模型协变量并拟合（with.mids 中变量直接在环境中可用）----
f_rhs <- paste("joint6 + age_bl + gender + edu_bl + rural_bl + marry_bl +",
               "nation_f + smoken_bl + drinkev_bl + exercise_f + sleep_bl +",
               "ln_hhcperc + ins_bl + pension_f")
fit_d <- with(imp, {
  ln_hhcperc <- log(hhcperc_bl)
  ln_hhcperc[!is.finite(ln_hhcperc)] <- NA  # 1 人消费为 0，coxph 自动剔除该行者
  exercise_f <- factor(ifelse(exercise_bl == 1, "yes", "no"), levels = c("no", "yes"))
  pension_f  <- factor(ifelse(pension_bl == 1, "yes", "no"), levels = c("no", "yes"))
  nation_f   <- factor(ifelse(nation_bl == 1, "minority", "han"), levels = c("han", "minority"))
  coxph(as.formula(paste("Surv(dtime_yr, died) ~", f_rhs)),
        cluster = householdID)
})
pool_d <- pool(fit_d)
sd_d <- summary(pool_d)
sd_d$HR <- round(exp(sd_d$estimate), 3)
sd_d$lo <- round(exp(sd_d$estimate - 1.96 * sd_d$std.error), 3)
sd_d$hi <- round(exp(sd_d$estimate + 1.96 * sd_d$std.error), 3)
write.csv(sd_d, "output/tables/Death_MI20_pooled_full.csv", row.names = FALSE, fileEncoding = "UTF-8")

fit_c <- with(imp, {
  ln_hhcperc <- log(hhcperc_bl)
  ln_hhcperc[!is.finite(ln_hhcperc)] <- NA
  exercise_f <- factor(ifelse(exercise_bl == 1, "yes", "no"), levels = c("no", "yes"))
  pension_f  <- factor(ifelse(pension_bl == 1, "yes", "no"), levels = c("no", "yes"))
  nation_f   <- factor(ifelse(nation_bl == 1, "minority", "han"), levels = c("han", "minority"))
  keep <- cog_impair_bl == 0
  coxph(as.formula(paste("Surv(ci_tm, cog_event) ~", f_rhs)),
        subset = keep, cluster = householdID)
})
pool_c <- pool(fit_c)
sd_c <- summary(pool_c)
sd_c$HR <- round(exp(sd_c$estimate), 3)
sd_c$lo <- round(exp(sd_c$estimate - 1.96 * sd_c$std.error), 3)
sd_c$hi <- round(exp(sd_c$estimate + 1.96 * sd_c$std.error), 3)
write.csv(sd_c, "output/tables/CI_MI20_pooled.csv", row.names = FALSE, fileEncoding = "UTF-8")

## ---- rvi/fmi 摘要 ----
## 注意：mice pool 对象的行名可能不可靠，统一用 pooled$term 列匹配
fmi_tab <- do.call(rbind, lapply(
  list(death = pool_d, cognition = pool_c), function(p) {
    pv <- p$pooled
    if (!"riv" %in% names(pv)) {
      cat("!! pooled 列名异常："); print(names(pv))
    }
    term <- if ("term" %in% names(pv)) as.character(pv$term) else rownames(pv)
    g <- grepl("joint6|exercise_f|pension_f|nation_f", term)
    if (sum(g) == 0) {
      cat("!! 模式未命中，实际 term 前 20：\n"); print(head(term, 20))
    }
    data.frame(model = "",
               riv_mean = round(mean(pv$riv[g]), 3), riv_max = round(max(pv$riv[g]), 3),
               fmi_mean = round(mean(pv$fmi[g]), 3), fmi_max = round(max(pv$fmi[g]), 3))
  }))
fmi_tab$model <- c("death_full", "cognition_riskset")
print(fmi_tab)
write.csv(fmi_tab, "output/tables/MI20_rvi_fmi_summary.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cat("06 完成（pool 与对照由 06b_mi_pooled_compare.R 输出）\n")
