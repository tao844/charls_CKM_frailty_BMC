## ============================================================
## 01_alt_cohort_flow.R —— 建队流程复核、纳入 vs 排除比较、随访/失访
## 回应：Major 6（37% 排除过大）；Minor（随访中位数、失访人数、风险集基线表）
## 在代码包 01_cohort_flow 基础上，按交接简报 §3.6 增补：
##   重建样本 11,982 与实际分析样本 10,556 的差距分解
##   （memrye==1 隐藏排除 + 未解释排除及其认知分布）
## 数据口径修正（2026-11-04 核查发现）：
##   ful 中 age/total_cognition/协变量均为字符型（"NA" 字符串），需数值化；
##   ful$ID 为字符型（带前导 0），bl$ID 为数值型，匹配前需统一。
## ============================================================
source("00_setup.R")

w1 <- ful[ful$wave == 1, ]
cat("2011 波全量：", nrow(w1), "人\n")
w1$age_num <- as.numeric(w1$age)
w1$cog_num <- as.numeric(w1$total_cognition)

## ---- 建队各步排除（与 Figure 1 / 简报 §3.6 核对）----
## 口径 2026-11-05：本表统一记录"每步排除人数"（差值），不再与起点人数混排；
## 排除链五步：648(<45岁) → 470(卒中) → 4,608(认知缺失/代理) → 143(memrye) → 1,283(其他)
excl_age    <- sum(is.na(w1$age_num) | w1$age_num < 45)
age_ok      <- !is.na(w1$age_num) & w1$age_num >= 45
excl_stroke <- sum(age_ok & w1$stroke == 1, na.rm = TRUE)
stroke_ok   <- is.na(w1$stroke) | w1$stroke != 1
excl_cog    <- sum(age_ok & stroke_ok & is.na(w1$cog_num))
n_rebuild   <- sum(age_ok & stroke_ok & !is.na(w1$cog_num))
## 认知缺失的二成分分解：总分 = 记忆 + 执行，按两分量缺失模式拆分（4,608 = 2,160 + 2,448）
cg <- w1[age_ok & stroke_ok & is.na(w1$cog_num), ]
mem_miss <- is.na(as.numeric(cg$memeory))
exe_miss <- is.na(as.numeric(cg$executive))
excl_cog_dual   <- sum(mem_miss & exe_miss)
excl_cog_single <- sum(xor(mem_miss, exe_miss))
stopifnot(excl_cog_dual + excl_cog_single == excl_cog)

w1$id_num <- suppressWarnings(as.numeric(w1$ID))
elig <- w1[age_ok & stroke_ok & !is.na(w1$cog_num), ]
excl_memrye <- sum(elig$memrye == 1, na.rm = TRUE)
excl_other  <- nrow(elig) - excl_memrye - nrow(bl)
stopifnot(n_rebuild - excl_memrye - excl_other == nrow(bl))

flow <- data.frame(
  step = c("2011波全量(起点)", "年龄<45岁或年龄缺失", "自报卒中",
           "认知缺失/代理访谈", "  其中:双成分缺失(记忆+执行均缺)",
           "  其中:单成分缺失(仅缺其一)", "自报记忆相关疾病(memrye==1)",
           "其他不满足最低派生要求", "最终分析样本(终点)"),
  excluded_n = c(NA, excl_age, excl_stroke, excl_cog, NA, NA,
                 excl_memrye, excl_other, NA),
  note = c("n=17,708", "", "", "4,608=2,160+2,448",
           paste0("n=", excl_cog_dual), paste0("n=", excl_cog_single), "", "",
           paste0("n=", nrow(bl))),
  remaining_n = c(nrow(w1), nrow(w1) - excl_age,
                  nrow(w1) - excl_age - excl_stroke, n_rebuild, NA, NA,
                  n_rebuild - excl_memrye, nrow(bl), nrow(bl)))
cat("逐步排除:", excl_age, excl_stroke, excl_cog,
    "(", excl_cog_dual, "+", excl_cog_single, ")", excl_memrye, excl_other, "\n")
print(flow)
write.csv(flow, "output/tables/Flow_steps_check.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("\n重建样本：", nrow(elig), "人；实际分析样本：", nrow(bl),
    "人；差距：", nrow(elig) - nrow(bl), "人\n")

cat("\n-- memrye==1 排除统计 --\n")
cat("全波1 memrye==1 例数：", sum(w1$memrye == 1, na.rm = TRUE), "\n")
cat("重建样本内 memrye==1 例数：", sum(elig$memrye == 1, na.rm = TRUE), "\n")
cat("进入分析样本者：", sum(elig$memrye == 1 & elig$id_num %in% bl$ID, na.rm = TRUE), "\n")

unexp <- elig[!elig$id_num %in% bl$ID, ]
cat("\n-- 差距分解 --\n")
cat("差距合计：", nrow(unexp), "人\n")
cat("其中 memrye==1：", sum(unexp$memrye == 1, na.rm = TRUE), "\n")
cat("未解释排除人数：", nrow(unexp) - sum(unexp$memrye == 1, na.rm = TRUE), "\n")
cog_un <- as.numeric(unexp$cog_num)
cat("未解释排除者认知分摘要（全差距 1,426 人，含 memrye==1）：\n")
qs_un <- c(mean = mean(cog_un, na.rm = TRUE),
           sd   = sd(cog_un, na.rm = TRUE),
           quantile(cog_un, c(.1, .25, .5, .75, .9), na.rm = TRUE))
print(round(qs_un, 2))
## 口径修正 2026-11-04：落盘表专指"未解释排除"（剔除 memrye==1 后的 1,283 人）
unexp_only <- unexp[is.na(unexp$memrye) | unexp$memrye != 1, ]
cog_uo <- as.numeric(unexp_only$cog_num)
cog_in <- bl$total_cognition
cat("分析样本认知分摘要（对照）：\n")
qs_in <- c(mean = mean(cog_in, na.rm = TRUE),
           sd   = sd(cog_in, na.rm = TRUE),
           quantile(cog_in, c(.1, .25, .5, .75, .9), na.rm = TRUE))
print(round(qs_in, 2))
## 落盘：未解释排除者（1,283，剔除 memrye==1）vs 分析样本的认知分布对照（Major 6 支撑表）
prof <- data.frame(group = c("未解释排除者(剔除memrye==1)", "分析样本对照"),
  n   = c(sum(!is.na(cog_uo)), sum(!is.na(cog_in))),
  mean = round(c(mean(cog_uo, na.rm = TRUE), mean(cog_in, na.rm = TRUE)), 2),
  sd   = round(c(sd(cog_uo, na.rm = TRUE), sd(cog_in, na.rm = TRUE)), 2),
  p10  = round(c(quantile(cog_uo, .1, na.rm = TRUE), quantile(cog_in, .1, na.rm = TRUE)), 2),
  p25  = round(c(quantile(cog_uo, .25, na.rm = TRUE), quantile(cog_in, .25, na.rm = TRUE)), 2),
  p50  = round(c(quantile(cog_uo, .50, na.rm = TRUE), quantile(cog_in, .50, na.rm = TRUE)), 2),
  p75  = round(c(quantile(cog_uo, .75, na.rm = TRUE), quantile(cog_in, .75, na.rm = TRUE)), 2),
  p90  = round(c(quantile(cog_uo, .90, na.rm = TRUE), quantile(cog_in, .90, na.rm = TRUE)), 2))
write.csv(prof, "output/tables/Unexplained_exclusions_profile.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
gap_tab <- data.frame(
  item = c("重建样本","实际分析样本","差距合计","其中 memrye==1","未解释排除"),
  n = c(nrow(elig), nrow(bl), nrow(unexp),
        sum(unexp$memrye == 1, na.rm = TRUE),
        nrow(unexp) - sum(unexp$memrye == 1, na.rm = TRUE)))
print(gap_tab)
write.csv(gap_tab, "output/tables/Gap_11982_vs_10556.csv", row.names = FALSE, fileEncoding = "UTF-8")

## ---- 纳入 vs 排除比较（Major 6）----
w1$included <- w1$id_num %in% bl$ID
cmp <- w1
vars <- c("age","gender","rural","smoken","drinkev","edu")
numv <- function(x) suppressWarnings(as.numeric(x))
tab_inc <- data.frame(
  variable = vars,
  included_n = sapply(vars, function(v) sum(!is.na(numv(cmp[[v]][cmp$included])))),
  included_mean = sapply(vars, function(v) round(mean(numv(cmp[[v]][cmp$included]), na.rm = TRUE), 2)),
  excluded_n = sapply(vars, function(v) sum(!is.na(numv(cmp[[v]][!cmp$included])))),
  excluded_mean = sapply(vars, function(v) round(mean(numv(cmp[[v]][!cmp$included]), na.rm = TRUE), 2)))
print(tab_inc)
write.csv(tab_inc, "output/tables/Included_vs_excluded.csv", row.names = FALSE, fileEncoding = "UTF-8")

## ---- 随访时间（中位数 IQR）与各波在访/失访 ----
fu <- aggregate(time ~ ID, data = lng, FUN = max)   # 每人最后随访年
rs$fu_time <- fu$time[match(rs$ID, fu$ID)]
cat("\n随访时间中位数(IQR)：",
    round(median(rs$fu_time, na.rm = TRUE), 1), " (",
    round(quantile(rs$fu_time, .25, na.rm = TRUE), 1), ",",
    round(quantile(rs$fu_time, .75, na.rm = TRUE), 1), ")\n")

## 各波在访人数 × 联合组（失访表，供 Major 4 IPC 与 Minor 回复）
lng$inwave <- !is.na(lng$total_cognition)
att <- with(lng[lng$ID %in% rs$ID, ],
            table(joint, wave, inwave))
att2 <- as.data.frame(att)
write.csv(att2, "output/tables/Attrition_by_group_wave.csv", row.names = FALSE, fileEncoding = "UTF-8")
print(att)

## ---- 风险集基线特征表（Minor：Table 1 主体改为风险集 n=8,916）----
summ <- function(x) {
  if (is.numeric(x)) paste0(round(mean(x, na.rm = TRUE), 2), " (",
                            round(sd(x, na.rm = TRUE), 2), ")")
  else paste0(table(x), collapse = "/")
}
chars <- c("age_bl","gender","edu_bl","rural_bl","marry_bl","nation_bl",
           "smoken_bl","drinkev_bl","cesd10_bl","ln_hhcperc","FI_bl","total_cognition")
t1 <- do.call(rbind, lapply(split(rs, rs$joint6), function(g)
  data.frame(group = as.character(g$joint6[1]),
             t(sapply(chars, function(v) summ(g[[v]]))))))
write.csv(t1, "output/tables/Table1_riskset.csv", row.names = FALSE, fileEncoding = "UTF-8")
print(t1)
cat("01_alt 完成\n")
