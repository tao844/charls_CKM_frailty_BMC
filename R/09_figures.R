## ============================================================
## 09_figures.R —— 三张图重绘（新队列）
## Figure 1 队列流程（Rebuild_flow_official.csv，17,705 起点）
## Figure 2 KM 曲线 + at-risk 表（真访谈月口径，六格新 n）
## Figure 3 认知森林图（FG SHR 五格 + LMM 斜率差）
## ============================================================
suppressMessages(library(dplyr)); suppressMessages(library(survival)); suppressMessages(library(ggplot2))
proj <- "D:/00 陶万清的工作/07 working stata/01 CKM_frailty_revision"
figdir <- file.path(proj, "figures"); dir.create(figdir, showWarnings = FALSE)
core <- readRDS(file.path(proj, "output/rebuild/bl_cov_step2.rds"))
ex   <- readRDS(file.path(proj, "output/rebuild/exposures_step3.rds"))

## ---- Figure 1 ----
fl <- read.csv(file.path(proj, "output/rebuild/Rebuild_flow_official.csv"), fileEncoding = "UTF-8")
main_steps <- fl[c(1:4, 9, 10), ]
png(file.path(figdir, "Figure1_YMJ_cohort_flow.png"), width = 2200, height = 1600, res = 200)
par(mar = c(1, 1, 3, 1))
n_steps <- nrow(main_steps)
plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
ys <- seq(0.92, 0.08, length.out = n_steps)
boxw <- 0.52; boxh <- 0.09
labels <- c("2011 wave full sample", "Excluded: age <45 or missing", "Excluded: self-reported stroke",
            "Excluded: missing cognition / proxy interview", "Excluded: self-reported memory-related disease (memrye = 1)",
            "Rebuilt analytic cohort")
remain <- main_steps$remaining_new
excl <- main_steps$new_n
for (i in seq_len(n_steps)) {
  if (i == 1) {
    rect(0.5 - boxw/2, ys[i] - boxh/2, 0.5 + boxw/2, ys[i] + boxh/2, col = "#E8EEF7", border = "#2C5F8A", lwd = 2)
    text(0.5, ys[i], sprintf("%s\nn = %s", labels[i], format(remain[i], big.mark = ",")), cex = 1.05, font = 2)
  } else {
    rect(0.5 - boxw/2, ys[i] - boxh/2, 0.5 + boxw/2, ys[i] + boxh/2,
         col = ifelse(i == n_steps, "#DFF0D8", "#F7E8E8"), border = ifelse(i == n_steps, "#3C763D", "#A94442"), lwd = 2)
    text(0.5, ys[i], sprintf("%s\nn = %s", labels[i], format(remain[i], big.mark = ",")), cex = 1.0)
    arrows(0.5, ys[i-1] - boxh/2, 0.5, ys[i] + boxh/2, length = 0.06, lwd = 1.4)
    text(0.66, (ys[i-1] + ys[i]) / 2, sprintf("excluded %s", format(excl[i], big.mark = ",")), cex = 0.85, col = "#A94442")
  }
}
mtext("Figure 1. Cohort flow (rebuilt from official CHARLS 2011 raw data)", 3, line = 1, font = 2, cex = 1.2)
text(0.5, 0.01, "Cognition-exclusion subgroup — by proxy status: proxy interview 1,235 / non-proxy item missing 2,237; by number of missing components: dual-component 1,536 / single-component 1,936. Subsequent steps: covariate derivation 13,256 → model sample 12,182 (Table 1 denominator 12,754).", cex = 0.72)
dev.off()
cat("Figure 1 done\n")

## ---- Figure 2 ----
d <- core %>% left_join(ex, by = "ID") %>%
  filter(!is.na(FI41), !is.na(CKM_Stage)) %>%
  mutate(joint6 = factor(interaction(CKM2, frail3, sep = " | "),
           levels = c("0-2 | robust", "0-2 | prefrail", "0-2 | frail",
                      "3-4 | robust", "3-4 | prefrail", "3-4 | frail"))) %>%
  filter(!is.na(dtime_yr), !is.na(died))
d$age_exit <- d$age + d$dtime_yr
km <- survfit(Surv(dtime_yr, died) ~ joint6, data = d)
ymin <- min(km$surv, na.rm = TRUE)
ylo <- max(0, floor((ymin - 0.03) * 100) / 100)
cat("KM 样本:", nrow(d), " 曲线最低生存率:", round(ymin, 3), " ylim 下限:", ylo, "\n")
cols <- c("#2166AC", "#67A9CF", "#D1E5F0", "#B2182B", "#EF8A62", "#FDDBC7")
labs <- c("0-2 robust", "0-2 prefrail", "0-2 frail", "3-4 robust", "3-4 prefrail", "3-4 frail")
png(file.path(figdir, "Figure2_YMJ_KM_mortality.png"), width = 2400, height = 2100, res = 200)
par(mar = c(8, 6, 3, 2), mfrow = c(2, 1), heights = c(3, 1.15))
plot(km, col = cols, lwd = 2.2, lty = c(1,1,1,2,2,2), xlab = "", ylab = "Survival probability", xaxt = "n", cex.lab = 1.1, ylim = c(ylo, 1), xlim = c(-3.6, 9))
axis(1, at = c(0, 2, 4, 7, 9), labels = c(0, 2, 4, 7, 9))
legend("bottomleft", labs, col = cols, lwd = 2.2, lty = c(1,1,1,2,2,2), bty = "n", cex = 0.95, ncol = 2)
mtext("Years since baseline (interview-date follow-up)", 1, line = 2.6, cex = 1.05)
mtext("Figure 2. Kaplan–Meier survival by CKM stage × frailty (rebuilt cohort)", 3, line = 1, font = 2, cex = 1.1, outer = FALSE, adj = 0)
par(mar = c(2, 6, 1, 2))
tt <- c(0, 2, 4, 7, 9)
nr <- t(sapply(levels(d$joint6), function(g) {
  dd <- d[d$joint6 == g, ]
  sapply(tt, function(t) sum(dd$dtime_yr >= t))
}))
plot.new(); plot.window(xlim = c(-3.6, 9), ylim = c(0, 7))
for (j in seq_len(6)) {
  yj <- 6.3 - (j - 1) * 1.05
  text(-0.35, yj, labs[j], adj = 1, cex = 0.92)
  for (k in seq_len(5)) text(tt[k], yj, format(nr[j, k], big.mark = ","), cex = 0.92)
}
for (k in seq_len(5)) text(tt[k], 6.85, paste0("year ", tt[k]), cex = 0.92, font = 2)
text(4.5, -0.75, "No. at risk", cex = 1.0, font = 2)
dev.off()
cat("Figure 2 done\n")

## ---- Figure 3 ----
fg <- data.frame(
  cell = c("0-2 prefrail", "0-2 frail", "3-4 robust", "3-4 prefrail", "3-4 frail"),
  est = c(0.990, 1.048, 1.096, 0.898, 0.920),
  lo  = c(0.878, 0.885, 0.900, 0.769, 0.770),
  hi  = c(1.117, 1.240, 1.333, 1.049, 1.100))
fg$cell <- factor(fg$cell, levels = rev(fg$cell))
p <- ggplot(fg, aes(x = cell, y = est)) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.7, color = "#2166AC") +
  coord_flip(ylim = c(0.7, 1.4)) +
  xlab("") + ylab("Subdistribution hazard ratio (95% CI)") +
  ggtitle("Figure 3. Incident cognitive impairment: Fine–Gray SHR by cell",
          subtitle = "Reference: CKM 0-2 robust; risk set n = 9,140 (1,908 events, 731 competing deaths)") +
  theme_bw(base_size = 13)
ggsave(file.path(figdir, "Figure3_YMJ_cognition_forest.png"), p, width = 8.5, height = 4.5, dpi = 200)
cat("Figure 3 done\n")

## at-risk 首列表数核对输出
print(nr[, 1])
cat("at-risk year0 六格:", paste(nr[, 1], collapse = " / "), "\n")
