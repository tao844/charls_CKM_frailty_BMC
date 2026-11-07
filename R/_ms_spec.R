suppressMessages(library(dplyr)); suppressMessages(library(survival)); suppressMessages(library(haven)); suppressMessages(library(jsonlite))
proj <- "D:/00 陶万清的工作/07 working stata/01 CKM_frailty_revision"
core <- readRDS(file.path(proj, "output/rebuild/bl_cov_step2.rds"))
ex   <- readRDS(file.path(proj, "output/rebuild/exposures_step3.rds"))
cov10 <- c("gender","edu_bl","rural_bl","marry_bl","smoken_bl","drinkev_bl",
           "sleep_bl","ln_hhcperc","ins_bl","pension_bl")
out <- list()

## ---------- 12,754 描述 + Table 1 ----------
d <- core %>% left_join(ex, by="ID") %>%
  filter(!is.na(FI41), !is.na(CKM_Stage), !is.na(dtime_yr), !is.na(died)) %>%
  mutate(age_exit = age + dtime_yr,
         joint6 = factor(interaction(CKM2, frail3, sep=" | "),
           levels=c("0-2 | robust","0-2 | prefrail","0-2 | frail",
                    "3-4 | robust","3-4 | prefrail","3-4 | frail")))
out$n12754 <- nrow(d); out$deaths12754 <- sum(d$died)
out$age_mean <- round(mean(d$age),1); out$age_sd <- round(sd(d$age),1)
out$male_pct <- round(mean(d$gender, na.rm=TRUE)*100,1)
msd <- function(x) sprintf("%.2f (%.2f)", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
pct1 <- function(x) sprintf("%.1f", mean(x, na.rm=TRUE)*100)
t1 <- list()
for (g in levels(d$joint6)) {
  gg <- d[d$joint6==g, ]
  t1[[g]] <- list(n = nrow(gg),
    age = msd(gg$age), male = msd(gg$gender), edu = msd(gg$edu_bl),
    rural = msd(gg$rural_bl), marry = msd(gg$marry_bl), smoke = msd(gg$smoken_bl),
    drink = msd(gg$drinkev_bl), sleep = msd(gg$sleep_bl),
    ln_hhcperc = msd(gg$ln_hhcperc), ins = msd(gg$ins_bl), pension = msd(gg$pension_bl),
    FI = msd(gg$FI41))
}
out$table1 <- t1
writeLines(toJSON(out, auto_unbox=TRUE, pretty=TRUE), file.path(proj, "output/rebuild/_ms_spec.json"))
cat("part1 ok\n")

## ---------- 死亡模型派生量 ----------
dd <- d %>% filter(complete.cases(.[, c("age", cov10)]))
m2 <- coxph(Surv(age, age_exit, died) ~ joint6 + gender + edu_bl + rural_bl +
            marry_bl + smoken_bl + drinkev_bl + sleep_bl + ln_hhcperc +
            ins_bl + pension_bl, cluster=householdID, data=dd, x=TRUE)
b <- coef(m2); V <- vcov(m2)
jb <- names(b)[grep("joint6", names(b))]
con <- function(w) { est <- sum(w*b[jb]); se <- sqrt(as.numeric(t(w) %*% V[jb,jb] %*% w));
  list(HR=round(exp(est),3), lo=round(exp(est-1.96*se),3), hi=round(exp(est+1.96*se),3), P=signif(2*pnorm(-abs(est/se)),3)) }
out$ws_34_frail  <- con(c(0,0,-1,0,1))
out$ws_34_prefrail <- con(c(0,0,-1,1,0))
out$tr_02 <- con(c(1,2,0,0,0)); out$tr_34 <- con(c(0,0,-3,1,2))
i34f <- jb[5]; i34r <- jb[3]; i02f <- jb[2]
th <- unname(b[c(i34f,i34r,i02f)]); grad <- c(exp(th[1]), -exp(th[2]), -exp(th[3]))
V3 <- V[c(i34f,i34r,i02f), c(i34f,i34r,i02f)]
reri <- exp(th[1]) - exp(th[2]) - exp(th[3]) + 1
se <- sqrt(as.numeric(t(grad) %*% V3 %*% grad))
out$reri <- list(est=round(reri,3), lo=round(reri-1.96*se,3), hi=round(reri+1.96*se,3))
cat("part2 ok\n")

## ---------- 认知主模型（教育分层界值，与 FG 同人群 9,140）----------
hh <- read_dta(file.path(proj, "CHARLS_official/H_CHARLS_D_Data (1)/H_CHARLS_D_Data.dta"))
hh$ID <- sprintf("%012.0f", as.numeric(as.character(hh$ID)))
mk <- function(i) data.frame(ID=hh$ID, cog=(as.numeric(hh[[paste0("r",i,"imrc")]])+
  as.numeric(hh[[paste0("r",i,"dlrc")]]))/2 + as.numeric(hh[[paste0("r",i,"ser7")]])+
  as.numeric(hh[[paste0("r",i,"orient")]]) + as.numeric(hh[[paste0("r",i,"draw")]]),
  t=as.numeric(hh[[paste0("r",i,"iwy")]]) + (as.numeric(hh[[paste0("r",i,"iwm")]])-0.5)/12, wave=i)
cl <- rbind(mk(2),mk(3),mk(4)) %>% filter(!is.na(cog), !is.na(t))
bl <- dd %>% select(ID, age, edu_bl, bl_time, died, dtime_yr, householdID, joint6, all_of(cov10)) %>%
  left_join(mk(1) %>% select(ID, cog_h=cog), by="ID")
thr <- bl %>% group_by(edu_bl) %>% summarise(thr=mean(cog_h, na.rm=TRUE)-sd(cog_h, na.rm=TRUE), .groups="drop")
bl <- bl %>% left_join(thr, by="edu_bl") %>% filter(!is.na(cog_h), !is.na(thr))
risk <- bl %>% filter(cog_h > thr)
out$cog_risk551 <- nrow(risk)
first_ev <- cl %>% right_join(risk %>% select(ID, thr), by="ID") %>%
  filter(!is.na(cog), cog <= thr) %>% group_by(ID) %>% summarise(t_ev=min(t), .groups="drop")
cg <- risk %>% left_join(first_ev, by="ID") %>%
  mutate(t_ev=ifelse(is.na(t_ev), Inf, t_ev - bl_time),
         ftime=pmin(t_ev, dtime_yr),
         fstatus=ifelse(t_ev < dtime_yr, 1, ifelse(died==1, 2, 0)),
         ftime=ifelse(ftime==Inf, dtime_yr, ftime))
out$cog_primary_n <- nrow(cg); out$cog_primary_events <- sum(cg$fstatus==1)
out$cog_competing <- sum(cg$fstatus==2)
mcc <- coxph(Surv(ftime, fstatus==1) ~ joint6 + age + gender + edu_bl + rural_bl +
             marry_bl + smoken_bl + drinkev_bl + sleep_bl + ln_hhcperc +
             ins_bl + pension_bl, cluster=householdID, data=cg)
cf <- summary(mcc)$coefficients; ci <- summary(mcc)$conf.int
cog_tab <- list()
for (k in rownames(cf)[grep("joint6", rownames(cf))]) {
  g <- sub("joint6", "", k)
  cog_tab[[g]] <- list(HR=round(cf[k,"exp(coef)"],3), lo=round(ci[k,"lower .95"],3),
                       hi=round(ci[k,"upper .95"],3), P=signif(cf[k,ncol(cf)],3),
                       events=as.integer(sum(cg$joint6==g & cg$fstatus==1)),
                       n=as.integer(sum(cg$joint6==g)))
}
out$cog_primary <- cog_tab
writeLines(toJSON(out, auto_unbox=TRUE, pretty=TRUE, na="string"), file.path(proj, "output/rebuild/_ms_spec.json"))
cat("part3 ok; 认知主模型 n=", out$cog_primary_n, " events=", out$cog_primary_events, "\n")
print(cog_tab)
