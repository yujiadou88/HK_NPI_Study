#
# R syntax to reproduce information for Figure 2 from:
#
# Ng S, Cowling BJ, Fang VJ, Chan KH, Ip DKM, et al.
# Effects of oseltamivir treatment on duration of clinical illness
# and viral shedding and household transmission of influenza virus.
# CID, 2010; 50:.
#
# Last updated by Vicky Fang and Sophie Ng
# February 5, 2010

source("../Oseltamivir_scripts/dataframe.r")  
require(chron)
library(gee)

tab4 <- matrix(rep(NA,364),ncol=13)
tab4[c(4,8,9,11,13,18,19,21,23,27),c(2,6,11)] <- 1


# Imputing age by relation, occupation and schooling
for (i in 1:nrow(trans.c)){
    if (trans.c$hhID[i]==7001 & trans.c$member[i]==2) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==7090 & trans.c$member[i]==3) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==7093 & trans.c$member[i]==2) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==7148 & trans.c$member[i]==2) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8148 & trans.c$member[i]==1) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8148 & trans.c$member[i]==2) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8148 & trans.c$member[i]==3) {trans.c$age[i] <-6}
    if (trans.c$hhID[i]==8273 & trans.c$member[i]==1) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8273 & trans.c$member[i]==2) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8283 & trans.c$member[i]==3) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8301 & trans.c$member[i]==1) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8317 & trans.c$member[i]==1) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8317 & trans.c$member[i]==2) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8317 & trans.c$member[i]==3) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8317 & trans.c$member[i]==4) {trans.c$age[i] <-5}
    if (trans.c$hhID[i]==8317 & trans.c$member[i]==5) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8398 & trans.c$member[i]==1) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8398 & trans.c$member[i]==2) {trans.c$age[i] <-18}
}
index.hh2 <- index.hh[c("hhID","agegp","male","vaccine","flu.type")]
names(index.hh2) <- c("hhID","indexagegp","indexmale","indexvaccine","flu.type")
trans.data <- merge(trans.c,index.hh2,all.x=TRUE)
trans.data$agegp <- cut(trans.data$age,c(-0.1,5,12,17,100))
trans.data$agegp <- factor(trans.data$agegp,labels=c("2youngerchildren0-5yo","3olderchildren6-12yo","4adolescents13-17yo","1adults18+yo"))
trans.data$vaccine[is.na(trans.data$vaccine)] <- 0

# oseltamivir start day
tamiflu <- merge(av[av$member==0&av$av=="tamiflu",],hchar[c("hhID","clinic_date","clinic_day")],by="hhID",all.x=TRUE)
tamiflu$onset.tami <- dates(as.character(tamiflu$dayfrom),format="d/m/y")-dates(as.character(tamiflu$clinic_date),format="d/m/y")+tamiflu$clinic_day
tamiflu$onset.tami[tamiflu$onset.tami>=2] <- 2
tamiflu$onset.tami <- tamiflu$onset.tami+1
trans.data <- merge(trans.data,tamiflu[c("hhID","onset.tami")],all.x=T)
trans.data$onset.tami[is.na(trans.data$onset.tami)] <- 0

# add intervention and study year
trans.data$year <- 1*(trans.data$hhID>8000)
trans.data <- merge(trans.data,hchar[1:2],all.x=T)

# Regression

# lab-confirmed transmission
lab.gee <- gee(labsedcase~factor(onset.tami)+factor(as.character(indexagegp))+indexmale+indexvaccine+factor(flu.type)
                     +factor(as.character(agegp))+male+vaccine+factor(intervention)+year,
                     id=factor(hhID),data=trans.data,corstr = "exchangeable", family="binomial")
results <- data.frame(beta=lab.gee$coef, se=sqrt(diag(lab.gee[[20]])),  row.names=names(lab.gee$coef))[-1,]
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results$p.value <- 2*pnorm(-1*abs(results$beta/results$se))
tab4[c(1:3,5:7,10,12,14:17,20,22,24:26,28),2:5] <- as.matrix(round(results, 2)[,-(1:2)])
tab4[,1] <- c(table(trans.data$onset.tami[trans.data$member==1])[c(2:4,1)], table(trans.data$indexagegp[trans.data$member==1]),
              table(trans.data$indexmale[trans.data$member==1]), table(trans.data$indexvaccine[trans.data$member==1]),
              table(trans.data$flu.type[trans.data$member==1]), table(trans.data$agegp), table(trans.data$male),
              table(trans.data$vaccine), table(trans.data$intervention[trans.data$member==1]),table(trans.data$year[trans.data$member==1]))

# Clinical secondary infection
clin.gee <- gee(clinicsedcase~factor(onset.tami)+factor(as.character(indexagegp))+indexmale+indexvaccine+factor(flu.type)
                     +factor(as.character(agegp))+male+vaccine+factor(intervention)+year,
                     id=factor(hhID),data=trans.data,corstr = "exchangeable", family="binomial")
results <- data.frame(beta=clin.gee$coef, se=sqrt(diag(clin.gee[[20]])), row.names=names(clin.gee$coef))[-1,]
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results$p.value <- 2*pnorm(-1*abs(results$beta/results$se))
tab4[c(1:3,5:7,10,12,14:17,20,22,24:26,28),6:9] <- as.matrix(round(results, 2)[,-(1:2)])

# Clinical secondary infection with laboratory confirmation
trans.data$labc2nd <- 1*(trans.data$labsedcase==1&trans.data$clinicsedcase==1)
cal.gee <- gee(labc2nd~factor(onset.tami)+factor(as.character(indexagegp))+indexmale+indexvaccine+factor(flu.type)
                     +factor(as.character(agegp))+male+vaccine+factor(intervention)+year,
                     id=factor(hhID),data=trans.data,corstr = "exchangeable", family="binomial")
results <- data.frame(beta=cal.gee$coef, se=sqrt(diag(cal.gee[[20]])), row.names=names(cal.gee$coef))[-1,]
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results$p.value <- 2*pnorm(-1*abs(results$beta/results$se))
tab4[c(1:3,5:7,10,12,14:17,20,22,24:26,28),10:13] <- as.matrix(round(results, 2)[,-(1:2)])

colnames(tab4) <- c("n","Lab-OR","95%CI_low","95%CI_up","p-value","Clinic-OR","95%CI_low","95%CI_up","p-value","Lab&Clinic-OR","95%CI_low","95%CI_up","p-value")
rownames(tab4) <- c("Oseltamivir<1d","Oseltamivir1-2d","Oseltamivir>2d","No antiviral","Age<5","Age6-12","Age13-17","Age18+",
                    "Female","Male","Not vaccine","Vaccine","FluA","FluB","Age<5","Age6-12","Age13-17","Age18+","Female","Male",
                    "Not vaccine","Vaccine","None","Mask","Hand hygien","Mask+HH","Year2007","Year2008")

tab4

#
# End of script
#





