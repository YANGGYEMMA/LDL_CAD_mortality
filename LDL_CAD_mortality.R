library(data.table)
library(stringr)
library(readxl)
library(dplyr)
library(tibble)
library(SUMnlmr)
library(survival)
library(survminer)
library(scales)
library(RNOmni)

setwd("D:/MR/Files")

results<-data.frame()#lmr#
results2<-data.frame()#fp#
results3<-data.frame()#pw#

####Ratio method (apob)####
test_data_all<-fread("D:/MR/LDL_CAD/R/Revise3/test_data_all.csv")

test_data_all$apob<-RankNorm(test_data_all$apob)
test_data_all$ldl<-RankNorm(test_data_all$ldl)
test_data_all$tg<-RankNorm(test_data_all$tg)

#analyses#
for (k in c('overall',"men","women","young","old")) {
  if(k=="overall"){test_data<-test_data_all}
  if(k=="men"){test_data<-test_data_all[test_data_all$sex_adj==1,]}
  if(k=="women"){test_data<-test_data_all[test_data_all$sex_adj==0,]}
  if(k=="young"){test_data<-test_data_all[test_data_all$ages<65,]}
  if(k=="old"){test_data<-test_data_all[!test_data_all$ages<65,]}
  
  attach(test_data)
  
  #bx#
  bx = lm(apob~apob_sum
          +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
          +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
          +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)$coef[2] 
  
  #by for chd#
  model = glm(chd~apob_sum
              +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
              +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
              +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20
              , family=binomial)
  
  by = model$coef[2] 
  byse = summary(model)$coef[2,2]
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"apob","chd","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  #by for dead (cox)#
  model = coxph(Surv(outagef,dead)~apob_sum
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  by = summary(model)$coef[1] 
  byse = summary(model)$coef[1,3] 
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"apob","dead","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  #by for cv_dead (cox)#
  model = coxph(Surv(outagef,cv_dead)~apob_sum
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  by = summary(model)$coef[1] 
  byse = summary(model)$coef[1,3] 
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"apob","cv_dead","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  #by for tumo_dead (cox)#
  model = coxph(Surv(outagef,tumo_dead)~apob_sum
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  by = summary(model)$coef[1] 
  byse = summary(model)$coef[1,3] 
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"apob","tumo_dead","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  #by for other_dead (cox)#
  model = coxph(Surv(outagef,other_dead)~apob_sum
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  by = summary(model)$coef[1] 
  byse = summary(model)$coef[1,3] 
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"apob","other_dead","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  detach(test_data)
  rm(test_data,bx,k)}

rm(test_data_all)
gc()
####Ratio method (ldl)####
test_data_all<-fread("D:/MR/LDL_CAD/R/Revise3/test_data_all.csv")

test_data_all$apob<-RankNorm(test_data_all$apob)
test_data_all$ldl<-RankNorm(test_data_all$ldl)
test_data_all$tg<-RankNorm(test_data_all$tg)

#analyses#
for (k in c('overall',"men","women","young","old")) {
  if(k=="overall"){test_data<-test_data_all}
  if(k=="men"){test_data<-test_data_all[test_data_all$sex_adj==1,]}
  if(k=="women"){test_data<-test_data_all[test_data_all$sex_adj==0,]}
  if(k=="young"){test_data<-test_data_all[test_data_all$ages<65,]}
  if(k=="old"){test_data<-test_data_all[!test_data_all$ages<65,]}
  
  attach(test_data)
  
  #bx#
  bx = lm(ldl~ldl_sum
          +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
          +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
          +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)$coef[2] 
  
  #by for chd#
  model = glm(chd~ldl_sum
              +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
              +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
              +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20
              , family=binomial)
  
  by = model$coef[2] 
  byse = summary(model)$coef[2,2]
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"ldl","chd","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  #by for dead (cox)#
  model = coxph(Surv(outagef,dead)~ldl_sum
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  by = summary(model)$coef[1] 
  byse = summary(model)$coef[1,3] 
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"ldl","dead","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  #by for cv_dead (cox)#
  model = coxph(Surv(outagef,cv_dead)~ldl_sum
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  by = summary(model)$coef[1] 
  byse = summary(model)$coef[1,3] 
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"ldl","cv_dead","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  #by for tumo_dead (cox)#
  model = coxph(Surv(outagef,tumo_dead)~ldl_sum
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  by = summary(model)$coef[1] 
  byse = summary(model)$coef[1,3] 
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"ldl","tumo_dead","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  #by for other_dead (cox)#
  model = coxph(Surv(outagef,other_dead)~ldl_sum
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  by = summary(model)$coef[1] 
  byse = summary(model)$coef[1,3] 
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"ldl","other_dead","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  detach(test_data)
  rm(test_data,bx,k)}

rm(test_data_all)
gc()
####Ratio method (tg)####
test_data_all<-fread("D:/MR/LDL_CAD/R/Revise3/test_data_all.csv")

test_data_all$apob<-RankNorm(test_data_all$apob)
test_data_all$ldl<-RankNorm(test_data_all$ldl)
test_data_all$tg<-RankNorm(test_data_all$tg)

#analyses#
for (k in c('overall',"men","women","young","old")) {
  if(k=="overall"){test_data<-test_data_all}
  if(k=="men"){test_data<-test_data_all[test_data_all$sex_adj==1,]}
  if(k=="women"){test_data<-test_data_all[test_data_all$sex_adj==0,]}
  if(k=="young"){test_data<-test_data_all[test_data_all$ages<65,]}
  if(k=="old"){test_data<-test_data_all[!test_data_all$ages<65,]}
  
  attach(test_data)
  
  #bx#
  bx = lm(tg~tg_sum
          +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
          +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
          +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)$coef[2] 
  
  #by for chd#
  model = glm(chd~tg_sum
              +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
              +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
              +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20
              , family=binomial)
  
  by = model$coef[2] 
  byse = summary(model)$coef[2,2]
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"tg","chd","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  #by for dead (cox)#
  model = coxph(Surv(outagef,dead)~tg_sum
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  by = summary(model)$coef[1] 
  byse = summary(model)$coef[1,3] 
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"tg","dead","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  #by for cv_dead (cox)#
  model = coxph(Surv(outagef,cv_dead)~tg_sum
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  by = summary(model)$coef[1] 
  byse = summary(model)$coef[1,3] 
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"tg","cv_dead","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  #by for tumo_dead (cox)#
  model = coxph(Surv(outagef,tumo_dead)~tg_sum
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  by = summary(model)$coef[1] 
  byse = summary(model)$coef[1,3] 
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"tg","tumo_dead","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  #by for other_dead (cox)#
  model = coxph(Surv(outagef,other_dead)~tg_sum
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  by = summary(model)$coef[1] 
  byse = summary(model)$coef[1,3] 
  
  beta.ratio = by/bx
  se.ratio = byse/bx
  z.ratio = abs(beta.ratio/se.ratio)
  p.ratio = 2*pt(z.ratio, df=Inf, lower.tail=FALSE)
  
  dat_tab<-data.frame(k,"tg","other_dead","GRS (unadj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","bx","k")])
  
  detach(test_data)
  rm(test_data,bx,k)}

rm(test_data_all)
gc()
####Two stage least squares (apob_adj)####
test_data_all<-fread("D:/MR/LDL_CAD/R/Revise3/test_data_all.csv")

test_data_all$apob<-RankNorm(test_data_all$apob)
test_data_all$ldl<-RankNorm(test_data_all$ldl)
test_data_all$tg<-RankNorm(test_data_all$tg)

#analyses#
for (k in c('overall',"men","women","young","old")) {
  if(k=="overall"){test_data<-test_data_all}
  if(k=="men"){test_data<-test_data_all[test_data_all$sex_adj==1,]}
  if(k=="women"){test_data<-test_data_all[test_data_all$sex_adj==0,]}
  if(k=="young"){test_data<-test_data_all[test_data_all$ages<65,]}
  if(k=="old"){test_data<-test_data_all[!test_data_all$ages<65,]}
  
  attach(test_data)
  
  #genetically predicted apob and tg using all snps and relevant weights#
  apob_fit = lm(apob~apobtg_apob_sum
                +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)$fitted
  
  tg_fit = lm(tg~apobtg_tg_sum
              +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
              +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
              +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)$fitted
  
  #chd#
  model = glm(chd~apob_fit+tg_fit
              +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
              +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
              +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20
              , family=binomial)
  
  beta.ratio = summary(model)$coef[2,1]
  se.ratio = summary(model)$coef[2,2]
  p.ratio = summary(model)$coef[2,4]
  
  dat_tab<-data.frame(k,"apob","chd","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","apob_fit","tg_fit","k")])
  
  #dead#
  model = coxph(Surv(outagef,dead)~apob_fit+tg_fit
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  beta.ratio = summary(model)$coef[1,1] 
  se.ratio = summary(model)$coef[1,3] 
  p.ratio = summary(model)$coef[1,5]
  
  dat_tab<-data.frame(k,"apob","dead","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","apob_fit","tg_fit","k")])
  
  #cv_dead#
  model = coxph(Surv(outagef,cv_dead)~apob_fit+tg_fit
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  beta.ratio = summary(model)$coef[1,1] 
  se.ratio = summary(model)$coef[1,3] 
  p.ratio = summary(model)$coef[1,5]
  
  dat_tab<-data.frame(k,"apob","cv_dead","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","apob_fit","tg_fit","k")])
  
  #tumo_dead#
  model = coxph(Surv(outagef,tumo_dead)~apob_fit+tg_fit
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  beta.ratio = summary(model)$coef[1,1] 
  se.ratio = summary(model)$coef[1,3] 
  p.ratio = summary(model)$coef[1,5]
  
  dat_tab<-data.frame(k,"apob","tumo_dead","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","apob_fit","tg_fit","k")])
  
  #other_dead#
  model = coxph(Surv(outagef,other_dead)~apob_fit+tg_fit
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  beta.ratio = summary(model)$coef[1,1] 
  se.ratio = summary(model)$coef[1,3] 
  p.ratio = summary(model)$coef[1,5]
  
  dat_tab<-data.frame(k,"apob","other_dead","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","apob_fit","tg_fit","k")])
  
  detach(test_data)
  rm(test_data,apob_fit,tg_fit,k)}

rm(test_data_all)
gc()
####Two stage least squares (ldl_adj)####
test_data_all<-fread("D:/MR/LDL_CAD/R/Revise3/test_data_all.csv")

test_data_all$apob<-RankNorm(test_data_all$apob)
test_data_all$ldl<-RankNorm(test_data_all$ldl)
test_data_all$tg<-RankNorm(test_data_all$tg)

#analyses#
for (k in c('overall',"men","women","young","old")) {
  if(k=="overall"){test_data<-test_data_all}
  if(k=="men"){test_data<-test_data_all[test_data_all$sex_adj==1,]}
  if(k=="women"){test_data<-test_data_all[test_data_all$sex_adj==0,]}
  if(k=="young"){test_data<-test_data_all[test_data_all$ages<65,]}
  if(k=="old"){test_data<-test_data_all[!test_data_all$ages<65,]}
  
  attach(test_data)
  
  #genetically predicted ldl and tg using all snps and relevant weights#
  ldl_fit = lm(ldl~ldltg_ldl_sum
               +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
               +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
               +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)$fitted
  
  tg_fit = lm(tg~ldltg_tg_sum
              +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
              +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
              +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)$fitted
  
  #chd#
  model = glm(chd~ldl_fit+tg_fit
              +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
              +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
              +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20
              , family=binomial)
  
  beta.ratio = summary(model)$coef[2,1]
  se.ratio = summary(model)$coef[2,2]
  p.ratio = summary(model)$coef[2,4]
  
  dat_tab<-data.frame(k,"ldl","chd","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","ldl_fit","tg_fit","k")])
  
  #dead#
  model = coxph(Surv(outagef,dead)~ldl_fit+tg_fit
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  beta.ratio = summary(model)$coef[1,1] 
  se.ratio = summary(model)$coef[1,3] 
  p.ratio = summary(model)$coef[1,5]
  
  dat_tab<-data.frame(k,"ldl","dead","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","ldl_fit","tg_fit","k")])
  
  #cv_dead#
  model = coxph(Surv(outagef,cv_dead)~ldl_fit+tg_fit
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  beta.ratio = summary(model)$coef[1,1] 
  se.ratio = summary(model)$coef[1,3] 
  p.ratio = summary(model)$coef[1,5]
  
  dat_tab<-data.frame(k,"ldl","cv_dead","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","ldl_fit","tg_fit","k")])
  
  #tumo_dead#
  model = coxph(Surv(outagef,tumo_dead)~ldl_fit+tg_fit
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  beta.ratio = summary(model)$coef[1,1] 
  se.ratio = summary(model)$coef[1,3] 
  p.ratio = summary(model)$coef[1,5]
  
  dat_tab<-data.frame(k,"ldl","tumo_dead","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","ldl_fit","tg_fit","k")])
  
  #other_dead#
  model = coxph(Surv(outagef,other_dead)~ldl_fit+tg_fit
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  beta.ratio = summary(model)$coef[1,1] 
  se.ratio = summary(model)$coef[1,3] 
  p.ratio = summary(model)$coef[1,5]
  
  dat_tab<-data.frame(k,"ldl","other_dead","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","ldl_fit","tg_fit","k")])
  
  detach(test_data)
  rm(test_data,ldl_fit,tg_fit,k)}

rm(test_data_all)
gc()
####Two state least squares (tg_adj)####
test_data_all<-fread("D:/MR/LDL_CAD/R/Revise3/test_data_all.csv")

test_data_all$apob<-RankNorm(test_data_all$apob)
test_data_all$ldl<-RankNorm(test_data_all$ldl)
test_data_all$tg<-RankNorm(test_data_all$tg)

#analyses#
for (k in c('overall',"men","women","young","old")) {
  if(k=="overall"){test_data<-test_data_all}
  if(k=="men"){test_data<-test_data_all[test_data_all$sex_adj==1,]}
  if(k=="women"){test_data<-test_data_all[test_data_all$sex_adj==0,]}
  if(k=="young"){test_data<-test_data_all[test_data_all$ages<65,]}
  if(k=="old"){test_data<-test_data_all[!test_data_all$ages<65,]}
  
  attach(test_data)
  
  #genetically predicted apob and tg using all snps and relevant weights#
  apob_fit = lm(apob~apobtg_apob_sum
                +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)$fitted
  
  tg_fit = lm(tg~apobtg_tg_sum
              +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
              +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
              +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)$fitted
  
  #chd#
  model = glm(chd~tg_fit+apob_fit
              +ages+ages2+sex_adj+ages*sex_adj+ages2*sex_adj
              +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
              +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20
              , family=binomial)
  
  beta.ratio = summary(model)$coef[2,1]
  se.ratio = summary(model)$coef[2,2]
  p.ratio = summary(model)$coef[2,4]
  
  dat_tab<-data.frame(k,"tg","chd","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","apob_fit","tg_fit","k")])
  
  #dead#
  model = coxph(Surv(outagef,dead)~tg_fit+apob_fit
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  beta.ratio = summary(model)$coef[1,1] 
  se.ratio = summary(model)$coef[1,3] 
  p.ratio = summary(model)$coef[1,5]
  
  dat_tab<-data.frame(k,"tg","dead","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","apob_fit","tg_fit","k")])
  
  #cv_dead#
  model = coxph(Surv(outagef,cv_dead)~tg_fit+apob_fit
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  beta.ratio = summary(model)$coef[1,1] 
  se.ratio = summary(model)$coef[1,3] 
  p.ratio = summary(model)$coef[1,5]
  
  dat_tab<-data.frame(k,"tg","cv_dead","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","apob_fit","tg_fit","k")])
  
  #tumo_dead#
  model = coxph(Surv(outagef,tumo_dead)~tg_fit+apob_fit
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  beta.ratio = summary(model)$coef[1,1] 
  se.ratio = summary(model)$coef[1,3] 
  p.ratio = summary(model)$coef[1,5]
  
  dat_tab<-data.frame(k,"tg","tumo_dead","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","apob_fit","tg_fit","k")])
  
  #other_dead#
  model = coxph(Surv(outagef,other_dead)~tg_fit+apob_fit
                +birth+birth2+sex_adj+birth*sex_adj+birth2*sex_adj
                +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
                +PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20)
  
  beta.ratio = summary(model)$coef[1,1] 
  se.ratio = summary(model)$coef[1,3] 
  p.ratio = summary(model)$coef[1,5]
  
  dat_tab<-data.frame(k,"tg","other_dead","GRS (adj)",beta.ratio,se.ratio,beta.ratio-1.96*se.ratio,beta.ratio+1.96*se.ratio,
                      p.ratio,NA,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                             ,"test_data_all","test_data","apob_fit","tg_fit","k")])
  
  detach(test_data)
  rm(test_data,apob_fit,tg_fit,k)}

rm(test_data_all)
gc()

####NLMR (apob)####
test_data_all<-fread("D:/MR/LDL_CAD/R/Revise3/test_data_all.csv")

#analyses#
for (k in c('overall',"men","women","young","old")) {
  if(k=="overall"){test_data<-test_data_all}
  if(k=="men"){test_data<-test_data_all[test_data_all$sex_adj==1,]}
  if(k=="women"){test_data<-test_data_all[test_data_all$sex_adj==0,]}
  if(k=="young"){test_data<-test_data_all[test_data_all$ages<65,]}
  if(k=="old"){test_data<-test_data_all[!test_data_all$ages<65,]}
  
  #format data chd#
  summ_data<-create_nlmr_summary(y = test_data$chd,
                                 x = test_data$apob,
                                 g = test_data$apob_sum,
                                 covar = matrix(data=c(test_data$ages,
                                                       test_data$ages2,
                                                       test_data$sex_adj,
                                                       test_data$ages*test_data$sex_adj,
                                                       test_data$ages2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref = 1)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"apob","chd","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"apob","chd","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  source("D:/MR/LDL_CAD/R/corrected_sumnlmr_surv.R")
  #format data dead#
  summ_data<-create_nlmr_summary(y = Surv(test_data$outagef,test_data$dead),
                                 x = test_data$apob,
                                 g = test_data$apob_sum,
                                 covar = matrix(data=c(test_data$birth,
                                                       test_data$birth2,
                                                       test_data$sex_adj,
                                                       test_data$birth*test_data$sex_adj,
                                                       test_data$birth2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref=1)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"apob","dead","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"apob","dead","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  #format data cv_dead#
  summ_data<-create_nlmr_summary(y = Surv(test_data$outagef,test_data$cv_dead),
                                 x = test_data$apob,
                                 g = test_data$apob_sum,
                                 covar = matrix(data=c(test_data$birth,
                                                       test_data$birth2,
                                                       test_data$sex_adj,
                                                       test_data$birth*test_data$sex_adj,
                                                       test_data$birth2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref=1)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"apob","cv_dead","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"apob","cv_dead","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  #format data tumo_dead#
  summ_data<-create_nlmr_summary(y = Surv(test_data$outagef,test_data$tumo_dead),
                                 x = test_data$apob,
                                 g = test_data$apob_sum,
                                 covar = matrix(data=c(test_data$birth,
                                                       test_data$birth2,
                                                       test_data$sex_adj,
                                                       test_data$birth*test_data$sex_adj,
                                                       test_data$birth2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref=1)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"apob","tumo_dead","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"apob","tumo_dead","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  #format data other_dead#
  summ_data<-create_nlmr_summary(y = Surv(test_data$outagef,test_data$other_dead),
                                 x = test_data$apob,
                                 g = test_data$apob_sum,
                                 covar = matrix(data=c(test_data$birth,
                                                       test_data$birth2,
                                                       test_data$sex_adj,
                                                       test_data$birth*test_data$sex_adj,
                                                       test_data$birth2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked"
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref=1)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"apob","other_dead","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"apob","other_dead","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3","test_data_all")])}

rm(test_data_all)

####NLMR (ldl)####
test_data_all<-fread("D:/MR/LDL_CAD/R/Revise3/test_data_all.csv")

#analyses#
for (k in c('overall',"men","women","young","old")) {
  if(k=="overall"){test_data<-test_data_all}
  if(k=="men"){test_data<-test_data_all[test_data_all$sex_adj==1,]}
  if(k=="women"){test_data<-test_data_all[test_data_all$sex_adj==0,]}
  if(k=="young"){test_data<-test_data_all[test_data_all$ages<65,]}
  if(k=="old"){test_data<-test_data_all[!test_data_all$ages<65,]}
  
  #format data chd#
  summ_data<-create_nlmr_summary(y = test_data$chd,
                                 x = test_data$ldl,
                                 g = test_data$ldl_sum,
                                 covar = matrix(data=c(test_data$ages,
                                                       test_data$ages2,
                                                       test_data$sex_adj,
                                                       test_data$ages*test_data$sex_adj,
                                                       test_data$ages2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref = 3.5)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"ldl","chd","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"ldl","chd","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  source("D:/MR/LDL_CAD/R/corrected_sumnlmr_surv.R")
  #format data dead#
  summ_data<-create_nlmr_summary(y = Surv(test_data$outagef,test_data$dead),
                                 x = test_data$ldl,
                                 g = test_data$ldl_sum,
                                 covar = matrix(data=c(test_data$birth,
                                                       test_data$birth2,
                                                       test_data$sex_adj,
                                                       test_data$birth*test_data$sex_adj,
                                                       test_data$birth2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref = 3.5)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"ldl","dead","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"ldl","dead","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  #format data cv_dead#
  summ_data<-create_nlmr_summary(y = Surv(test_data$outagef,test_data$cv_dead),
                                 x = test_data$ldl,
                                 g = test_data$ldl_sum,
                                 covar = matrix(data=c(test_data$birth,
                                                       test_data$birth2,
                                                       test_data$sex_adj,
                                                       test_data$birth*test_data$sex_adj,
                                                       test_data$birth2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref = 3.5)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"ldl","cv_dead","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"ldl","cv_dead","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  #format data tumo_dead#
  summ_data<-create_nlmr_summary(y = Surv(test_data$outagef,test_data$tumo_dead),
                                 x = test_data$ldl,
                                 g = test_data$ldl_sum,
                                 covar = matrix(data=c(test_data$birth,
                                                       test_data$birth2,
                                                       test_data$sex_adj,
                                                       test_data$birth*test_data$sex_adj,
                                                       test_data$birth2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref = 3.5)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"ldl","tumo_dead","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"ldl","tumo_dead","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  #format data other_dead#
  summ_data<-create_nlmr_summary(y = Surv(test_data$outagef,test_data$other_dead),
                                 x = test_data$ldl,
                                 g = test_data$ldl_sum,
                                 covar = matrix(data=c(test_data$birth,
                                                       test_data$birth2,
                                                       test_data$sex_adj,
                                                       test_data$birth*test_data$sex_adj,
                                                       test_data$birth2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref = 3.5)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"ldl","other_dead","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"ldl","other_dead","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3","test_data_all")])}

rm(test_data_all)

####NLMR (tg)####
test_data_all<-fread("D:/MR/LDL_CAD/R/Revise3/test_data_all.csv")

#analyses#
for (k in c('overall',"men","women","young","old")) {
  if(k=="overall"){test_data<-test_data_all}
  if(k=="men"){test_data<-test_data_all[test_data_all$sex_adj==1,]}
  if(k=="women"){test_data<-test_data_all[test_data_all$sex_adj==0,]}
  if(k=="young"){test_data<-test_data_all[test_data_all$ages<65,]}
  if(k=="old"){test_data<-test_data_all[!test_data_all$ages<65,]}
  
  #format data chd#
  summ_data<-create_nlmr_summary(y = test_data$chd,
                                 x = test_data$tg,
                                 g = test_data$tg_sum,
                                 covar = matrix(data=c(test_data$ages,
                                                       test_data$ages2,
                                                       test_data$sex_adj,
                                                       test_data$ages*test_data$sex_adj,
                                                       test_data$ages2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref = 1.5)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"tg","chd","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"tg","chd","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  source("D:/MR/LDL_CAD/R/corrected_sumnlmr_surv.R")
  #format data dead#
  summ_data<-create_nlmr_summary(y = Surv(test_data$outagef,test_data$dead),
                                 x = test_data$tg,
                                 g = test_data$tg_sum,
                                 covar = matrix(data=c(test_data$birth,
                                                       test_data$birth2,
                                                       test_data$sex_adj,
                                                       test_data$birth*test_data$sex_adj,
                                                       test_data$birth2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref = 1.5)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"tg","dead","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"tg","dead","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  #format data cv_dead#
  summ_data<-create_nlmr_summary(y = Surv(test_data$outagef,test_data$cv_dead),
                                 x = test_data$tg,
                                 g = test_data$tg_sum,
                                 covar = matrix(data=c(test_data$birth,
                                                       test_data$birth2,
                                                       test_data$sex_adj,
                                                       test_data$birth*test_data$sex_adj,
                                                       test_data$birth2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref = 1.5)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"tg","cv_dead","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"tg","cv_dead","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  #format data tumo_dead#
  summ_data<-create_nlmr_summary(y = Surv(test_data$outagef,test_data$tumo_dead),
                                 x = test_data$tg,
                                 g = test_data$tg_sum,
                                 covar = matrix(data=c(test_data$birth,
                                                       test_data$birth2,
                                                       test_data$sex_adj,
                                                       test_data$birth*test_data$sex_adj,
                                                       test_data$birth2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref = 1.5)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"tg","tumo_dead","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"tg","tumo_dead","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  #format data other_dead#
  summ_data<-create_nlmr_summary(y = Surv(test_data$outagef,test_data$other_dead),
                                 x = test_data$tg,
                                 g = test_data$tg_sum,
                                 covar = matrix(data=c(test_data$birth,
                                                       test_data$birth2,
                                                       test_data$sex_adj,
                                                       test_data$birth*test_data$sex_adj,
                                                       test_data$birth2*test_data$sex_adj,
                                                       test_data$PC1,
                                                       test_data$PC2,
                                                       test_data$PC3,
                                                       test_data$PC4,
                                                       test_data$PC5,
                                                       test_data$PC6,
                                                       test_data$PC7,
                                                       test_data$PC8,
                                                       test_data$PC9,
                                                       test_data$PC10,
                                                       test_data$PC11,
                                                       test_data$PC12,
                                                       test_data$PC13,
                                                       test_data$PC14,
                                                       test_data$PC15,
                                                       test_data$PC16,
                                                       test_data$PC17,
                                                       test_data$PC18,
                                                       test_data$PC19,
                                                       test_data$PC20
                                 ),ncol=25),
                                 family = "binomial",
                                 controlsonly = FALSE, 
                                 q = 10,
                                 strata_method= "ranked" 
                                 # ,report_het = TRUE
  ) 
  
  # fitting a fractional polynomial model
  model<- with(summ_data$summary, frac_poly_summ_mr(bx=bx,
                                                    by=by, 
                                                    bxse=bxse, 
                                                    byse=byse, 
                                                    xmean=xmean,
                                                    family="binomial",
                                                    fig=TRUE,
                                                    ref = 1.5)
  )
  
  summary<-summary(model)
  dat_tab2<-data.frame(k,"tg","other_dead","ranked",summary$powers,summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,summary$coefficients$pval,
                       summary$p_tests$fp_d1_d2,summary$p_tests$fp,summary$p_tests$quad,summary$p_tests$Q)
  colnames(dat_tab2)<-c("sex","exp","out","method","powers","beta","se","lo","up","pval",
                        "pval(degree)","pval(non-linear)","pval(quadratic)","pval(cochran_Q)")
  results2<-rbind(results2,dat_tab2)
  
  #fitting a piecewise linear model#
  model2 <-with(summ_data$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                     ci="bootstrap_se",
                                                     nboot=1000, 
                                                     fig=FALSE,
                                                     family="binomial",
                                                     ci_fig="ribbon")
  )
  summary<-summary(model2)
  
  dat_tab3<-data.frame(k,"tg","other_dead","ranked",summ_data$summary$bx,summ_data$summary$xmean,
                       summary$coefficients$beta,summary$coefficients$se,
                       summary$coefficients$lci,summary$coefficients$uci,
                       summary$coefficients$pval)
  colnames(dat_tab3)<-c("sex","exp","out","method","bx","mean","beta","se","lo","up","pval")
  results3<-rbind(results3,dat_tab3)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","results3","test_data_all")])}

rm(test_data_all)
