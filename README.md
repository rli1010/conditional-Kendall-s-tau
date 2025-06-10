# conditional-Kendall-s-tau
R codes for calculating the conditional Kendall's tau with bivariate continuous outcomes that are uncensored or subject to censoring

Code developed by Xiangyu Liu, Jing Ning, and Ruosha Li.

Liu, X, Ning, J, Cheng, Y, Huang, X, Li, R. A flexible and robust method for assessing conditional association and conditional concordance. Statistics in Medicine. 2019; 38: 3656– 3668.

Example codes to two example datasets:

path='H:\\R code share\\bivKendall\\'  #reset to local file folder

##the .dll file is for 64-bit windows##
dyn.load(paste(path,"quantKendall.dll",sep=''))

library(quantreg)

##analysis for uncensored data##
##conducts 200 bootstrap by default##
source(paste(path,"quantKendall.comp.R",sep=''))
simda_comp=read.csv(paste(path,"simda_comp.csv",sep=''))

fit1=quantKendall.comp(data=simda_comp,Y1~Z1+Z2,Y2~Z1+Z2)
print(fit1)

###analysis for censored data###
source(paste(path,"quantKendall.cens.R",sep=''))
simda_cens=read.csv(paste(path,"simda_cens.csv",sep=''))
c(mean(simda_cens$eta1),mean(simda_cens$eta2))

fit2=quantKendall.cens(simda_cens, Surv(Y1,eta1)~Z1+Z2, Surv(Y2,eta2)~Z1+Z2, upper=c(0.8,0.8),n.bt=200)
print(fit2)
