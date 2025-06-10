# conditional-Kendall-s-tau
R codes for calculating the conditional Kendall's tau with bivariate continuous outcomes that are uncensored or subject to censoring

Example codes to two example datasets:<br><br>


##the .dll file is for 64-bit windows##<br><br>
#reset path to local file folder where the files are saved<br><br>

path='H:\\R code share\\bivKendall\\'  
dyn.load(paste(path,"quantKendall.dll",sep=''))<br>
library(quantreg)<br>

##analysis for uncensored data##<br>
##conducts 200 bootstrap by default##<br>

source(paste(path,"quantKendall.comp.R",sep=''))<br>
simda_comp=read.csv(paste(path,"simda_comp.csv",sep=''))<br>

fit1=quantKendall.comp(data=simda_comp,Y1 ~ Z1+Z2,Y2 ~ Z1+Z2)<br>
print(fit1)

###analysis for censored data###<br>
source(paste(path,"quantKendall.cens.R",sep=''))<br>
simda_cens=read.csv(paste(path,"simda_cens.csv",sep=''))<br>
c(mean(simda_cens$eta1),mean(simda_cens$eta2))<br>

fit2=quantKendall.cens(simda_cens, Surv(Y1,eta1)~Z1+Z2, Surv(Y2,eta2)~Z1+Z2, upper=c(0.8,0.8),n.bt=200)<br>
print(fit2)<br>




Code developed by Xiangyu Liu, Jing Ning, and Ruosha Li.<br>
Liu, X, Ning, J, Cheng, Y, Huang, X, Li, R. A flexible and robust method for assessing conditional association and conditional concordance. Statistics in Medicine. 2019; 38: 3656– 3668.

