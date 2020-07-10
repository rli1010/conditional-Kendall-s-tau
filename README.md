# conditional-Kendall-s-tau
R codes for calculating the conditional Kendall's tau with bivariate survival data
Code developed by Xiangyu Liu and Ruosha Li

Liu, X, Ning, J, Cheng, Y, Huang, X, Li, R. A flexible and robust method for assessing conditional association and conditional concordance. Statistics in Medicine. 2019; 38: 3656– 3668.

#####################

path='H:\\bivKendall\\'
##the .dll file is for 64-bit windows##
dyn.load(paste(path,"quantKendall.dll",sep=''))

##analysis for uncensored data##
##conducts 200 bootstrap by default##
source(paste(path,"quantKen_nocensor.r",sep=''))
simda_comp=read.csv(paste(path,"simda_comp.csv",sep=''))
quantKendall.comp(data=simda.comp,Y1~Z1+Z2,Y2~Z1+Z2)
