rm(list=ls())
source("./Rfun/req_packs.R")
source("./Rfun/req_data.R")
source("./Rfun/req_funs.R")
source("./Rfun/as_consts.R")

# tt time vector with input data (MATLAB format)
# transform from matlab to R with 'chron' this way:
# chron(tt[i], origin=c(month = 12, day = 31, year =1999))

M_start <- tt[1] 
M_stop <- tt[length(tt)]
res<-Rlake_main(M_start, M_stop)

#
# output of function: zz,Az,Vz,tt,Qst,Kzt,Tzt,Czt,Szt,Pzt,Chlzt,
# PPzt,DOPzt,DOCzt,Qzt_sed,lambdazt, P3zt_sed,P3zt_sed_sc,
# His,DoF,DoM,MixStat,Wt,w_chl_zt
#
# obtain results stored in 'res' by appending the object names on 'res'
# example:  '>res$Tzt'

xx<-chron(sort(rep(M_start:M_stop, 20)))
yy<-rep(0:19, length(M_start:M_stop))
rbw<-rainbow(30, start=0, end=0.6)[30:1]
x11(8,10); par(mfrow=c(3, 1))
plot(xx,-yy, pch="'", cex=3, col=rbw[cut(res$Tzt, 30, lab=F)], main="Temperature profile", xlab="date", ylab="depth")
plot(xx,-yy, pch="'", cex=3, col=rbw[cut(res$Chlzt, 30, lab=F)], main="Chl-a profile", xlab="date", ylab="depth")
plot(xx,-yy, pch="'", cex=3, col=rbw[cut(res$Pzt, 30, lab=F)], main="PO4 profile", xlab="date", ylab="depth")
