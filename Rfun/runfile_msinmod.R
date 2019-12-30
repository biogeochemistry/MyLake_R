source("Rfun/rho.R")
source("Rfun/convection_v12.R")
require(R.matlab)
ies80<-c(6.536332e-9,-1.120083e-6,1.001685e-4,-9.09529e-3,6.793952e-2,999.842594)


% start main module
[zz,Az,Vz,tt,Qst,Kzt,Tzt,Czt,Szt,Pzt,Chlzt,PPzt,DOPzt,DOCzt,Qzt_sed,lambdazt,...
       P3zt_sed,P3zt_sed_sc,His,DoF,DoM,MixStat,Wt,w_chl_zt]...
          = simple_STEIN1(m_start,m_stop,initfile,'lake',inputfile,'timeseries', parafile,'lake'); 


# testing convection_v12
obj<-readMat("subdat_1.mat")
attach(obj)
# call conv-function
obj2<-convection_v12(Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz,Tprof.prev,
	Vz,Cw,f.par,lambdaz.wtot.avg,zz,swa.b0,tracer.switch,1)