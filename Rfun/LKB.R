LKB<-function(Reu){
# LKB: computes rougness Reynolds numbers for temperature and humidity
# [Ret,Req]=LKB(Reu) computes the roughness Reynolds for temperature
# and humidity following Liu, Katsaros and Businger (1979), J. Atmos.
# Sci., 36, 1722-1735.

######################################################################
# 8/28/98: version 1.1
# 8/5/99: version 1.2
######################################################################

   Ret <- Reu; Ret[]<-1
   Req <- .292*Ret[]
j <- which(Reu>.11 & Reu<=.825)
   Ret[j] <- 1.376*Reu[j]^0.929
   Req[j] <- 1.808*Reu[j]^0.826
j <- which(Reu>.825 & Reu<=3)
   Ret[j] <- 1.026/Reu[j]^0.599
   Req[j] <- 1.393/Reu[j]^0.528
j <- which(Reu>3 & Reu<=10)
   Ret[j] <- 1.625/Reu[j]^1.018
   Req[j] <- 1.956/Reu[j]^0.870
j <- which(Reu>10 & Reu<=30)
   Ret[j] <- 4.661/Reu[j]^1.475
   Req[j] <- 4.994/Reu[j]^1.297
j <- which(Reu>30)
   Ret[j] <- 34.904/Reu[j]^2.067
   Req[j] <- 30.790/Reu[j]^1.845
return(list(Ret, Req))
}