# === MyLake model, version 1.1, 16.03.04 ===
# by Tom Andersen & Tuomo Saloranta, NIVA 2004
#
# Module for calculating new temperature profile in sediments
# Code checked by TSA, 16.03.2004

sedimentheat_v11<-function(Tzy_sed, K_sed, dt){

# Grid resolution in sediment is 0.2 m in the 0-2 layer and 0.5 m in the 2-10 m layer 
# (totally 10 + 16 = 26 layers)
dz_s <- c(NA, rep(0.2, 9), 0.35, rep(0.5,15), NA); #length = 27

alpha_sed <-  dt/(dz_s^2)
Ns<-rep(NA,dim(Tzy_sed)[1])
Nw<-rep(NA,dim(Tzy_sed)[2])

Tzy_sed_updated<-matrix(0,dim(Tzy_sed))

for (i in 1:Nw){ #do for every water layer in contact with sediment
    
   az <- K_sed * alpha_sed[1:(length(K_sed)-1)];
   bz <- K_sed * alpha_sed[2:length(K_sed)]
   az[1] <- 0; #top sediment boundary condition
   bz[1] <- 0; #top sediment boundary condition
   bz[length(bz)] <- 0; #deep sediment boundary condition
   
   Gi <- c(-bz, (1 + az + bz), -az)
   Fi <- spdiags(Gi,-1:1,Ns,Ns)';
   Tzy_sed_updated(:,i) = Fi \ Tzy_sed(:,i);
return [Tzy_sed_updated] 