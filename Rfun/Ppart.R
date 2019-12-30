Ppart<-function(vf,TIP,Psat,Fmax,rho_sed,Fstable){
# Function for calculating the partitioning between
# dissolved and inorganic particle bound phosphorus. 
# Based on Langmuir isotherm approach 
#vf:    volume fraction of suspended inorganic matter (m3 m-3); S/rho_sed OR (1-porosity)
#TIP:   Total inorganic phosphorus (mg m-3)
#Psat, mg m-3 - Langmuir half-saturation parameter
#Fmax, mg kg-1  - Langmuir scaling parameter
#rho_sed, kg m-3 - Density of dry inorganic sediment mass
#Fstable, mg kg-1 - Inactive P conc. in inorg. particles


N<-length(TIP)

Pdiss <- rep(NA, l=N)

for (w in 1:N){
 a <- vf[w]-1
 b <- TIP[w] + (vf[w]-1)*Psat - vf[w]*rho_sed*(Fmax+Fstable)
 c <- Psat*TIP[w] - vf[w]*rho_sed*Fstable*Psat 
Pdiss[w] <-  suppressWarnings(max(as.numeric(1/polyroot(c(a, b, c))))) # CHECK!!!
}

Pfpart = (TIP - (1-vf)*Pdiss)/(rho_sed*vf) #inorg. P conc. in sediment particles(mg kg-1 dry w.) 

return(list(Pdiss, Pfpart))
}