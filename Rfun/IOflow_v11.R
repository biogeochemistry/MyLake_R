# === MyLake model, version 1.1, 16.03.04 ===
# by Tom Andersen & Tuomo Saloranta, NIVA 2004
#
# Module for calculating river inflow and its effect on profiles of temperature and other state variables
# Code checked by TSA, 16.3.2004


# INPUTS:
#   vertical arrays:    z is model grid, Vz is layer volume, Tz layer property (e.g. temperature);
#   scalars:            lvlD is the grid depth level *above* which the inflow settles (m),
#                       Iflw is the inflow volume (m3/day), and
#                       T_Iflw is the property of the inflow, and dz the grid stepsize (layer thickness)
# OUTPUTS: Cz is the new property profile after inflow

IOflow_v11<-function(dz, z, Vz, Tz, lvlD, Iflw, T_Iflw){

dum<-which(z>=lvlD) #same as lvlD but in grid level numbers
lvlG<-dum[1]
Cz<-Tz

if (lvlD==0){ #if inflow is lighter than surface water, mix it with first layer
Cz[1] <- (Vz[1]*Tz[1]+Iflw*T_Iflw)/(Vz[1]+Iflw);

}else{    #otherwise add the inflow to the appropriate
        #depth level and "lift" the water column above

 Trev_pour<-c(Tz[1:lvlG-1], T_Iflw, 0)[(1+lvlG):1]   #up-side down Tz to be poured in (plus an extra zero)       
 Vzrev_pour<-c(Vz[1:lvlG-1], Iflw, 0)[(1+lvlG):1]    #up-side down Vz to be poured in (plus an extra zero)     
 Vzrev_fill<-c(Vz[1:lvlG-1])[(lvlG-1):1]             #up-side down Vz to be filled in

 for (i in 1:length(Vzrev_fill)){
  inxA<-which(cumsum(Vzrev_pour)>Vzrev_fill[i])
  inxB<-inxA[1] #index of the last layer to be poured in (partly)
  ShakerV<-c(Vzrev_pour[1:(inxB-1)], Vzrev_fill[i]-sum(Vzrev_pour[1:(inxB-1)]))
  ShakerT<-c(Trev_pour[1:inxB])
  Cz[lvlG-i]<-sum(ShakerV*ShakerT)/sum(ShakerV) #new property after mixing
  Vzrev_pour[1:inxB]<-Vzrev_pour[1:inxB]-ShakerV #subtract poured volumes from the "reserves"
 } # loop
} # if...

return(Cz)
}