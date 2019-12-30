blwhf<-function(Ts,Ta,rh,Fc){
# - RPT: ecluded switch, default='berliand'
# - reduced input-vrbls to 4, using bulkf_default
# - included 'satvap(Ta)' into script
#
# BLWHF: estimates net long-wave heat flux using bulk formulas
# lwest = BLWHF(Ts,Ta,rh,Fc,bulk_f) estimates the net longwave 
# heat flux into the ocean using one of a number of bulk formulas 
# evaluated by Fung et al (1984), Rev. Geophys. Space Physics, 
# 22,177-193. 
#
# INPUT:   Ts - sea surface temperature [C]
#          Ta - air temperature  [C]
#          rh - relative humidity  [#]
#          Fc - cloudiness correction factor F(C) (=1 for clear sky) 
#                  (see cloudcor.m for more information)
#          bulk_f - bulk formulas to be used 
#                  (see Fung et al, 1984) for details). 
#                   Options are:
#                     'brunt'
#                     'berliand'  - probably the best one (default)
#                     'clark'
#                     'hastenrath'
#                     'efimova'
#                     'bunker'
#                     'anderson'
#                     'swinbank' - does not use rh
#
# OUTPUT:  lwest - net downward longwave flux [W/m^2]

######################################################################
# 8/31/98: version 1.1 (contributed by RP)
# 8/5/99: version 2.0
######################################################################


# compute vapour pressure from relative humidity (using formulas 
# from Gill, 1982)

#function 'satvap(Ta)' embedded here
ew<-10^((0.7859+0.03477*Ta)/(1+0.00412*Ta))
fw<-1 + 1e-6*P_default*(4.5+0.0006*Ta^2)
ew<-fw*ew

rw <- eps_air*ew/(P_default-ew)
r <- (rh/100)*rw
e_a <- r*P_default/(eps_air+r)

# convert to K
ta <- Ta+CtoK
ts <- Ts+CtoK
lwest <- -emiss_lw*sigmaSB*ta^4*(0.39 - 0.05*sqrt(e_a))*Fc - 4*emiss_lw*sigmaSB*ta^3*(Ts - Ta)

return(lwest)
}