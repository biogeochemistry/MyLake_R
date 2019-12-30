heatflux_v12<-function(T,GR,CC,Ta,Ra,Pa,Ua,Tw,lat,lon,Hs,Hi,alb_melt_ice,alb_melt_snow,alb_table){

#RPT: function 'stressve' implemented; cdnve now called from 'heatflux'
# === MyLake model, version 1.2, 15.03.05 ===
# by Tom Andersen & Tuomo Saloranta, NIVA 2005
#
# Module for calculating heat fluxes and other physical parameters 
# Code checked by TSA, 03.03.2005
# Last modified by TSA, 03.03.2005

# Daily water surface heat balance fluxes
# 		Based on MATLAB Air-Sea Toolbox
#
# Output fluxes:
# Qsw	: Shortwave (solar) radiation flux (W m-2)
# Qlw	: Net longwave radiation flux (W m-2)
# Qsl	: Net sensible and latent heat flux (W m-2)
# tau	: Wind stress (N m-2)
# Dayfrac	: Fraction of daylight (-)
# Dayfracheating	: Fraction of day when sun is above a threshold angle
# (alt_trsh)  (-); used to partition heat flux calculations into night and daytime
#
#
# Input data:
# t   : time 
# GR  : Global radiation (MJ m-2 day-1))
# CC  : Fraction cloud cover (-)
# Ta  : Air temperature (C, at 2 m height)
# Pa  : Air pressure (mbar==hPa)
# Ra  : Relative air humidity (#, at 2 m height)
# Ua  : Wind speed (m/s, at 10 m height)
# Tw  : Surface temperature (C)
# lat : Latitude (decimal degrees)
# lon : Longitude (decimal degrees)
# Hs  : Snow water equivalent (m)
# Hi  : Ice thickness (m)
# alb_melt_ice : Albedo for melting ice
# alb_melt_snow : Albedo for melting snow

# Air-Sea Toolbox: sensible + latent heat flux

A <- hfbulktc_speed(Ua,10,Ta,2,Ra,2,Pa,Tw) #function modified (speed-up) by TSA
Qsl <- A[[1]] + A[[2]]

#    Eevapor<-Eevapor+A[[2]]; #testing evaporation amounts# Evapor=0 (from main call)

# Air_sea toolbox: long-wave radiation flux
Qlw <- blwhf(Tw,Ta,Ra,cloudcor(CC,lat))
    
# Air_sea toolbox: short-wave radiation flux
# Albedo corrected, measured solar radiation
    
#dv<-chron(731276, origin=c(month = 12, day = 31, year =1999)) # for testing; should be mar 1 2002
#dv<-chron(T, origin=c(month = 12, day = 31, year =1999))
dv<-T
yr <- dv[1]
resol_t<-24*60/5 #at how many times solar radiation is calculated per day (24*60/n; every n minutes)
obj<-soradna1(dv, resol_t, 0, lat); #alt in degrees
alt<-obj[[1]]; rad<-obj[[2]]

 inx<-which(alt<0)
 alt[inx]<-0 #transform negative altitudes to zeros
 Dayfrac <- 1 - length(inx)/length(alt) # Calculate fraction of daylight
 
 alt_trsh<-15
 inx2<-which(alt > alt_trsh) #find all altitude-angles above alt_trsh
 Dayfracheating<-length(inx2)/length(alt) #fraction of day when sun is above alt_trsh
 
 #Calculate cloud correction Reed (1977) (used if no Global Radiation is inputted)
 if (CC<0.3){CloudCorr<-1
 }else{CloudCorr<-(1-0.62*CC+0.0019*max(alt))}
 
if(is.na(GR)){   # if Global radiation input is missing calculate total transmissivity from formulaes 
    
      #Trnsmiss=0.95*sin(alt*pi/180)*CloudCorr # by Iqbal (1993??) and Reed (1977) (Varies over day, and zeros at night)
      Trnsmiss<- rep( (0.377 + 0.00513*max(alt))*CloudCorr, length=length(alt))
      # empirically determined from Breisjo & Vansjo data

}else{ #if Global radiation input is given calculate total transmissivity based on the observations  
    
    Qmm <- mean(rad)				# Calculated Max. solar radiation - no atmosphere
    Qmo <- (1000000/(24*60*60)) * GR	# Observed solar radiation (W/m^2)
    Trnsmiss<- rep(min((Qmo/Qmm),1), length(alt))   #Total atmospheric transmissivity (Constant over day & night)
        
}

if(Hs>0) {alb <- alb_melt_snow    
}else{if (Hi>0) {alb <- alb_melt_ice
}else alb <- albedo_mod(Trnsmiss, alt, alb_table)} #modified albedo function to save execution time


Qma <- (1 - alb) * rad * Trnsmiss


Qma[which(is.na(Qma))] <- 0

Qsw <- mean(Qma)

#wind stress (N/m2)
if (Ua>0){
obj<-cdnve(Ua, 10)
cd<-obj[[1]]; u10<-obj[[2]]
tau<-rho_air*(cd*u10^2)
}else{
tau<-0
}

return(list(Qsw,Qlw,Qsl,tau,Dayfrac,Dayfracheating))
}