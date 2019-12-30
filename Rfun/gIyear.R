# Functions to convert degrees to radians and vice versa.
# These two functions are not used in gIyear. Use'em as you like.
d2r <- function(deg){ (pi/180)*deg }
r2d <- function(rad){ (180/pi)*rad }

# If Tlk is not supplied, defaults are used (splinefunction).
# If Tlk is supplied, it must be for each day in a year (or sequence in question).
# Daynr sequence is per default set to 365 days with steps of 1 day.
# Hour sequence is per default set to 0-24 hours with a step of 10 minutes.
# Latitude must be given in radians.
# z is meters above sea level.
gIyear <- function(lat, z, daynr = 1:365, hseq = seq(0,24,1/6), Tlk = NULL){
 if(length(lat) != 1 || length(z) != 1 ) stop("Expected lat and z to be of length 1.")
 if(is.null(Tlk)){
   message("--- Using the default monthly Tlk values. ---")
   message("(2.2,2.4,3.0,3.2,3.5,3.6,3.8,3.9,3.5,3.0,2.4,2.2)")

   # Slightly modified from the SoDa database to get a good polynomial fit.
   Tlk <- c(2.2,2.4,3.0,3.2,3.5,3.6,3.8,3.9,3.5,3.0,2.4,2.2)
   year <- as.numeric(format(Sys.Date(),"%Y"))
   x <- seq(ISOdate(year,1,15), ISOdate(year,12,15), by="month")
   x <- as.numeric(x - ISOdate(year,1,1))
   Tlk.spline <- splinefun(x,Tlk, method="periodic")
   Tlk <- Tlk.spline(daynr)
   rm(year, x, Tlk.spline)
 }

 # -- Beam irradiance --
 # Day angle in radians given the day number (1-365(366))
 j <- (2*pi*daynr)/365.25
 # Correction factor for the eccentricity of the earth orbit
 e <- 1 + 0.03344 * cos(j - 0.048869)
 # Space irradiance
 sI <- 1367 * e
 # Solar declination per day over a year.
 declin <- asin(0.3978*sin(j-1.4+(0.0355*sin(j-0.0489))))
 rm(e,j)
 # Hour angle over a whole day, repeated over length(daynr).
 ha <- 0.261799*(hseq-12) # Hour angle from midnight to midnight
 ha <- matrix(ha, nrow=length(daynr), ncol=length(hseq), byrow=TRUE)
 # Sun's height above the horizon.
 sh <- cos(lat)*cos(declin)*cos(ha) + sin(lat)*sin(declin)
 ash <- asin(sh)
 # Relative optical air mass
 dh0 <- 0.061359*((0.1594+1.123*ash+0.065656*ash^2)/(1+28.9344*ash+277.3971*ash^2))
 h0ref <- ash+dh0
 m <- exp(-z/8434.5)/(sin(h0ref)+0.50572*(h0ref+6.07995)^-1.6364)
#--- TO DO: m can't be negative. Find out where to truncate and which values to set.
 # Relative optical thickness
 hrows <- which(m > 20)
 lrows <- which(m <= 20)
 ot <- matrix(nrow=dim(m)[1], ncol=dim(m)[2])
 ot[lrows] <- 1/(6.6296+1.7513*m[lrows]-0.1202*m[lrows]^2+0.0065*m[lrows]^3-0.00013*m[lrows]^4)
 ot[hrows] <- 1/(10.4+0.718*m[hrows])
 # Calculate beam irradiance on a horizontal surface.
 bI <- sI * exp(-0.8662*Tlk*m*ot) * sh
 bI[bI<0] <- 0 # Truncating negative values to zero.
 #return(bI)
 #--- END beam irradiance ---

 #--- Diffuse irradiance ---
 # Transmission function
 tn <- -0.015843+0.030543*Tlk+0.0003797*Tlk^2
 # Solar altitude function
 A1 <- 0.26463-0.061581*Tlk+0.0031408*Tlk^2
 x <- A1 * tn
 A1[x<0.0022] <- 0.0022/tn[x<0.0022]
 A2 <- 2.04020+0.018945*Tlk-0.011161*Tlk^2
 A3 <- -1.3025+0.039231*Tlk+0.0085079*Tlk^2
 fd <- A1 + A2*sh + A3*sh^2
 rm(A2, A3)# RPT: keep A1
 # diffuse irradiance
 dI <- sI * tn * fd
 dI[dI<0] <- 0
 #--- END diffuse irradiance ---

 #--- Global irradiance ---
 gI <- bI + dI
 #--- End global irradiance ---
 SoAl<-60*asin(sh) # RPT solar altitude [deg] ??
 return(list(sh, gI))# RPT

}
