psittc<-function(zet){
# PSITTC: computes potential temperature profile following TOGA/COARE.
# y=PSITTC(zet) computes the turbulent potential temperature profile 
# function given zet = (z/L), L the Monin-Obukoff length scale, following 
# Edson and Fairall TOGA COARE code (version 2.0), as modified to use
# Rogers' weighting factor to combine the Dyer and free convective 
# forms for unstable conditions. 

######################################################################
# 8/28/98: version 1.1
# 8/5/99: version 1.2
######################################################################

c13 <- 1.0/3.0
sq3 <- sqrt(3.0)

# stable conditions
y <- -4.7*zet

# unstable conditions
j <- which(zet<0)
zneg <- zet[j]

# nearly stable (standard functions)
 x <- (1-16.0*zneg)^0.25
 y1 <- 2.0*log((1+x^2)/2)

# free convective limit
 x <- (1-12.87*zneg)^c13
 y2 <- 1.5*log((x^2+x+1)/3.0) - sq3*atan((2*x+1)/sq3) + pi/sq3

# weighted sum of the two
 F <- 1.0/(1+zneg^2)
 y[j] <- F*y1 + (1-F)*y2
return(y)
}