cloudcor<-function(C,lat){
# RPT: excluded options, 'bunker' by default
# CLOUDCOR: computes cloud correction factor for bulk long-wave flux.
# Fc=CLOUDCOR(C,optns,lat) computes the cloud correction factor F(C)
# as a function of the cloud fraction C for bulk long-wave flux formulae.
# In general, these are functions of the form
#             1 - a_n*C^n
# Since the coefficients and powers depend a lot on the dominant cloud
# type which may vary from region to region and season to season, it is
# not clear which parametrization is best (see Fung et al (1984), 
# Rev. of Geophys. and Space Phys., 22, 177-193). 
#
# The particular parametrization used here depends on the second input
# variable, for which no default is given to emphasize the fact that you
# really need to understand what you are doing here!
#
# optns = [a1 a2] = use a correction factor of [1-a1*C-a2*C^2].
#
# There are several "built-in" formulae (from Fung et al) that all have
# a latitude-dependence of some kind.
#
# optns = 'clarke',lat = Clarke (1974) corrections for abs(lat)<50.
#       = 'bunker',lat = Bunker (1976) corrections for N Atlantic.
#
# INPUT:   C - cloud fraction
#          optns - see above for details
#          lat - latitude [deg] - required for "built-in" formulae only
#
# OUTPUT:  Fc - correction factor used as input to BLWHF
######################################################################
# 3/12/98: version 1.1 (contributed by RP)
# 8/5/99: version 2.0
######################################################################

	ls1<-c(0, 7, 15, 25, 35, 45, 55, 65, 75)
	ls2<-c(0.5,0.52,0.59,0.63,0.68,0.72,0.76,0.8,0.84)
	a1<-ls2[which.max(1/(-ls1+lat))]
#      if(abs(lat)>75) a1 <- 0.84
#      if(abs(lat)>65) a1 <- 0.80
#      if(abs(lat)>55) a1 <- 0.76
#      if(abs(lat)>45) a1 <- 0.72
#      if(abs(lat)>35) a1 <- 0.68
#      if(abs(lat)>25) a1 <- 0.63
#      if(abs(lat)>15) a1 <- 0.59
#      if(abs(lat)>7)  a1 <- 0.52
#      else            a1 <- 0.50

Fc  <-  1 - a1*C
return(Fc)
}