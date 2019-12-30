albedo_mod<-function(trans,sunalt,albedot1){

# ALBEDO: computes sea surface albedo following Payne (1972).
# alb=ALBEDO(trans,sunalt) computes the sea surface albedo from the
# atmospheric transmittance and sun altitude by linear interpolation 
# using Table 1 in Payne (1972), J. Atm. Sci., 29, 959-970. Assumes 
# trans and sunalt both matrices of same size. Table 1 is called 
# albedot1.mat.
#
# INPUT:   trans - atmospheric transmittance 
#          sunalt - sun altitude [deg] 
# 
# OUTPUT:  alb - albedo 

########################################################################
# 3/10/96: version 1.0
# 7/24/98: version 1.1 (rev. to handle out-of-range input values by RP)
# 8/5/99:  version 2.0
########################################################################


# create axes
x<-seq(0, 90, 2)
y<-seq(0, 1, 0.05)
xx<-sort(rep(x, length(y)))
yy<-rep(y, length(x))
alb<-rep(NA, length(trans)) 
k<-which(sunalt>0 & is.finite(trans) & trans<=1.01)

# interpolate
k<-k[order(sunalt[k])]
dip<-diag(interp(xx,yy,albedot1,xo=sunalt[k],yo=trans[k])$z)
alb[sort(k)]<-dip[order(k)]
return(alb)
}