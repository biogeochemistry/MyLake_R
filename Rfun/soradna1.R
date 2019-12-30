soradna1<-function(dat, res,long,lat){
# RPT: function works only for one spot: long, lat must be scalars
# SORADNA1: computes no-sky solar radiation and solar altitude.
# [z,sorad]=SORADNA1(yd,yr,long,lat) computes instantaneous values of 
# solar radiation and solar altitude from yearday, year, and position 
# data. It is put together from expressions taken from Appendix E in the
# 1978 edition of Almanac for Computers, Nautical Almanac Office, U.S.
# Naval Observatory. They are reduced accuracy expressions valid for the
# years 1800-2100. Solar declination computed from these expressions is
# accurate to at least 1'. The solar constant (1368.0 W/m^2) represents a 
# mean of satellite measurements made over the last sunspot cycle (1979-1995) 
# taken from Coffey et al (1995), Earth System Monitor, 6, 6-10.  Assumes 
# yd is either a column or row vector, the other input variables are scalars,
# OR yd is a scalar, the other inputs matrices.
#
#  INPUT:   dt   - date in chron format
#           res  - resolution(steps per day)
#           long - longitude (west is positive!) [deg] 
#           lat  - latitude  [deg]
#old INPUT:   yd   - decimal yearday (e.g., 0000Z Jan 1 is 0.0)
#old          yr   - year (e.g., 1995)
#
#  OUTPUT:  z    - solar altitude [deg]
#           sorad- no atmosphere solar radiation  [W/m^2]

######################################################################
# 3/8/97: version 1.0
# 8/28/98: version 1.1 (vectorized by RP)
# 8/5/99: version 2.0
#######################################################################


dat<-dat+seq(0,1,1/res)[-(res+1)]# make time series for focal day



SC<-seconds(dat)
MN<-minutes(dat)
HR<-hours(dat)
D<-M<-Y<-rep(NA, length=length(MN))
D[]<-as.numeric(unlist(strsplit(as.character(dates(dat)[1]), "/"))[2])
M[]<-as.numeric(unlist(strsplit(as.character(dates(dat)[1]), "/"))[1])
Y[]<-as.numeric(unlist(strsplit(as.character(dates(dat)[1]), "/"))[3])
if(Y[1]<1900) Y[]<-Y[1]+2000

# convert to new variables
LONG<-rep(long, length=length(SC))
LAT<-rep(lat, length=length(SC))


# constants
   DTR<-3.14159265/180
   RTD<-1/DTR

# compute Universal Time in hours
   UT<-HR+(MN+SC/60)/60

# compute Julian ephemeris date in days (Day 1 is 1 Jan 4713 B.C.=-4712 Jan 1)
  JED<-367*Y-floor(7*(Y+floor((M+9)/12))/4)+floor(275*M/9)+D+1721013 + UT/24

# compute interval in Julian centuries since 1900
    T<-(JED-2415020.0)/36525

# compute mean anomaly of the sun
    G<-358.475833+35999.049750*T-.000150*T^2
   NG<-floor(G/360)
    G<-(G-NG*360)*DTR

# compute mean longitude of sun
    L<-279.696678+36000.768920*T+.000303*T^2
   NL<-floor(L/360)
    L<-(L-NL*360)*DTR

# compute mean anomaly of Jupiter
  JUP<-225.444651+2880.0*T+154.906654*T
 NJUP<-floor(JUP/360)
  JUP<-(JUP-NJUP*360)*DTR

# compute longitude of the ascending node of the moon's orbit
   NM<-259.183275-1800*T-134.142008*T+.002078*T^2
  NNM<-floor(NM/360)
   NM<-(NM-NNM*360+360)*DTR

# compute mean anomaly of Venus
    V<-212.603219+58320*T+197.803875*T+.001286*T^2
   NV<-floor(V/360)
    V<-(V-NV*360)*DTR

# compute sun theta
 THETA<-.397930*sin(L)+.009999*sin(G-L)+.003334*sin(G+L)-
     .000208*T*sin(L)+.000042*sin(2*G+L)-.000040*cos(L)-
     .000039*sin(NM-L)-.000030*T*sin(G-L)-.000014*sin(2*G-L)-
     .000010*cos(G-L-JUP)-.000010*T*sin(G+L)

# compute sun rho
  RHO<-1.000421-.033503*cos(G)-.000140*cos(2*G)+
     .000084*T*cos(G)-.000033*sin(G-JUP)+.000027*sin(2*G-2*V)

# compute declination
   DECL<-asin(THETA/sqrt(RHO))

# compute equation of time (in seconds of time) (L in degrees)
    L <-276.697+0.98564734*(JED-2415020.0)
    L <-(L - 360*floor(L/360))*DTR
  EQT <--97.8*sin(L)-431.3*cos(L)+596.6*sin(2*L)-1.9*cos(2*L)+
         4.0*sin(3*L)+19.3*cos(3*L)-12.7*sin(4*L)
  EQT <- EQT/60
    L <- L*RTD

# compute local hour angle
  GHA <- 15*(UT-12) + 15*EQT/60
  LHA <- GHA - LONG

# compute radius vector
   RV<-sqrt(RHO)

# compute solar altitude
   SZ<-sin(DTR*LAT)*sin(DECL)+cos(DTR*LAT)*cos(DECL)*cos(DTR*LHA)
    z<-RTD*asin(SZ)

# compute solar radiation outside atmosphere
sorad<-rep(0, length=length(z))
ii<-which(z>0)
sorad[ii]<-(Solar_const/RV[ii]^2)*sin(DTR*z[ii])
return(list(z, sorad))
}
