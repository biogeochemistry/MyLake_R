# RPT: removed nargin>8 options
# As hfbulktc.m but in order to speed up the function, the iteration loop is 
# changed from "if" to "while" structure by TSA 05.09.03

hfbulktc_speed<-function(ur,zr,Ta,zt,rh,zq,Pa,Ts){
# HFBULKTC: computes sensible and latent heat fluxes and other variables.
# A=HFBULKTC(ur,zr,Ta,zt,rh,zq,Pa,Ts,sal,dlw,dsw,nsw) computes the following:
#
#         Hs      = sensible heat flux INTO ocean [W/m^2]
#         Hl      = latent heat flux INTO ocean [W/m^2]
#         Hl_webb = Webb correction to latent heat flux INTO ocean [W/m^2]
#         stress  = wind stress [N/m^2]
#         U_star  = velocity friction scale [m/s]
#         T_star  = temperature scale [deg C]
#         Q_star  = humidity scale [kg/kg]
#         L       = Monin-Obukhov length [m]
#         zetu    = zr/L
#         CD      = drag coefficient
#         CT      = temperature transfer coefficient (Stanton number)
#         CQ      = moisture transfer coefficient (Dalton number)
#         RI      = bulk Richardson number
#         Dter    = cool-skin temperature difference (optional output) [C]; 
#                   positive if surface is cooler than bulk (presently no 
#                   warm skin permitted by model)
#                    
# Based on the following buoy input data:
#
#           ur     = wind speed [m/s] measured at height zr [m] 
#           Ta     = air temperature [C] measured at height zt [m]
#           rh     = relative humidity [#] measured at height zq [m]
#           Pa     = air pressure [mb]
#           Ts     = sea surface temperature [C]
#
# where ur, Ta, rh, Pa, Ts, zr, zt, and zq may be either row or column vectors; and rh, Pa, 
# zr, zt, and zq may also be fixed scalars.
#
# Output variables are given as column vectors in A:
#
# 1) without cool-skin correction:
#
#   A=[Hs Hl Hl_webb stress U_star T_star Q_star L zetu CD CT CQ RI]
#
# 2) with cool-skin correction: 
#
#   A=[Hs Hl Hl_webb stress U_star T_star Q_star L zetu CD CT CQ RI Dter];

# Code follows Edson and Fairall TOGA COARE code (version 2.0), modified 
# to include Rogers' weighting factor for unstable conditions.  Code does
# include gustiness, and assumes that the marine boundary layer height is
# known and constant over time for simiplicity. zr/L is limited to 
# be <=3.0 to ensure that the code converges to nonzero stress and heat 
# flux values for strongly stable conditions.  The bulk Richardson number
# is computed between the sea surface and zr as a diagnostic about whether
# turbulent boundary layer theory is applicable.  Code does not include 
# warm layer effects to modify Ts.  See Fairall et al (1996), J. Geophys. 
# Res., 101, 3747-3764, for description of full TOGA COARE code and 
# comparison with data. 

######################################################################
# 8/19/98: version 1.1 (rewritten by RP to remove inconsistencies in 
#          virtual and real temperatures, improve loop structure, 
#          correct gustiness component of stress computation) 
# 4/9/99: version 1.2 (rewritten by AA to clarify some variable names
#         and include cool-skin effect and Webb correction to latent 
#         heat flux added to output matrix)
# 8/5/99: version 2.0
######################################################################


M <- length(ur)

# [] RPT excluded


# initialize various constants


tol<-.001;    # tolerance on Re changes to make sure soln has converged.

onethird <- 1/3;
o61<-1/eps_air-1;   # 0.61 (moisture correction for temperature)
visc <- 1.326e-5*(1 + 6.542e-3*Ta + 8.301e-6*Ta^2 - 4.84e-9*Ta^3) # viscosity

Qsats<-Qsat_coeff*qsat(Ts,Pa)      # saturation specific humidity; the Qsat_coeff
                                  # value is set in routine as_consts.m
Q <- (0.01 *rh) *qsat(Ta,Pa)      # specific humidity of air [kg/kg]
T <- Ta+CtoK                      # convert to K
Tv<-T *(1 + o61*Q);               # air virtual temperature
rho<-(100*Pa) /(gas_const_R*Tv)   # air density
Dt<-(Ta+0.0098 *zt)-Ts            # adiabatic temperature difference
Dq<-Q-Qsats                       # humidity difference

# compute initial neutral scaling coefficients
S<-sqrt(ur^2 + min_gustiness^2)
cdnhf<-sqrt((cdntc(S,zr,Ta))[[1]]) # Smith's neutral cd as first guess
z0t<-7.5*10^(-5)
ctnhf<-kappa/log(zt/z0t)

z0q<-z0t
cqnhf<-kappa/log(zq/z0q)

U_star <- cdnhf*S      # (includes gustiness)
T_star <- ctnhf*Dt     # 
Q_star <- cqnhf*Dq     #


Dter   <- 0
Dqer   <- 0
# pause
Reu=0;Ret=0;Req=0;
# begin iteration loop to compute best U_star, T_star, and Q_star
                #for iter1=1:80;
cvrgnce<-1;      #modified by TSA
ii <- 0;         #modified by TSA
while (cvrgnce){ #modified by TSA
    ReuO<-Reu; RetO<-Ret; ReqO<-Req; # Save old values
    
    # Compute Monin-Obukov length (NB - definition given as eqn (7)
    # of Fairall et al (1996) probably wrong, following, e.g.
    # Godfrey and Bellars (1991), JGR, 96, 22043-22048 and original code)
    bs<-g*(T_star*(1 + o61*Q) + o61*T*Q_star)/Tv 
    L<-(U_star^2)/(kappa*bs)
    # set upper limit on zr/L = 3.0 to force convergence under 
    # very stable conditions. Assume that zr, zt and zq comparable.
    index_limit  <- (L<zr/3 & L>0)
    L[index_limit]=zr[index_limit]/3
    
    zetu<-zr/L;  # nondimensionalized heights
    zett<-zt/L;
    zetq<-zq/L;

    # surface roughness
    z0<-(Charnock_alpha/g)*U_star^2 + R_roughness*visc/U_star

    # compute U_star correction for non-neutral conditions
    cdnhf<-kappa/(log(zr/z0)-psi(zetu, "utc"))
    U_star<-cdnhf*S
  
    Reu<-z0*U_star/visc;   # roughness Reynolds #
    Ret<-LKB(Reu)[[1]];  # compute other roughness Reynolds #s
    Req<-LKB(Reu)[[2]];

    # compute t and q roughness scales from roughness R#s
    z0t<-visc*Ret/U_star
    z0q<-visc*Req/U_star

    # compute new transfer coefficients at measurement heights
    cthf<-kappa/(log(zt/z0t)-psi(zett, "ttc"))
    cqhf<-kappa/(log(zq/z0q)-psi(zetq, "ttc"))

    # compute new values of T_star, Q_star
    T_star<-cthf*(Dt + Dter)
    Q_star<-cqhf*(Dq + Dqer)

    # estimate new gustiness
    Ws<-U_star*(-CVB_depth/(kappa*L))^onethird
    wg<-rep(min_gustiness, length=M)
    j<-which(zetu<0)                 # convection in unstable conditions only
    wg[j]<-max(min_gustiness,beta_conv*Ws[j]) # set minimum gustiness
    S<-sqrt(ur^2 + wg^2)

    # [] excluded coolskin loop
    cvrgnce <- abs(Reu-ReuO)>tol | abs(Ret-RetO)>tol | abs(Req-ReqO)>tol; #modified by TSA

    ii = ii + 1
    if (ii>80){warning("Iteration did not converge (hfbulktc_speed)!"); break}
} # end of iteration loop

# compute latent heat
Le <-(2.501-0.00237*(Ts-Dter))*10^6

# compute fluxes into ocean
Hs <- rho*cp*U_star*T_star
Hl <- rho*Le*U_star*Q_star

# compute transfer coefficients at measurement heights
CD <- (U_star/S)^2
CT <- U_star*T_star/(S*(Dt + Dter)) # Stanton number
CQ <- U_star*Q_star/(S*(Dq + Dqer)) # Dalton number

# to compute mean stress, we don't want to include the effects
# of gustiness which average out (in a vector sense).
stress <- rho*CD*S*ur

# compute bulk Richardson number (as a diagnostic) - the "T"
# is probably not quite right - assumes T \ approx. Ts (good enough though)
RI <- g*zr*((Dt + Dter) + o61*T*(Dq + Dqer))/(Tv*S^2)

# compute Webb correction to latent heat flux into ocean
W <- 1.61*U_star*Q_star + (1 + 1.61*Q)*U_star*T_star/T # eqn. 21
Hl_webb <- rho*Le*W*Q # eqn. 22, Fairall et al. (1996), JGR, 101, p3751.

return(list(Hs, Hl, Hl_webb, stress, U_star, T_star, Q_star, L, zetu, CD, CT, CQ, RI))
}
##### end hfbulktc_speed


