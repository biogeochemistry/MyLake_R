qsat<-function(Ta,Pa){

#if nargin==1,
#  source("Rfun/as_consts")
#  Pa=P_default # pressure in mb
#end;


ew = 6.1121*(1.0007 + 3.46e-6 * Pa) * exp((17.502 * Ta) / (240.97 + Ta)) # in mb
q  = 0.62197 * (ew / (Pa-0.378 * ew))						   #mb -> kg/kg
return(q)
}
