cdntc<-function(sp,z,Ta){

#CHECK: simplified - works well for 1-dim input


tol<-.00001 # iteration endpoint

visc <- 1.326e-5*(1 + 6.542e-3*Ta + 8.301e-6*Ta^2 - 4.84e-9*Ta^3)

# remove any sp==0 to prevent division by zero
sp[which(sp==0)]<-.1

# initial guess
ustaro<-array(0,length(sp))
ustarn<-.036*sp
iter<-0
# iterate to find z0 and ustar
while(abs(ustarn-ustaro)>tol){
  ustaro<-ustarn
  z0<-Charnock_alpha*ustaro^2/g + R_roughness*visc/ustaro
  ustarn<-sp*(kappa/log(z/z0))
  iter<-iter+1
}

sqrcd<-kappa/log((10)/z0)
cd<-sqrcd^2

u10=ustarn/sqrcd
return(list(cd,u10))
}
