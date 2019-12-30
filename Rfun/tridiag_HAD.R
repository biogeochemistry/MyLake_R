tridiag_HAD_v11<-function(Kz,U,Vz,Az,dz,dt){

#if (U<0)
#    RPT obsolete: error('only positive (downward) velocities allowed')
#end

U[U==0] <- .Machine$double.eps #set Vz next to nothing (=2.2204e-016) in order to avoid division by zero

      #Nz=length(Vz); #number of grid points/layers
      n <- length(Vz)
	y <- array(0, c(n,n))
theta <- U*(dt/dz)
#RPT insert '.' for array instead of scalar multiplicattion
az <- theta*(1 + (1/(exp( (U*Vz)/(Kz*Az) ) - 1)))                  #coefficient for i-1
cz <- theta/(exp( (U*Vz)/(c(Kz[-1], NA)*c(Az[-1], NA)) ) - 1)  #coefficient for i+1
bz <- 1 + az + cz                                                       #coefficient for i

#Boundary conditions, surface

  az[1] <- 0
  #cz(1) remains unchanged
  bz[1] <- 1 + theta[1] + cz[1] # RPT changed into 'theta(1)'

#Boundary conditions, bottom

#az(end) remains unchanged
  cz[length(cz)] <- 0
  bz[length(bz)] <- 1 + az[length(az)]


	y[0 + 1:(n - 1) * (n + 1)] <- -cz[-length(bz)]	# superdiagonal
	y[1 + 0:(n - 1) * (n + 1)] <- bz	# diaygonal
	y[2 + 0:(n - 2) * (n + 1)] <- -az[-1] 	# subdiagonal
	return(y)
}