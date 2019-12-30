tridiag_DIF <- function(Kz, Vz, Az, dz, dt) {

# Returns a tridiagonal square matrix with vector b on diagonal
# and vectors a and c on sub- and super-diagonals, respectively
# Vectors a and c are assumed to be of length 1 less than length(b)
	
#	Nz<-length(Vz)  #number of grid points/layers
      n <- length(Vz)
	y <- array(0, c(n,n))

# Linearized heat conservation equation matrix (diffusion only)
   az <- (dt/dz) * Kz * (Az / Vz)                                        #coefficient for i-1
   cz <- (dt/dz) * c(Kz[-1], NA) * (c(Az[-1], NA) / Vz)            #coefficient for i+1
   bz <- 1 + az + cz                                                       #coefficient for i+1
#Boundary conditions, surface

  az[1] <- 0
  #cz(1) remains unchanged 
  bz[1]<- 1 + az[1] + cz[1]
   

#Boundary conditions, bottom

  #az(end) remains unchanged 
  cz[length(cz)] <- 0
  bz[length(bz)] <- 1 + az[length(az)] + cz[length(cz)]
	
	y[0 + 1:(n - 1) * (n + 1)] <- -cz[-length(bz)]	# superdiagonal
	y[1 + 0:(n - 1) * (n + 1)] <- bz	# diaygonal
	y[2 + 0:(n - 2) * (n + 1)] <- -az[-1] 	# subdiagonal
	return(y)
}
