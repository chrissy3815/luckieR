## Convert state matrix to state vector
vec <- function(nmat) matrix(nmat,ncol=1) 

## Convert state vector to state matrix
unvec <- function(nvec,nrow=NULL,ncol=NULL){ 
	if(is.null(nrow)) return(matrix(nvec,ncol=ncol)); 
	if(is.null(ncol)) return(matrix(nvec,nrow=nrow)); 
}

## Convert iteration array to iteration matrix 
flatten <- function(A4) {
	dim(A4)<-rep(sqrt(length(A4)),2)
	return(A4)
}

## Convert iteration matrix to iteration array
unfold <- function(A2,dim) {
    dim(A2) <- dim; 
    return(A2)
} 

######################################################################
#  Function to make the B matrix B[i,j] is the probability that a
#  class-j individual has i-1 kids.  We assume Poisson-distributed
#  number of offspring.
#
#  maxKids = maximum number of kids = 20 by default
#############################################################

mk_B = function (maxKids=20, F, Fdist) {
  bigmz = ncol(F)
  B = matrix (0, maxKids+1, bigmz)
  if (Fdist == "Poisson") {
    for (z in 1:bigmz)
      B[,z] = dpois (0:maxKids, lambda=sum(F[,z]))
  } else if (Fdist == "Bernoulli") {
    for (z in 1:bigmz)
      B[,z] = dbinom (0:1, size=1, prob=sum(F[,z]))
  } else {
    stop ("mk_B: I don't recognize that option for Fdist.\n")
  }

  return (B)
}

######################################################################
# Function to take B and M matrices, and compute the transition 
# probabilities from size-class i and j total kids, to all size classes
# and l total kids. This returns a vector of zeros if (l-j) is < 0
# or above the assumed maxnumber of kids per year,  ncol(B) - 1. 
######################################################################
p_xT <- function(l, i, j, B, M) {
  bigmz <- ncol(M); maxKids <- nrow(B)-1;
  newKids <- (l-j); 
  if((newKids < 0) | (newKids > maxKids)) {
    return(rep(0, bigmz))
  }else{
    return(M[,i]*B[newKids+1,i])
  }
}

##############################################################################
## Function to make the 2D iteration matrix A for a size-kids model
## based on the M and B matrices summarizing a size-structured
## IPM. Apart from B and M the only input is mT, dimension for T 
## (so range of T is 0 to mT-1). The 4-D iteration array K is also returned. 
##
## Iteration matrix is modifited so individuals who get to the maximum
## values of T in the matrix stay there, but continue to grow/shrink/die
#############################################################################
make_AxT <- function(B, M, mT) {
  bigmz=ncol(M); Kvals=array(0,c(bigmz,mT,bigmz,mT));  
  for(z in 1:bigmz){ # initial size
    for(k in 1:mT){ # initial T 
      for(kp in 1:(mT-1)){ # final T
        Kvals[,kp,z,k]=p_xT(kp,z,k,B,M)
      }
      for (zp in 1:bigmz)
        ## make last kids class absorbing
        Kvals[zp,mT,z,k] = M[zp,z] - sum(Kvals[zp,,z,k])
    }
  }
  ## make kids-class mT absorbing: stay there with prob=1
  Kvals[1:bigmz,1:mT,1:bigmz,mT] <- 0; 
  Kvals[1:bigmz,mT,1:bigmz,mT] <- M; 
  A <- Kvals; dim(A) <- c(bigmz*mT,bigmz*mT); 

  return(list(A=A,K=Kvals)) 
}  

############################################################################
## Function to make megamatrix from a list of "fast" matrices
## (e.g. list of survival/growth transition matrices or fecundity
## matrices for different environments) and a "slow" matrix, e.g. the
## environment transition matrix.  The cross-classified states cycle
## through all the states in the relevant fast matrix for the first
## state in slowMatrix, then again for the second state in slowMatrix,
## etc.
##
## fastMatrixList: a list of matrices.  There should be a fast matrix
## for every state in slow matrix.  E.g. a fecundity matrix for each
## environment.
##
## slowMatrix: a matrix
############################################################################

makeM = function (fastMatrixList, slowMatrix) {  
  dimSlow = dim(slowMatrix)[1]
  dimFast = dim(fastMatrixList[[1]])[1]

  if (length(fastMatrixList) != dimSlow) 
    stop("makeM: There must be a fast matrix for each state in the slow matrix.\n")
    if (dim(slowMatrix)[2] != dimSlow)
    stop ("makeM: Matrices must be square.\n")
  if (dim(fastMatrixList[[1]])[2] != dimFast)
    stop ("makeM: Matrices must be square.\n")

  M = matrix (NA, dimFast*dimSlow, dimFast*dimSlow)
  
  for (i in 1:dimSlow) {
    for (j in 1:dimSlow) {
      M[(i-1)*dimFast + 1:dimFast, (j-1)*dimFast + 1:dimFast] =
        fastMatrixList[[j]]*slowMatrix[i,j]
    }
  }

  return (M)
}

############################################################################
## Function to make an array of megamatrices, one for each trait
## value, from a list of "fast" matrices (e.g. list of survival/growth
## transition matrices or fecundity matrices for different
## environments) and a "slow" matrix, e.g. the environment transition
## matrix.  The cross-classified states cycle through all the states
## in the relevant fast matrix for the first state in slowMatrix, then
## again for the second state in slowMatrix, etc.
##
## fastMatrixArray: an array of matrices.  The first index indicates
## trait, the second index indicates environment (or whatever the slow
## state is), the last two indices are for the fast transition
## matrices associated with each trait and environment (or whatever
## the slow state is).
############################################################################

makeMArray = function (fastMatrixArray, slowMatrix) {
  numTraits = dim(fastMatrixArray)[1]
  dimSlow = dim(fastMatrixArray)[2]
  dimFast = dim(fastMatrixArray)[3]

  if (dim(slowMatrix)[1] != dimSlow) 
    stop("makeM: There must be a fast matrix for each state in the slow matrix.\n")
  if (dim(slowMatrix)[2] != dimSlow)
    stop ("makeM: Matrices must be square.\n")
  if (dim(fastMatrixArray)[4] != dimFast)
    stop ("makeM: Matrices must be square.\n")

  Marray = array (NA, dim=c(numTraits, dimFast*dimSlow, dimFast*dimSlow))
      
  for (x in 1:numTraits) {
    for (i in 1:dimSlow) {
      for (j in 1:dimSlow) {
        Marray[x, (i-1)*dimFast + 1:dimFast,
        (j-1)*dimFast + 1:dimFast] =
          fastMatrixArray[x,j,,]*slowMatrix[i,j]
      }
    }
  }
  
  return (Marray)
}
