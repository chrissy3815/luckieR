###########################################################################
## R translation of formulas for mean/sd of conditional lifespan,from
##     Cochran, M. E. and S. Ellner.  1992.  Simple methods for 
##     calculating age-based life history parameters for stage-
##     structured populations. Ecological Monographs 62:345-364.
##
## Original coding by Erin Souder, merge into luckieR by Steve Ellner 
###########################################################################

######################### NOTATION CONVENTIONS ############################
##      luckieR M = StageCoach P (state transitions) 
##      luckieR F = StageCoach B (sexual births) + F (clonal births, fission) 
########################################################################### 

#### For testing 
library(readxl); 
if( Sys.info()["user"]=="Ellner") setwd("c:/repos/luckieR/stagecoach"); 
M = as.matrix(read_xlsx("Caswell_P.xlsx")); 
F = as.matrix(read_xlsx("Caswell_B.xlsx")); 


M1 = matrix(c(0, 0,  0,  0,
              1, 0,  0,  0,
			  0, 0.5, 0, 0,
			  0, 0.5, 1, 0), 4,4,byrow=TRUE) 
N1 = solve(diag(4)-M1); colSums(N1); 			  
## 3.5 2.5 2.0 1.0   

########################################################################### 
## Compute mean life expectancy given current state, for all states,     ##
##		given by row sums of the fundamental matrix                      ##
##  Arguments:                                                           ##
##  		M: state transition matrix M                                 ##
##  Return Value/s:                                                      ##
##        Mean lifespan expectancy for each state i  					 ##
########################################################################### 
life_expectancy = function(M){ 
  matDim = nrow(M); 
  N <- try(solve(diag(matDim) - M), silent = TRUE)
  if (inherits(N, "try-error")) {
        LE = rep(NA,matDim)
		} else {
        LE = colSums(N); 
  }
  return(LE)
}  

########################################################################
##  Function to calculate modified transition matrix D_i for state i, ##
##      where individuals entering state i stay there for one time    ##
##      step and then die. This is a utility for meantime() and       ##
##      generally not of value to users.                              ##
##  Equation 8 in Cochran and Ellner 1992                             ##
##  Argument/s:                                                       ##
##     M: state transition matrix M                                   ##
##  Return Value/s:                                                   ##
##     Array of D of transition matrices, such that D[,,i] is D_i     ##
########################################################################
Di_mat <- function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}
  
  results <- array(0, dim = c(dim(M),ncol(M)))
  for (i  in 1:ncol(M)) {
      results[,,i] = M; 
      results[,i,i] = 0; 
  }
  return(results) 
}

############################################################################
##  Function to compute the average time to reach state i from state j,   ##
##         conditional on that happening                                  ## 
##  Equation 9 in Cochran and Ellner 1992                                 ##
##  Argument/s:                                                           ##
##     M: state transition matrix M                                       ##
##  Return Value/s:                                                       ##
##     Matrix, (i,j) entry = conditional mean time to reach i from j,     ##
##       with NA if the transition is not possible.                       ##
##  NOTE: 																  ##
##     We subtract 1 from the CE92 formula. In CE92, all times are the    ##
##	   largest possible value consistent with the data, in particular: a  ##
##	   state-j newborn *could have* taken one time step to get where it   ##
##     now is, so the CE92 formula gives 1 on the diagonal of the result. ##
##     Here we put zeros on the diagonal: if first seen in state j, it    ##
##     too you 'no time' to get there.                                    ##
##     The CE92 convention gives annuals a lifespan of 2 years. That      ##
##     was generally considered to be a bad idea.                         ##
############################################################################
mean_conditional_time <- function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}
  
  D <- Di_mat(M)
  
  results <- matrix(0, nrow(M), ncol(M))
  for (i in 1:nrow(M)) {
    Imat <- diag(x=1, dim(D[,,i])[1])
    den <- solve(Imat - D[,,i])  
    num <- den %*% den
    results[i,] <- (num[i,]/den[i,])-1
  }

  results[is.nan(results)] <- NA
  # NA used when because can't get to this state from previous state
  
  return(results)
}

SD_conditional_time <- function(M){
    
  ######################################################################################
  ##  Function to determine the SD of the time to reach stage i from stage j          ##
  ##     conditional on that happening.                                               ##  
  ##  Equation 10 in Cochran and Ellner 1992                                          ##
  ##  Argument/s:                                                                     ##
  ##  P: Survival matrix P                                                            ##
  ##  Return Value/s:                                                                 ##
  ##  The standard deviation for the time to reach stage i from stage j,              ##
  ##           conditional on that happening 
  ##  Author/s:                                                                       ##
  ##  Dr. Stephen Ellner                                                              ##
  ##  Dr. Simone Blomberg                                                             ##
  ##  Erin Souder                                                                     ##
  ##  Date: 03/01/2023                                                                ##
  ######################################################################################
  
  if(!is.matrix(P)) {stop("This is not a matrix")} 
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  
  #OUTPUT DOESN'T MATCH STAGECOACH
  D <- Di_mat(P)
  m0 <- meantime(P)
  m1 <- m0 + 1
  # m1 used to match equations from paper
  results <- matrix(0, nrow = nrow(P), ncol= ncol(P))
  # Make new matrix 
  
  for (i in 1:dim(results)[1]) {
    Imat <- diag(x=1, dim(D[,,i])[1])
    num <- (Imat + D[,,i]) %*% solve((Imat - D[,,i]) %*% (Imat - D[,,i]) %*% (Imat - D[,,i])) 
    den <- solve(Imat - D[,,i]) 
    results[i,] <- (num[i,]/den[i,])
    results[i,] <- sqrt(results[i,] - (m1[i,]) ^ 2)
  }
  results[is.nan(results)] <- NA
  results[m0 < 0] <- NA
  # Transition not defined for this so it's NA
 return(results)
  
}
















#########################################################################
##  Function to determine mean of lifespan conditional on reaching     ##
##      state i before death, for individuals born in state j.         ##
##  Equation 6 in Cochran and Ellner 1992                              ##
##  Arguments:                                                         ##
##  		M: state transition matrix M                               ##
##  Return Value/s:                                                    ##
##        Conditional mean lifespan for each state i                   ##
#########################################################################
conditional_mean_lifeSpan <- function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}
  
  LE <- life_expectancy(M)
  MT <- meantime(M)
  
  results <- matrix(0,nrow = nrow(M),ncol = ncol(M))
  for (i in 1:nrow(M)) {
	results[i,]=MT[i,] + LE[i]
  }
  return(results)
}  

