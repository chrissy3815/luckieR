#############################################################################
## R translations of formulas for mean/var of conditional lifespan, and 
## mean/var of conditional times to reach state i from state j, from 
##     Cochran, M. E. and S. Ellner.  1992.  Simple methods for 
##     calculating age-based life history parameters for stage-
##     structured populations. Ecological Monographs 62:345-364.
##
## Original coding by Erin Souder, merge into luckieR by Steve Ellner 
## NOTE: 																   
##     CE92 used the convention that all times (lifetime, transition time, etc.) 
##     are the largest possible value consistent with the census data. 
##     But that convention gives annuals a lifespan of 2 years, which was        
##     generally considered to be a bad idea, so we abandon it here. 

##     Your lifespan is the number of times that you were alive at a census 
##     time (rather than that number plus one). 
## 
##     For transition times, a state-j newborn *could have* lived up to one 
## 	   step prior to being censused in state j, so the CE92 formula gives 1 as 
##	   the conditional mean time to reach state j starting from state j. 
##     The code here gives instead 0: if you are first censuses in state j,  
## 	   we assume that it took you 'no time' to get there.
##############################################################################                                          
                  
##############################################################################

######################### NOTATION CONVENTIONS  ##############################
##      luckieR M = StageCoach P (state transitions) 
##      luckieR F = StageCoach B (sexual births) + F (clonal births, fission) 
############################################################################## 

### The defined functions are:
##			mean_lifespan, var_lifespan
##			Di_mat (a utility function)
## 			mean_conditional_time
##			mean_conditional_lifespan
##			var_conditional_time
##			var_conditional_lifespan (NOT FINISHED) 

#### For testing 
library(readxl); 
if( Sys.info()["user"]=="Ellner") setwd("c:/repos/luckieR/stagecoach"); 
M = as.matrix(read_xlsx("Caswell_P.xlsx")); 
F = as.matrix(read_xlsx("Caswell_B.xlsx")); 


M1 = matrix(c(0, 0,    0,  0,
              1, 0,    0,  0,
			  0, 0.5,  0,  0,
			  0, 0.5,  1,  0), 4,4,byrow=TRUE) 
N1 = solve(diag(4)-M1); colSums(N1); 			  
## 3.5 2.5 2.0 1.0  

M2 = matrix(c(0, 0,    0,  0,
              1, 0,    0,  0,
			  0, 0.25, 0,  0,
			  0, 0.25, 1,  0), 4,4,byrow=TRUE) 
N2 = solve(diag(4)-M2); colSums(N2); 			  
## 2.75 1.75 2.00 1.00
 
########################################################################### 
## Compute mean remaining lifespan given current state, for each state,       
##		given by row sums of the fundamental matrix. "Remaining lifespan"
##      includes the current year --- it is the number of years when the
##		individual is alive to be censused, including the current year and
##      all subsequent years, but not prior years (if there were any). 
##  Arguments:                                                            
##  		M: state transition matrix M                                  
##  Return Value/s:                                                       
##        Mean lifespan expectancy for each state i  					  
## 
##  Similar to life_expect_mean in Rage, but returns a vector of values
##  for each state, rather than one value based on a mixing distribution
##  across initial states 
########################################################################### 
mean_lifespan = function(M){ 
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}

  matDim = nrow(M); 
  N <- try(solve(diag(matDim) - M), silent = TRUE)
  if (inherits(N, "try-error")) {
        LE = rep(NA,matDim)
		} else {
        LE = colSums(N); 
  }
  return(LE)
}  
 
## mean_lifespan(M1);
## [1] 3.5 2.5 2.0 1.0
## Give the same results as life_expect_mean(M1,start=j) for j = 1,2,3,4

########################################################################### 
## Compute variance of remaining lifespan given current state. See 
##      mean_lifespan for the definition of "remaining lifespan". 
## Equation 5 in Cochran and Ellner 1992, which is equivalent 
##      to eqn. 5.12 in Caswell 2001 (proof: multiply each formula
##      for B on the right by N%*%N, and the result for both is (I+M).
##  Arguments:                                                            
##  		M: state transition matrix M                                  
##  Return Value/s:                                                       
##        Variance of remaining lifespan given current state i  					  
## 
##  Similar to life_expect_var in Rage, but returns a vector of values
##  for each state, rather than one value based on a mixing distribution
##  across initial states.  
########################################################################### 
var_lifespan = function(M){ 
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}
  
  eye = diag(nrow(M)); 
  N = solve(eye-M); 
  B = (eye + M)%*%N%*%N;
  results = colSums(B) - mean_lifespan(M)^2; 
  return(results)
}  

########################################################################
##  Function to calculate modified transition matrices D_i for each    
##      state i, where individuals entering state i stay there one     
##      time step and then die. This function is used by meantime()    
##      and would not be accessed directly by users.                   
##  Equation 8 in Cochran and Ellner 1992                              
##  Argument/s:                                                        
##     M: state transition matrix M                                    
##  Return Value/s:                                                    
##     Array of D of transition matrices, such that D[,,i] is D_i      
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
##  Function to compute the average time to reach state i from state j,    
##         conditional on that happening                                    
##  Equation 9 in Cochran and Ellner 1992                                  
##  Argument/s:                                                            
##     M: state transition matrix M                                        
##  Return Value/s:                                                        
##     Matrix, (i,j) entry = conditional mean time to reach i from j,      
##       with NA if the transition from j to i is not possible.                        
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

#########################################################################
##  Function to determine mean of lifespan conditional on reaching      
##      state i before death, for individuals born in state j.          
##  Equation 6 in Cochran and Ellner 1992                               
##  Arguments:                                                          
##  	M: state transition matrix M                                
##  Return Value/s:   
##     Matrix, (i,j) entry = mean lifespan conditional on reaching
##       state i before death, for an individual born in state j, 
##       with NA if the transition from j to i is not possible.                 
#########################################################################
mean_conditional_lifespan <- function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}
  
  LE <- mean_lifespan(M)
  MT <- mean_conditional_time(M)
  
  results <- matrix(0,nrow = nrow(M),ncol = ncol(M))
  for (i in 1:nrow(M)) {
	results[i,]=MT[i,] + LE[i]
  }
  return(results)
}  

######################################################################################
##  Function to determine the variance of the time to reach stage i from stage j           
##     conditional on that happening.                                                  
##  Equation 10 in Cochran and Ellner 1992                                           
##  Argument/s:                                                                      
##  	M: state transition matrix                                                             
##  Return Value/s:                                                                  
##      Matrix, (i,j) entry = variance of the time to reach i from j, conditional     
##       on that happening, with NA if the transition is not possible.  
##
##   Note: the difference in time conventions (shorter here by 1 than in CE92)
##     does not affect the variance, so here we use the CE92 convention and formula.
######################################################################################
var_conditional_time <- function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}
  
  D  = Di_mat(M)
  tau = mean_conditional_time(M) + 1 ## CE92 convention 
  
  results <- matrix(0, nrow(M), ncol(M))
  Imat <- diag(nrow(M))
  
  for (i in 1:nrow(M)) {
    Ni = solve(Imat - D[,,i])
    num <- (Imat + D[,,i]) %*% Ni %*% Ni %*% Ni 
    results[i,] <- (num[i,]/Ni[i,]) - tau[i,]^2; 
  }
  results[is.nan(results)] <- NA
  return(results)
  
}

  
#########################################################################
##  Function to calculate variance of lifespan conditional on reaching      
##      state i before death, for individuals born in state j.          
##  Equation 7 in Cochran and Ellner 1992                               
##  Arguments:                                                          
##  	M: state transition matrix M                                
##  Return Value/s:   
##     Matrix, (i,j) entry = variance of lifespan conditional on reaching
##       state i before death, for an individual born in state j, 
##       with NA if the transition from j to i is not possible.                 
#########################################################################
var_conditional_lifespan  = function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}
  
  Omega <- var_conditional_time(M)
  
  l <- life_expectancy_SD(P)
  results <- matrix(0,nrow = nrow(P), ncol = ncol(P))
  for (i in 1:nrow(P)) {
    
      results[i,] <- l[i] ^ 2 + m[i,] ^ 2  
    }
    
  return(sqrt(results))
  
}





