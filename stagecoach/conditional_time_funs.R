#############################################################################
## R translations of formulas for mean/var of conditional lifespan, and 
## mean/var of conditional times to reach state i from state j, from 
##     Cochran, M. E. and S. Ellner. 1992.  Simple methods for 
##     calculating age-based life history parameters for stage-
##     structured populations. Ecological Monographs 62:345-364.
##
## Original coding by Erin Souder, revised and extended for luckieR 
## by Steve Ellner 
##
## NOTE on time conventions: 																   
##     CE92 used the convention that all times (lifetime, transition time, etc.) 
##     are the supremum of values that are consistent with the census data. 
##     See figure 2 of CE92. That convention gives annuals a lifespan of 
## 	   2 years, which was generally considered to be a bad idea, 
##     so we abandon it here. 
##
##     An individual's lifespan is the number of times that they were alive  
##     at a census time (rather than that number plus one). 
## 
##     For transition times, a state-j newborn *could have* lived up to one 
## 	   step prior to being censused in state j, so the CE92 formula gives 1 as 
##	   the conditional mean time to reach state j starting from state j. 
##     The code here gives instead 0: if you are first censused in state j,  
## 	   we assume that it took you 'no time' to get there.
##############################################################################                                          
                  
##############################################################################

######################### NOTATION CONVENTIONS  ##############################
##      luckieR M = StageCoach P (state transitions) 
##      luckieR F = StageCoach B (sexual births) + F (clonal births, fission) 
############################################################################## 

### The defined functions are:
##			mean_lifespan, var_lifespan
## 			mean_conditional_times, var_conditional_times 
##			mean_conditional_lifespan, var_conditional_lifespan 
##			first_breed -- prob. of breeding before death, mean & var of age 
##                         at first breeding, conditional on breeding 
##			wrappers for population moments: pop_mean_var, pop_mu3, pop_skew 
##			Di_mat (a utility function)

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
##  across initial states. 
########################################################################### 
mean_lifespan = function(M){ 
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}

  matDim = nrow(M); 
  N = try(solve(diag(matDim) - M), silent = TRUE)
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
mean_conditional_times <- function(M){
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
  MT <- mean_conditional_times(M)
  
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
##     does not affect the variance, so we use the CE92 convention and formula.
######################################################################################
var_conditional_times <- function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}
  
  D  = Di_mat(M)
  tau = mean_conditional_times(M) + 1 ## use CE92 convention 
  
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
##
#########################################################################
var_conditional_lifespan = function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}
  
  # variance of conditional time to reach i from j 
  var_tau_ij = var_conditional_times(M)   
  
  # variance of remaining lifespan starting at i 
  var_Omega_i = var_lifespan(M) 
  
  nx = nrow(M); results = matrix(NA,nx,nx);  
  for(i in 1:nx){
  for(j in 1:nx){
	results[i,j] = var_tau_ij[i,j] + var_Omega_i
  }}   
  return(results) 
}

#########################################################################
##  Function to compute probability of breeding at least once, and 
##  mean and variance of age at first breeding conditional on breeding,
##  as a function of initial state. This requires that the model 
##  specification includes state-dependent breeding probability, not just
##  state-dependent mean fecundity. 
##
##  Probability of breeding at least once from now until death must be 
##  positive for all states. If not, remove those states from the model
##  (can't start in any of them, going into them becomes death) and 
##  apply this function to the remaining states. 
##
##  Assumptions and coding follow Ellner et al. (2016) IPM monograph,
##  chapter 3. The state transition matrix M is assumed to have the
##  form M = p_b M_b + (1-p_b) M_0 to allow for costs of reproduction,
##  but M_b and M_0 can be equal if there are no costs. 
##  Breeding does not necessarily imply producing any new recruits,
##  for example clutch size conditional on breeding could be Poisson. 
##
##  This function complements mature_age() in Rage, which computes
##  the mean age at first entering a state with positive mean fecundity. 
##
##  Arguments:                                                          
##  	M_0: transition matrix for non-breeders
##		p_b: probability breeding for each state.                                 
##  Return Value/s:   
##	   B: prob. to breed at once, condional on initial state	(vector)
##     abar_R: mean age at first breeding, conditional initial state (vector)        
#########################################################################
first_breed = function(M0,pb) {
  if(!is.matrix(M0)) {stop("This is not a matrix")}
  if(nrow(M0)!=ncol(M0)) {stop("This is not a square matrix")}
  if(!is.numeric(M0)) {stop("This is not a numeric matrix")}
  if(length(pb)!=nrow(M0)) {stop("Length mismatch, breeding probability")} 

    nx = nrow(M0)
	
	## survival and growth without breeding
	P0 = matrix(NA,nx,nx)
	for(i in 1:nx) P0[i,] = (1-pb)*M0[i,]
	N0 = solve(diag(nx)-P0); 
	
	## IPM book, page 70 
	B = as.numeric(matrix(pb,1,nx)%*%M0) 
	
	## IPM book, page 71 
	Pb = matrix(NA,nx,nx) 
	for(z in 1:nx) Pb[,z] = P0[,z]*B/B[z] 
	Nb = solve(diag(nx)-Pb); 
	abar_R = colSums(Nb)-1; 
	# -1 because age at birth = 0 by assumption. A 'lifespan' of 2 in Pb 
	# means that you breed at your second census time, which is age=1. 
	
	## Caswell formula, only involves fundamental matrix
	var_R = colSums(2*(Nb%*%Nb) - Nb) - colSums(Nb)^2;  
	
	return(list(p_breed=B,mean_age = abar_R, var_age = var_R))
}	
	
#########################################################################
##  Function to compute the population mean and variance of some 
##  attribute X, given the vectors of mean and variance conditional on 
##  individual 'type' Z. This substitutes for all of the individual  
##  formulas in Table 2 of Cochran and Ellner 1992. 
##  Arguments:                                                          
##  	mean_by_type: mean of X for each value of Z (vector)
##		var_by_type: variance of X for each value of Z (vector)                                 
##		mixdist: frequency distribution of Z within the population. 
##  Return Value/s:   
##     population_mean, population_variance: scalars  
##            
#########################################################################
pop_mean_var = function(mean_by_type,var_by_type,mixdist){
		if(min(mixdist)<0) {stop("Not a valid mixing distribution")} 
		if(sum(mixdist)!=1) {stop("Not a valid mixing distribution")}
		if(min(var_by_type)<0) {stop("Not a valid variances vector")}
		nZ = length(mean_by_type) 
		if(length(var_by_type)!=nZ) {stop("Length mismatch, mean and var")}
		if(length(mixdist)!=nZ) {stop("Length mismatch, mean and mixing distribution")}
		
		pop_mean_X = sum(mixdist*mean_by_type); 
		EXsq_by_type = var_by_type + mean_by_type^2 
		pop_mean_Xsq = sum(mixdist*EXsq_by_type);
		pop_var_X = pop_mean_Xsq - (pop_mean_X)^2; 
		
		return(list(pop_mean = pop_mean_X, pop_var = pop_var_X))
}		

#########################################################################
##  Function to compute the population third central moment of some 
##  attribute X, given the vectors of mean, variance and mu3 conditional  
##  individual 'type' Z. This is done using the order-3 case of the  
##  Law of Total Cumulance 
##  Arguments:                                                          
##  	mean_by_type: mean of X for each value of Z (vector) 
##		var_by_type: variance of X for each value of Z (vector)                                 
##      mu3_by_type: 3rd central moment of X for each value of Z (vector)
##		mixdist: frequency distribution of Z within the population. 
##  Return Value/s:   
##     population_mu3 (scalar)  
##            
#########################################################################
pop_mu3 = function(mean_by_type,var_by_type,mu3_by_type,mixdist){
		if(min(mixdist)<0) {stop("Not a valid mixing distribution")} 
		if(sum(mixdist)!=1) {stop("Not a valid mixing distribution")}
		if(min(var_by_type)<0) {stop("Not a valid variances vector")}
		nZ = length(mean_by_type) 
		if(length(var_by_type)!=nZ) {stop("Length mismatch, mean and var")}
		if(length(mu3_by_type)!=nZ) {stop("Length mismatch, mean and mu3")}
		if(length(mixdist)!=nZ) {stop("Length mismatch, means and mixing distribution")}
		
		term1 = sum(mixdist*mu3_by_type)    # E[mu3(X|Z)] 
		
		mean_mean_by_type = sum(mixdist*mean_by_type); 
		term2 = sum(mixdist* ((mean_by_type - mean_mean_by_type)^3) );  
		
		
		Emean_times_var = sum(mixdist*(mean_by_type*var_by_type)) 
		Emean_times_Evar = sum(mixdist*mean_by_type)*sum(mixdist*var_by_type); 
		term3 = 3*(Emean_times_var - Emean_times_Evar) 

		result = term1+term2+term3
		return(result)
}		

#########################################################################
##  Function to compute the population skewness of some attribute
##  attribute X, given the vectors of mean, variance and skemness   
##  conditional on individual 'type' Z. Done by calling pop_mu3
##  and pop_mean_var.  
##  Arguments:                                                          
##  	mean_by_type: mean of X for each value of Z (vector) 
##		var_by_type: variance of X for each value of Z (vector)                                 
##      skew_by_type: skewness of X for each value of Z (vector)
##		mixdist: frequency distribution of Z within the population. 
##  Return Value/s:   
##     population skewness (scalar)  
##            
#########################################################################
pop_skew = function(mean_by_type,var_by_type,skew_by_type,mixdist){
		if(min(mixdist)<0) {stop("Not a valid mixing distribution")} 
		if(sum(mixdist)!=1) {stop("Not a valid mixing distribution")}
		if(min(var_by_type)<0) {stop("Not a valid variances vector")}
		nZ = length(mean_by_type) 
		if(length(var_by_type)!=nZ) {stop("Length mismatch, mean and var")}
		if(length(skew_by_type)!=nZ) {stop("Length mismatch, mean and mu3")}
		if(length(mixdist)!=nZ) {stop("Length mismatch, means and mixing distribution")}
		
		mu3_by_type = skew_by_type*(var_by_type^(3/2)); 
		mu3 = pop_mu3(mean_by_type,var_by_type,mu3_by_type,mixdist)
		mu2 = pop_mean_var(mean_by_type,var_by_type,mixdist)$pop_var 
		
		result = mu3/(mu2^(3/2)); 
		return(result)
}