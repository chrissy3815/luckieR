###                     (I) OVERVIEW AND CONVENTIONS

#############################################################################
## R translations of formulas for mean/var of conditional lifespan, and
## mean/var of conditional times to reach state i from state j, from
##     Cochran, M. E. and S. Ellner. 1992.  Simple methods for
##     calculating age-based life history parameters for stage-
##     structured populations. Ecological Monographs 62:345-364.
##
## Original coding by Erin Souder supervised by Simone Blomberg,
## Revised and extended for luckieR by Steve Ellner
##
## NOTE on time conventions:
##     CE92 used the convention that all times (lifetime, transition time, etc.)
##     are the supremum of values that are consistent with the census data -
##     see figure 2 of CE92. That convention gives annuals a lifespan of
## 	   2 years, which was generally considered to be a bad idea, so we
##     use a different convention here.
##
##     1. An individual's lifespan is the number of times that they are
##     alive at a census time (rather than that number plus one).
##
##     2. For transition times, a state-j newborn *could have* lived up to one
## 	   step prior to being censused in state j, so the CE92 formula gives 1 as
##	   the conditional mean time to reach state j starting from state j.
##     The functions here give instead 0: if you are first censused in state j,
## 	   we assume that it took you 'no time at all' to get there.
##############################################################################


######################### NOTATION CONVENTIONS  ##############################
##      luckieR M = StageCoach P (state transitions)
##      luckieR F = StageCoach B (sexual births) + F (clonal births = fission)
##############################################################################

##############################################################################
## The defined functions are:
## 			mean_conditional_times, var_conditional_times
##			mean_conditional_lifespan, var_conditional_lifespan
##			first_breed -- prob. of breeding before death, mean & var of age
##                         at first breeding (conditional on breeding)
##			wrappers for calculating population moments from state-specific
##				moments and a mixing distribution: pop_mean_var, pop_mu3, pop_skew
##			Di_mat, a utility function of no independent interest (computes
##				modified kernel in which entering state i is followed by death).
##############################################################################

#' Modified Transition Matrices \eqn{D_i} for Each Focal State
#'
#' Construct, for each state \eqn{i}, a modified transition matrix \eqn{D_i}
#' in which individuals entering state \eqn{i} remain there for one time step
#' and then die. This is a utility function implementing eq. 8 in
#' Cochran and Ellner (1992), and is not intended to be called directly by
#; package users.
#'
#' @details
#' Given a square, numeric state transition matrix \eqn{M} (with columns
#' representing destination states and rows representing origin states),
#' this function returns a 3D array \eqn{D} where the third dimension indexes
#' the focal state \eqn{i}. For each \eqn{i}, \eqn{D_{\,\cdot\,\cdot\,i}} is a
#' copy of \eqn{M} with the \eqn{i}-th column set to 0, so that individuals
#' transition into state \eqn{i} spend one time step there and then die.
#'
#' @param M A square, numeric state transition matrix for a matrix projection
#' model.
#'
#' @return
#' A 3D numeric array \code{D} of dimension \code{c(nrow(M), ncol(M), ncol(M))},
#' where \code{D[,, i]} is the modified transition matrix \eqn{D_i} associated
#' with focal state \eqn{i}.
#'
#' @references
#' Cochran, M. E., & Ellner, S. (1992). Simple methods for calculating age‐based
#' life history parameters for stage‐structured populations. *Ecological Monographs*,
#' 62(3), 345–364.
#'
#' @examples
#' # A simple 2-state (juvenile, adult) model
#' M <- matrix(c(0.0, 0.0,
#'               0.6, 0.7),
#'               nrow = 2, byrow = TRUE)
#'
#' D <- Di_mat(M)
#' dim(D)            # 2 x 2 x 2
#' D[, , 1]          # D_1: column 1 zeroed
#' D[, , 2]          # D_2: column 2 zeroed
#'
#' @keywords internal
#' @noRd
Di_mat <- function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M) != ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}

  results <- array(0, dim = c(dim(M), ncol(M)))
  for (i in 1:ncol(M)) {
    results[,, i] <- M
    results[, i, i] <- 0
  }
  return(results)
}

#' Compute conditional time to arriving in a state
#'
#' Based on Equation 9 from Cochran and Ellner (1992), this function will
#' calculate the expected time for an individual to transition from a starting
#' state \eqn{j} to an arrival state \eqn{i}, conditional on arriving in state
#' \eqn{i}.
#'
#' @param M An \eqn{n \times n} matrix of transition rates among the \eqn{n}
#'   states in the matrix population model
#'
#' @returns An \eqn{n \times n} matrix where the \eqn{(i,j)} element represents
#'   the conditional mean time to reach i from j, with NA if the transition from
#'   j to i is not possible.
#' @export
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' meanCondTimes<- meanConditionalTimes(Umat)
meanConditionalTimes <- function(M){
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

#' Mean Conditional Lifespan
#'
#' Calculates the mean lifespan for individuals starting in state \eqn{j},
#' conditional on them reaching state \eqn{i} before death. This calculation is
#' based on Equation 6 in Cochran and Ellner (1992).
#'
#' @param M An \eqn{n \times n} matrix of transition rates among the \eqn{n}
#'   states in the matrix population model
#'
#' @returns An \eqn{n \times n} matrix where the \eqn{(i,j)} element represents
#'   the remaining expected lifespan of individuals starting in state \eqn{j},
#'   conditional on them reaching state \eqn{i} before death, with NA if the
#'   transition from j to i is not possible.
#' @export
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' meanCondLifespan<- meanConditionalLifespan(Umat)
meanConditionalLifespan <- function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}

  LE <- meanLifespan(M)
  MT <- meanConditionalTimes(M)

  results <- matrix(0,nrow = nrow(M),ncol = ncol(M))
  for (i in 1:nrow(M)) {
	results[i,]=MT[i,] + LE[i]
  }
  return(results)
}

#' Variance in conditional time to reach state i from state j
#'
#' Calculates the variance in lifespan for individuals starting in state
#' \eqn{j}, conditional on them reaching state \eqn{i} before death. This
#' calculation is based on Equation 10 in Cochran and Ellner (1992).
#'
#' @param M An \eqn{n \times n} matrix of transition rates among the \eqn{n}
#'   states in the matrix population model
#'
#' @returns An \eqn{n \times n} matrix where the \eqn{(i,j)} element represents
#'   the variance in remaining expected lifespan of individuals starting in
#'   state \eqn{j}, conditional on them reaching state \eqn{i} before death,
#'   with NA if the transition from j to i is not possible.
#' @export
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' varCondTimes<- varConditionalTimes(Umat)
varConditionalTimes <- function(M){
  ##   Note: the difference in time conventions (shorter here by 1 than in CE92)
  ##     does not affect the variance, so we use the CE92 convention and formula.
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}

  D  = Di_mat(M)
  tau = meanConditionalTimes(M) + 1 ## use CE92 convention

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

#' Variance of Conditional Lifespan
#'
#' Calculates the variance in remaining lifespan for individuals starting in
#' state \eqn{j}, conditional on them reaching state \eqn{i} before death. This
#' calculation is based on Equation 7 in Cochran and Ellner (1992).
#'
#' @param M An \eqn{n \times n} matrix of transition rates among the \eqn{n}
#'   states in the matrix population model
#'
#' @returns An \eqn{n \times n} matrix where the \eqn{(i,j)} element represents
#'   the variance of remaining lifespan of individuals starting in state
#'   \eqn{j}, conditional on them reaching state \eqn{i} before death, with NA
#'   if the transition from j to i is not possible.
#' @export
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' varCondLifespan<- varConditionalLifespan(Umat)
varConditionalLifespan = function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}

  # variance of conditional time to reach i from j
  var_tau_ij = varConditionalTimes(M)

  # variance of remaining lifespan starting at i
  var_Omega_i = varLifespan(M)

  nx = nrow(M); results = matrix(NA,nx,nx);
  for(i in 1:nx){
  for(j in 1:nx){
	results[i,j] = var_tau_ij[i,j] + var_Omega_i[i]
  }}
  return(results)
}

#' Probability of breeding and age at first breeding
#'
#' Calculates the probability of breeding at least once, and the mean and
#' variance of age at first breeding conditional on breeding, as a function of
#' initial state. In order to use this function, the model specification must
#' include state-dependent breeding probability, not just state-dependent mean
#' fecundity. The probability of breeding at least once from now until death
#' must be positive for all states. If not, remove those states from the model
#' (can't start in any of them, going into them becomes death) and apply this
#' function to the remaining states.
#'
#' @param M0 An \eqn{s \times s} matrix of transition rates among the \eqn{s}
#'   non-breeding states in the population
#' @param pb A vector of length \eqn{s} containing the probability of breeding
#'   (in the next time step) for each state.
#'
#' @returns A list of:
#' * `p_breed` A vector of length \eqn{s} containing the probability of breeding
#'   before death, given the starting state.
#' * `mean_age` A vector of length \eqn{s} containing the mean age of first
#'   reproduction, given the starting state.
#' * `var_age` A vector of length \eqn{s} containing the variance of age of first
#'   reproduction, given the starting state.
#'
#' @export
#'
#' @seealso [probRepro()] if you do not have values for `pb`,
#' [Rage::mature_age()] which computes the mean age at first entering a state
#' with positive mean fecundity.
#'
#' @details In the `luckieR` package, we also provide the `probRepro()`
#'   function, which can be used in scenarios where the distribution of
#'   probability of breeding is not known. Instead, `probRepro` infers the
#'   probability of breeding by assuming that the mean fecundity (provided in
#'   `Fmat`) is the mean of a Poisson, Bernoulli, or fixed reproductive process.
#'
#'   Assumptions and coding follow Ellner et al. (2016) IPM monograph, chapter
#'   3. The state transition matrix M is assumed to have the form M = p_b M_b +
#'   (1-p_b) M_0 to allow for costs of reproduction, but M_b and M_0 can be
#'   equal if there are no costs. Breeding does not necessarily imply producing
#'   any new recruits, for example clutch size conditional on breeding could be
#'   Poisson.
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' pb<- c(0, 0.1, 0.4) # probability of breeding in the next time step, given current state
firstBreed = function(M0,pb) {
  if(!is.matrix(M0)) {stop("This is not a matrix")}
  if(nrow(M0)!=ncol(M0)) {stop("This is not a square matrix")}
  if(!is.numeric(M0)) {stop("This is not a numeric matrix")}
  if(length(pb)!=nrow(M0)) {stop("Length mismatch, breeding probability vector
                                 should have the same length as the dimension of
                                 the transition matrix")}

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

#' Population mean and variance of any attribute
#'
#' Computes the population mean and variance of some attribute X, given the
#' vectors of mean and variance conditional on individual 'type' Z. This
#' substitutes for all of the individual formulas in Table 2 of Cochran and
#' Ellner (1992). For example, if multiple offspring types are possible, this
#' can be used to calculate mean and variance over the possible starting states.
#'
#' @param mean_by_type A vector of length \eqn{z} containing the mean of
#'   attribute X for individuals of types 1 to \eqn{z}
#' @param var_by_type A vector of length \eqn{z} containing the variance of
#'   attribute X for individuals of types 1 to \eqn{z}
#' @param mixdist A vector of length \eqn{z} indicating the frequency
#'   distribution of individuals across the \eqn{z} types
#'
#' @returns A list containing scalar values `pop_mean_X` and `pop_var_X` which
#'   give the mean and variance, respectively, of attribute X across the
#'   population
#' @export
#'
#' @seealso [meanLifespan()], [meanLRO()], [varLifespan()], and [varLRO()],
#'   which should give the same results when the user provides a mixing
#'   distribution.
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' expectedVal<- meanLifespan(Umat)
#' varVal<- varLifespan(Umat)
#' mixdist<- c(0.3, 0.7, 0)
#' popVals<- popMeanVar(expectedVal, varVal, mixdist)
popMeanVar = function(mean_by_type,var_by_type,mixdist){
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


