## These functions can be used to perform the luck calculations as presented in
## Hernandez et al. 2024 - The natural history of luck https://doi.org/10.1111/ele.14390



#' The mean of lifetime reproductive output (LRO)
#'
#' Calculates the expected value of lifetime reproductive output, or lifetime
#' reproductive success, across individuals in a population whose dynamics can
#' be described a matrix population model. The survival and transition rates
#' among stages must be provided in the `Umat` and the mean per capita
#' reproductive output (per time step of the population model) must be provided
#' in the `Fmat`. The default behavior is to return the expected value of LRO
#' for individuals starting in all states of the population. This function also
#' allows the user to provide a "mixing distribution" of starting states, and/or
#' to provide offspring weights such that offspring of different types
#' contribute differently to reproductive success.
#'
#' @param Umat An \eqn{n \times n} matrix of transition rates among the \eqn{n}
#'   states in the matrix population model
#' @param Fmat An \eqn{n \times n} matrix of per capita rates of reproduction
#'   rates by any of the \eqn{n} states into offspring (new individuals of any
#'   of the \eqn{n} states)
#' @param mixdist Optional, a vector of length \eqn{n} that contains the
#'   proportion of individuals starting in each state, to calculate a weighted
#'   average LRO over the entire population
#' @param offspring_weight Optional, a vector of length \eqn{n} that contains
#'   the relative value of offspring of different types. If included, the `Fmat`
#'   is scaled column-wise by these weights, such that the rewards accrued
#'   differ across offspring types. For more details and an example of usage,
#'   see Hernandez et al. 2024. The natural history of luck: a synthesis study
#'   of structured population models. Ecology Letters, 27, e14390. Available
#'   from: https://doi.org/10.1111/ele.14390
#'
#' @returns Value(s) of expected lifetime reproductive output. If
#'   `mixdist=NULL`, then `meanLRO` returns a vector containing the expected
#'   future LRO for individuals starting in each of the \eqn{n} states. If a
#'   mixing distribution is provided, then `meanLRO` returns a single value,
#'   calculated as a weighted average of the offspring starting states.
#' @export
#' @seealso [calcDistOffspringCohort()]
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' Fmat<- matrix(c(0.2, 0, 0, 1, 0.5, 0, 5, 3, 0), ncol=3)
#' expectedVal<- meanLRO(Umat, Fmat)
#'
#' # A standard choice for the mixing distribution is the distribution of
#' # offspring types in a cohort produced by the population at its stable
#' # distribution.
#' Amat<- Umat+Fmat
#' mixdist<- calcDistOffspringCohort(Amat, Fmat)
#' expectedVal<- meanLRO(Umat, Fmat, mixdist)
#'
#' # In Hernandez et al. 2024, we used the probability of surviving to reproduce
#' # to rescale reproductive rewards when individuals can produce multiple
#' # offspring types.
#' off_wts<- probRepro(Umat, Fmat, repro_var='poisson')
#' expectedVal<- meanLRO(Umat, Fmat, offspring_weight=off_wts)
meanLRO<- function(Umat, Fmat, mixdist=NULL, offspring_weight=NULL){
  # check that Umat and Fmat are the same size, and that they are both square:
  if (dim(Umat)[1]==dim(Fmat)[1] & dim(Umat)[2]==dim(Fmat)[2]){
    if (dim(Umat)[1]!=dim(Umat)[2]){
      warning('Umat and Fmat are not square matrices.')
    }
  } else {warning('Umat and Fmat are not the same size.')}
  Nclasses<- dim(Umat)[1]

  if (!is.null(offspring_weight)){
    Fmat<- sweep(Fmat, MARGIN = 1, offspring_weight, FUN="*")
  }

  ## Calculate Ex(R | birth size z)
  N<- fundamental_matrix(Umat)
  expRCond_z<- rep(1,Nclasses)%*%Fmat%*%N

  if(!is.null(mixdist)){
    expR<- expRCond_z%*%mixdist
    return(expR)
  } else{
    return(expRCond_z)
  }
}

#' Probability of reproducing at least once, based on starting state
#'
#' Calculates the probability of reproducing at least once, for individuals
#' starting in all possible starting states. In Hernandez et al. 2024, we used
#' the probability of surviving to reproduce to rescale reproductive rewards
#' when individuals can produce multiple offspring types.
#'
#' @param Umat An \eqn{n \times n} matrix of transition rates among the \eqn{n}
#'   states in the matrix population model
#' @param Fmat An \eqn{n \times n} matrix of per capita rates of reproduction
#'   rates by any of the \eqn{n} states into offspring (new individuals of any
#'   of the \eqn{n} states)
#' @param repro_var The form of variance in reproductive output, with possible
#'   values of `"poisson"`, `"bernoulli"` or `"fixed"`. This is required to
#'   calculate the probability of breeding this year from the values of mean per
#'   capita reproductive output.
#'
#' @returns A vector of length \eqn{n} containing the probability of reproducing
#'   at least once, for individuals starting in each of the \eqn{n} states
#' @export
#' @importFrom stats ppois
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' Fmat<- matrix(c(0.2, 0, 0, 1, 0.5, 0, 5, 3, 0), ncol=3)
#' probVals<- probRepro(Umat, Fmat, repro_var='poisson')
probRepro<- function(Umat, Fmat, repro_var = 'poisson'){
  # take the column sum of reproduction:
  betabar<- colSums(Fmat) # the average offspring production

  # Breeding probability values, p_b:
  if (repro_var %in% c("poisson", "Poisson")){
    p_b<- 1-ppois(0, lambda=betabar) # poisson: prob(x>0)
  } else if (repro_var %in% c("Bernoulli", "bernoulli")){
    p_b<- betabar # bernoulli: prob(occurrence)
  } else if (repro_var %in% c("fixed", "Fixed")){
    p_b<- ifelse(betabar>0, 1, 0) # all individuals in that size reproduce the exact number of offspring shown
  }

  P_0<- sweep(Umat, MARGIN = 2, (1-p_b), '*') # survival matrix with reproduction as an absorbing state
  N_0<- fundamental_matrix(P_0)

  # probability of reproducing at least once:
  B<- p_b%*%N_0
  return(B)
}
