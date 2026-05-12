## These functions can be used to perform the luck calculations as presented in
## Hernandez et al. 2024 - The natural history of luck https://doi.org/10.1111/ele.14390

#' Calculate mean lifetime reproductive output (LRO)
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

#' Calculate variance in lifetime reproductive output
#'
#' Calculates the variance of lifetime reproductive output, or lifetime
#' reproductive success, across individuals in a population whose dynamics can
#' be described a matrix population model. The survival and transition rates
#' among stages must be provided in the `Umat` and the mean per capita
#' reproductive output (per time step of the population model) must be provided
#' in the `Fmat`. The default behavior is to assume Poisson-distributed annual
#' reproduction among individuals in a given stage class. The default is to
#' return the variance of LRO for individuals starting in all states of the
#' population. This function also allows the user to provide a "mixing
#' distribution" of starting states, and/or to provide offspring weights such
#' that offspring of different types contribute differently to reproductive
#' success.
#'
#' @param Umat An \eqn{n \times n} matrix of transition rates among the \eqn{n}
#'   states in the matrix population model
#' @param Fmat An \eqn{n \times n} matrix of per capita rates of reproduction
#'   rates by any of the \eqn{n} states into offspring (new individuals of any
#'   of the \eqn{n} states)
#' @param repro_var The form of variance in per capita reproductive output, with
#'   possible values of `"poisson"`, `"bernoulli"` or `"fixed"`.
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
#' @returns Value(s) of variance in lifetime reproductive output. If
#'   `mixdist=NULL`, then `varLRO` returns a vector containing the variance in
#'   future LRO for individuals starting in each of the \eqn{n} states. If a
#'   mixing distribution is provided, then `varLRO` returns a single value of
#'   variance, calculated over the offspring starting states using the law of
#'   total variance.
#'
#' @export
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' Fmat<- matrix(c(0.2, 0, 0, 1, 0.5, 0, 5, 3, 0), ncol=3)
#' variance<- varLRO(Umat, Fmat, repro_var='poisson')
#'
#' # A standard choice for the mixing distribution is the distribution of
#' # offspring types in a cohort produced by the population at its stable
#' # distribution.
#' Amat<- Umat+Fmat
#' mixdist<- calcDistOffspringCohort(Amat, Fmat)
#' variance_overall<- varLRO(Umat, Fmat, repro_var='poisson', mixdist=mixdist)
#'
#' # In Hernandez et al. 2024, we used the probability of surviving to reproduce
#' # to rescale reproductive rewards when individuals can produce multiple
#' # offspring types.
#' off_wts<- probRepro(Umat, Fmat, repro_var='poisson')
#' variance<- varLRO(Umat, Fmat, repro_var='poisson', offspring_weight=off_wts)
varLRO<- function(Umat, Fmat, repro_var = 'poisson', mixdist=NULL, offspring_weight=NULL){
  # quick check that Umat and Fmat are the same size, and that they are both square:
  if (dim(Umat)[1]==dim(Fmat)[1] & dim(Umat)[2]==dim(Fmat)[2]){
    if (dim(Umat)[1]!=dim(Umat)[2]){
      warning('Umat and Fmat are not square matrices.')
    }
  } else {warning('Umat and Fmat are not the same size.')}
  Nclasses<- dim(Umat)[1]

  # Add offspring weights:
  if (!is.null(offspring_weight)){
    Fmat<- sweep(Fmat, MARGIN = 1, offspring_weight, FUN="*")
  }

  # Build the Markov Chain model:
  mortrow<- 1-colSums(Umat) # probability of mortality
  Pmat<- rbind(Umat, mortrow, deparse.level = 0) # add the mortality row
  Pmat<- cbind(Pmat, 0) # add a column of 0's for the dead individuals
  Pmat[Nclasses+1, Nclasses+1]<- 1 # Death is an absorbing state, individuals cannot leave it

  # Z matrix operator:
  Zmat<- cbind(diag(1, Nclasses, Nclasses, names=FALSE), 0)
  # column matrix of 1's:
  OneVec<- rep(1, Nclasses+1)
  # The fundamental matrix:
  Nmat<- fundamental_matrix(Umat)
  # Take the sum of offspring:
  eff<- colSums(Fmat) #eff is the stage-specific offspring production

  # Calculate the raw moments of the reward matrix:
  R1<- OneVec %*% t(c(eff, 0)) # first moment is the stage-specific reproductive output
  if (repro_var %in% c("poisson", "Poisson")){
    R2<- R1 + (R1*R1)
    R3<- R1 + 3*(R1*R1) + (R1*R1*R1)
  } else if (repro_var %in% c("fixed", "Fixed")){
    R2<- R1*R1
    R3<- R1*R1*R1
  } else if (repro_var %in% c("bernoulli", "Bernoulli")){
    R2<- R1
    R3<- R1
  } else {
    stop("Reproductive random variable type not recognized. The available options are Poisson, Fixed, and Bernoulli")
  }
  # Rtilde1 is the first raw moment of the reward matrix for only transient states
  Rtilde1<- Zmat%*%R1%*%t(Zmat)
  # Rtilde2 is the second raw moment of the reward matrix for only transient states
  Rtilde2<- Zmat%*%R2%*%t(Zmat)

  # Calculate the raw moments of LRO conditional on starting state:
  mu_prime1<- t(Nmat)%*%Zmat %*% t(Pmat*R1)%*%OneVec
  mu_prime2<- t(Nmat)%*% (Zmat %*% t(Pmat*R2)%*%OneVec + 2*t(Umat*Rtilde1) %*% mu_prime1)
  # mu_prime3<- t(Nmat)%*% (Zmat %*% t(Pmat*R3)%*%OneVec + 3*t(Umat*Rtilde2) %*% mu_prime1 + 3*t(Umat*Rtilde1) %*% mu_prime2)

  ## The first raw moment of LRO conditional on starting state is the expected value:
  expRCond_z<- mu_prime1

  ## Var(R | birth size z)
  varRCond_z<- mu_prime2 - (mu_prime1*mu_prime1)

  if(is.null(mixdist)){
    return(varRCond_z)
  } else{
    # variance in LRO due to differences along trajectories:
    varR_within<- t(mixdist)%*%varRCond_z
    # variance in LRO due to differences among starting states:
    varR_between<- t(mixdist)%*%(expRCond_z^2) - (t(mixdist)%*%expRCond_z)^2
    # total variance in LRO, given the mixing distribution:
    varR<- varR_within + varR_between
    return(varR)
  }
}

#' Calculate skew in lifetime reproductive output
#'
#' Calculates the skew of lifetime reproductive output, or lifetime reproductive
#' success, across individuals in a population whose dynamics can be described a
#' matrix population model. The survival and transition rates among stages must
#' be provided in the `Umat` and the mean per capita reproductive output (per
#' time step of the population model) must be provided in the `Fmat`. The
#' default behavior is to assume Poisson-distributed annual reproduction among
#' individuals in a given stage class. The default is to return the skewness of
#' LRO for individuals starting in all states of the population. This function
#' also allows the user to provide a "mixing distribution" of starting states,
#' and/or to provide offspring weights such that offspring of different types
#' contribute differently to reproductive success.
#'
#' @param Umat An \eqn{n \times n} matrix of transition rates among the \eqn{n}
#'   states in the matrix population model
#' @param Fmat An \eqn{n \times n} matrix of per capita rates of reproduction
#'   rates by any of the \eqn{n} states into offspring (new individuals of any
#'   of the \eqn{n} states)
#' @param repro_var The form of variance in per capita reproductive output, with
#'   possible values of `"poisson"`, `"bernoulli"` or `"fixed"`.
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
#' @returns Value(s) of skewness in lifetime reproductive output. If
#'   `mixdist=NULL`, then `skewLRO` returns a vector containing the skewness in
#'   future LRO for individuals starting in each of the \eqn{n} states. If a
#'   mixing distribution is provided, then `skewLRO` returns a single value of
#'   skewness, calculated over the offspring starting states using the law of
#'   total cumulance
#' @export
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' Fmat<- matrix(c(0.2, 0, 0, 1, 0.5, 0, 5, 3, 0), ncol=3)
#' skewness<- skewLRO(Umat, Fmat, repro_var='poisson')
#'
#' # A standard choice for the mixing distribution is the distribution of
#' # offspring types in a cohort produced by the population at its stable
#' # distribution.
#' Amat<- Umat+Fmat
#' mixdist<- calcDistOffspringCohort(Amat, Fmat)
#' skew_overall<- skewLRO(Umat, Fmat, repro_var='poisson', mixdist=mixdist)
#'
#' # In Hernandez et al. 2024, we used the probability of surviving to reproduce
#' # to rescale reproductive rewards when individuals can produce multiple
#' # offspring types.
#' off_wts<- probRepro(Umat, Fmat, repro_var='poisson')
#' skewness<- skewLRO(Umat, Fmat, repro_var='poisson', offspring_weight=off_wts)
skewLRO<- function(Umat, Fmat, repro_var = 'poisson', mixdist=NULL, offspring_weight=NULL){
  # quick check that Umat and Fmat are the same size, and that they are both square:
  if (dim(Umat)[1]==dim(Fmat)[1] & dim(Umat)[2]==dim(Fmat)[2]){
    if (dim(Umat)[1]!=dim(Umat)[2]){
      warning('Umat and Fmat are not square matrices.')
    }
  } else {warning('Umat and Fmat are not the same size.')}
  Nclasses<- dim(Umat)[1]

  # Add offspring weights:
  if (!is.null(offspring_weight)){
    Fmat<- sweep(Fmat, MARGIN = 1, offspring_weight, FUN="*")
  }

  # Build the Markov Chain model:
  mortrow<- 1-colSums(Umat) # probability of mortality
  Pmat<- rbind(Umat, mortrow, deparse.level = 0) # add the mortality row
  Pmat<- cbind(Pmat, 0) # add a column of 0's for the dead individuals
  Pmat[Nclasses+1, Nclasses+1]<- 1 # Death is an absorbing state, individuals cannot leave it

  # Z matrix operator:
  Zmat<- cbind(diag(1, Nclasses, Nclasses, names=FALSE), 0)
  # column matrix of 1's:
  OneVec<- rep(1, Nclasses+1)
  # The fundamental matrix:
  Nmat<- fundamental_matrix(Umat)
  # Take the sum of offspring:
  eff<- colSums(Fmat) #eff is the stage-specific offspring production

  # Calculate the raw moments of the reward matrix:
  R1<- OneVec %*% t(c(eff, 0)) # first moment is the stage-specific reproductive output
  if (repro_var %in% c("poisson", "Poisson")){
    R2<- R1 + (R1*R1)
    R3<- R1 + 3*(R1*R1) + (R1*R1*R1)
  } else if (repro_var %in% c("fixed", "Fixed")){
    R2<- R1*R1
    R3<- R1*R1*R1
  } else if (repro_var %in% c("bernoulli", "Bernoulli")){
    R2<- R1
    R3<- R1
  } else {
    stop("Reproductive random variable type not recognized. The available options are Poisson, Fixed, and Bernoulli")
  }

  # Rtilde1 is the first raw moment of the reward matrix for only transient states
  Rtilde1<- Zmat%*%R1%*%t(Zmat)
  # Rtilde2 is the second raw moment of the reward matrix for only transient states
  Rtilde2<- Zmat%*%R2%*%t(Zmat)

  # Calculate the raw moments of LRO conditional on starting state:
  mu_prime1<- t(Nmat)%*%Zmat %*% t(Pmat*R1)%*%OneVec
  mu_prime2<- t(Nmat)%*% (Zmat %*% t(Pmat*R2)%*%OneVec + 2*t(Umat*Rtilde1) %*% mu_prime1)
  mu_prime3<- t(Nmat)%*% (Zmat %*% t(Pmat*R3)%*%OneVec + 3*t(Umat*Rtilde2) %*% mu_prime1 + 3*t(Umat*Rtilde1) %*% mu_prime2)

  ## The first raw moment of LRO conditional on starting state is the expected value:
  expRCond_z<- mu_prime1

  ## Second central moment: mu2 = Var(R | birth size z)
  varRCond_z<- mu_prime2 - (mu_prime1*mu_prime1)

  # third central moment of LRO, conditional on Z:
  mu3Cond_z<- mu_prime3 - 3*(mu_prime2*mu_prime1) + 2*(mu_prime1*mu_prime1*mu_prime1)
  # Skewness, conditional on Z:
  skewCond_z<- varRCond_z^(-3/2)*mu3Cond_z

  if(is.null(mixdist)){
    return(skewCond_z)
  } else{ # law of total cumulance for third central moment
    # expected value of third central moment:
    mu3_within<- t(mu3Cond_z) %*% mixdist
    # third central moment of the expected value:
    mu3_between<- t((mu_prime1-(t(mu_prime1)%*%mixdist)[1])^3) %*% mixdist
    # covariance of expected value and variance:
    covExpVar<- t(mu_prime1*varRCond_z)%*%mixdist - (t(mu_prime1)%*%mixdist) * (t(varRCond_z)%*%mixdist)
    # total third central moment:
    mu3_total<- mu3_within + mu3_between + 3*covExpVar

    # total variance:
    var_total<- t(varRCond_z)%*%mixdist + t(mu_prime1^2)%*%mixdist - (t(mu_prime1)%*%mixdist)^2

    # total skewness:
    skew_total<- (var_total)^(-3/2)*mu3_total
    return(skew_total)
  }
}

#' Calculate mean lifespan
#'
#' Calculates the expected value of lifespan across individuals in a population
#' whose dynamics can be described a matrix population model. The survival and
#' transition rates among stages must be provided in the `Umat`. The default
#' behavior is to return the expected remaining lifespan for individuals
#' starting in all states of the population. This function also allows the user
#' to provide a "mixing distribution" of starting states. Note that in the
#' Markov Chain with Rewards framework, moments of lifespan are a special case
#' of moments of LRO, where individuals accrue fixed rewards of 1 for each
#' timestep that they survive.
#'
#' @param Umat An \eqn{n \times n} matrix of transition rates among the \eqn{n}
#'   states in the matrix population model
#' @param mixdist Optional, a vector of length \eqn{n} that contains the
#'   proportion of individuals starting in each state, to calculate a weighted
#'   average lifespan over the possible starting states
#'
#' @returns Value(s) of expected remaining lifespan in terms of number of model
#'   timesteps. If `mixdist=NULL`, then `meanLifespan` returns a vector
#'   containing the expected remaining lifespan for individuals starting in each
#'   of the \eqn{n} states. If a mixing distribution is provided, then
#'   `meanLifespan` returns a single value, calculated as a weighted average
#'   over the offspring starting states.
#' @export
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' expectedVal<- meanLifespan(Umat)
#'
#' # A standard choice for the mixing distribution is the distribution of
#' # offspring types in a cohort produced by the population at its stable
#' # distribution.
#' Fmat<- matrix(c(0.2, 0, 0, 1, 0.5, 0, 5, 3, 0), ncol=3)
#' Amat<- Umat+Fmat
#' mixdist<- calcDistOffspringCohort(Amat, Fmat)
#' expectedVal<- meanLifespan(Umat, mixdist=mixdist)
meanLifespan<- function(Umat, mixdist=NULL){
  # quick check that Umat is square:
  if (dim(Umat)[1]!=dim(Umat)[2]){
    warning('Umat and Fmat are not square matrices.')
  }
  Nclasses<- dim(Umat)[1]

  ## Calculate Ex(R | current state)
  N<- fundamental_matrix(Umat)
  expLCond_z<- rep(1,Nclasses)%*%N

  if(!is.null(mixdist)){
    expL<- expLCond_z%*%mixdist
    return(expL)
  } else{
    return(expLCond_z)
  }
}

#' Calculate variance of lifespan
#'
#' Calculates the variance of lifespan across individuals in a population whose
#' dynamics can be described a matrix population model. The survival and
#' transition rates among stages must be provided in the `Umat`. The default
#' behavior is to return the variance of remaining lifespan for individuals
#' starting in all states of the population. This function also allows the user
#' to provide a "mixing distribution" of starting states. Note that in the
#' Markov Chain with Rewards framework, moments of lifespan are a special case
#' of moments of LRO, where individuals accrue fixed rewards of 1 for each
#' timestep that they survive.
#'
#' @param Umat An \eqn{n \times n} matrix of transition rates among the \eqn{n}
#'   states in the matrix population model
#' @param mixdist Optional, a vector of length \eqn{n} that contains the
#'   proportion of individuals starting in each state, to calculate a weighted
#'   average lifespan over the possible starting states
#'
#' @returns Value(s) of variance in remaining lifespan. If `mixdist=NULL`, then
#'   `varLifespan` returns a vector containing the variance in remaining
#'   lifespan for individuals starting in each of the \eqn{n} states. If a
#'   mixing distribution is provided, then `varLifespan` returns a single value
#'   of variance in remaining lifespan, calculated over the starting states
#'   according to the law of total variance.
#' @export
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' variance<- varLifespan(Umat)
#'
#' # A standard choice for the mixing distribution is the distribution of
#' # offspring types in a cohort produced by the population at its stable
#' # distribution.
#' Fmat<- matrix(c(0.2, 0, 0, 1, 0.5, 0, 5, 3, 0), ncol=3)
#' Amat<- Umat+Fmat
#' mixdist<- calcDistOffspringCohort(Amat, Fmat)
#' variance<- varLifespan(Umat, mixdist=mixdist)
varLifespan<- function(Umat, mixdist=NULL){
  # quick check that Umat is square:
  if (dim(Umat)[1]!=dim(Umat)[2]){
    warning('Umat and Fmat are not square matrices.')
  }
  Nclasses<- dim(Umat)[1]

  # calculate the fundamental matrix
  N<- fundamental_matrix(Umat)

  ## Calculate Ex(R | current state)
  expLCond_z<- meanLifespan(Umat, mixdist = NULL)

  ## Var(L | current state) using eqn. 5.12 from Hal's book:
  eT<- matrix(data=1, ncol=Nclasses, nrow=1) # column vector of 1's
  varLCond_z<- eT %*% (2*N%*%N - N) - (expLCond_z)^2

  if(is.null(mixdist)){
    return(varLCond_z)
  } else{
    # variance in LRO due to differences along trajectories:
    varL_within<- varLCond_z %*% mixdist
    # variance in LRO due to differences among starting states:
    varL_between<- t(mixdist)%*%t(expLCond_z^2) - (t(mixdist)%*%t(expLCond_z))^2
    # total variance in lifespan, given the mixing distribution:
    varL<- varL_within + varL_between
    return(varL)
  }
}

#' Calculate skewness of lifespan
#'
#' Calculates the skewness of lifespan across individuals in a population whose
#' dynamics can be described a matrix population model. The survival and
#' transition rates among stages must be provided in the `Umat`. The default
#' behavior is to return the skewness of remaining lifespan for individuals
#' starting in all states of the population. This function also allows the user
#' to provide a "mixing distribution" of starting states. Note that in the
#' Markov Chain with Rewards framework, moments of lifespan are a special case
#' of moments of LRO, where individuals accrue fixed rewards of 1 for each
#' timestep that they survive.
#'
#' @param Umat An \eqn{n \times n} matrix of transition rates among the \eqn{n}
#'   states in the matrix population model
#' @param mixdist Optional, a vector of length \eqn{n} that contains the
#'   proportion of individuals starting in each state, to calculate a weighted
#'   average lifespan over the possible starting states
#'
#' @returns Value(s) of skewness in remaining lifespan. If `mixdist=NULL`, then
#'   `skewLifespan` returns a vector containing the skewness in remaining
#'   lifespan for individuals starting in each of the \eqn{n} states. If a
#'   mixing distribution is provided, then `skewLifespan` returns a single value
#'   of skewness in remaining lifespan, calculated over the starting states
#'   according to the law of total cumulance.
#' @export
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' skewness<- skewLifespan(Umat)
#'
#' # A standard choice for the mixing distribution is the distribution of
#' # offspring types in a cohort produced by the population at its stable
#' # distribution.
#' Fmat<- matrix(c(0.2, 0, 0, 1, 0.5, 0, 5, 3, 0), ncol=3)
#' Amat<- Umat+Fmat
#' mixdist<- calcDistOffspringCohort(Amat, Fmat)
#' skewness<- skewLifespan(Umat, mixdist=mixdist)
skewLifespan<- function(Umat, mixdist=NULL){
  # quick check that Umat is square:
  if (dim(Umat)[1]!=dim(Umat)[2]){
    warning('Umat and Fmat are not square matrices.')
  }
  Nclasses<- dim(Umat)[1]

  # calculate the fundamental matrix
  N<- fundamental_matrix(Umat)
  # row vector of 1's
  eT<- matrix(data=1, ncol=Nclasses, nrow=1)
  # identity matrix
  eye<- diag(x=1, nrow=Nclasses, ncol=Nclasses, names=FALSE)

  # calculate the moments of lifespan
  eta1<- meanLifespan(Umat)
  eta2<- eta1%*%(2*N-eye)
  eta3<- eta1%*%(6*N%*%N-6*N+eye)

  # calculate the central moments:
  eta_hat2<- eta2-eta1*eta1
  eta_hat3<- eta3 - 3*eta1*eta2 + 2*eta1*eta1*eta1

  # calculate skew:
  skewCond_z<- eta_hat2^(-3/2)*eta_hat3

  if(is.null(mixdist)){
    return(skewCond_z)
  } else{ # law of total cumulance for third central moment
    # expected value of third central moment:
    eta_hat3_within<- eta_hat3 %*% mixdist
    # third central moment of the expected value:
    eta_hat3_between<- ((eta1-(eta1%*%mixdist)[1])^3) %*% mixdist
    # covariance of expected value and variance:
    covExpVar<- (eta1*eta_hat2)%*%mixdist - (eta1%*%mixdist) * (eta_hat2%*%mixdist)
    # total third central moment:
    eta_hat3_total<- eta_hat3_within + eta_hat3_between + 3*covExpVar

    # total variance:
    var_total<- eta_hat2%*%mixdist + (eta1)^2%*%mixdist - (eta1%*%mixdist)^2

    # total skewness:
    skew_total<- (var_total)^(-3/2)*eta_hat3_total
    return(skew_total)
  }
}
