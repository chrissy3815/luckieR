## Other matrix utilities


#' Calculate the fundamental matrix
#'
#' The fundamental matrix, \eqn{\mathbf{N}}, of a Markov Chain contains the
#' expected number of visits to each of the transient states before
#' transitioning into an absorbing state. It is given by:
#' \deqn{\mathbf{N}=(\mathbf{I}-\mathbf{P})^{-1},} where \eqn{\mathbf{P}} is the
#' matrix of transition probabilities among the transient states, and the power
#' \eqn{-1} indicates the matrix inverse.
#'
#' @param P An \eqn{n \times n} matrix of transition rates among the \eqn{n}
#'   transient states in the matrix population model. If the absorbing state is
#'   death, then this matrix is the `Umat` of transition rates among the living
#'   stages
#'
#' @returns An \eqn{n \times n} matrix where each \eqn{ij} entry is the expected
#'   number of time steps spent in state \eqn{i} given that an individual starts
#'   in state \eqn{j}
#' @export
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' Nmat<- fundamental_matrix(Umat)
fundamental_matrix = function(P){
  solve( diag(ncol(P)) - P)
}

#' Calculate the stable state distribution
#'
#' The stable state distribution is the proportional distribution of the
#' population across states (i.e., ages, stages, or sizes) when the population
#' is growing at its asymptotically-constant equilibrium growth rate
#' \eqn{\lambda}. The asymptotic population growth rate \eqn{\lambda} is given
#' by the eigenvalue with the largest real part. The stable state distribution
#' is the right eigenvector that corresponds to \eqn{\lambda}, normalized to sum
#' to 1.
#'
#' @param Amat An \eqn{n \times n} matrix containing all of the transition and
#'   reproductive rates that project the population from time \eqn{t} to
#'   \eqn{t+1}. Generally given by `Umat+Fmat`.
#'
#' @returns A vector of length \eqn{n} that contains the proportion of
#'   individuals in each of the \eqn{n} states when the population is growing at
#'   its asymptotic population growth
#' @export
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' Fmat<- matrix(c(0.2, 0, 0, 1, 0.5, 0, 5, 3, 0), ncol=3)
#' Amat<- Umat+Fmat
#' ssd<- stable_dist(Amat)
stable_dist<- function(Amat){
  # Calculate the eigen values and vectors:
  eigz<- eigen(Amat)
  # Find the index of lambda:
  I<- which(Re(eigz$values)==max(Re(eigz$values)))
  # Calculate the stable population distribution:
  wmean<- Re(eigen(Amat)$vectors[,I]) # right eigenvector of the input matrix
  # Rescale to sum to 1 (proportions):
  wmean<- wmean/sum(wmean)
  return(wmean)
}

#' The distribution of offspring types in a cohort produced at the stable
#' population distribution
#'
#' If adults in a population can produce multiple types of offspring, then life
#' can start in multiple states. Therefore, calculations of population-level
#' mean, variance, and skewness in lifespan or lifetime reproductive output must
#' take into account this variation in starting states. We often call this the
#' "mixing distribution" and a reasonable standard mixing distribution is the
#' distribution of offspring types in a cohort produced when the population is
#' at its stable population distribution.
#'
#' @param Amat An \eqn{n \times n} matrix containing all of the transition and
#'   reproductive rates that project the population from time \eqn{t} to
#'   \eqn{t+1}. Generally given by `Umat+Fmat`.
#' @param Fmat An \eqn{n \times n} matrix of per capita rates of reproduction
#'   rates by any of the \eqn{n} states into offspring (new individuals of any
#'   of the \eqn{n} states)
#'
#' @returns A vector of length \eqn{n} that contains the proportion of
#'   individuals in each of the \eqn{n} states in a cohort of offspring
#' @export
#'
#' @examples
#' Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
#' Fmat<- matrix(c(0.2, 0, 0, 1, 0.5, 0, 5, 3, 0), ncol=3)
#' Amat<- Umat+Fmat
#' mixdist<- calcDistOffspringCohort(Amat, Fmat)
calcDistOffspringCohort<- function(Amat, Fmat){
  wmean<- stable_dist(Amat)
  # Multiply the stable distribution by Fmat to get a cohort of offspring:
  offspring<- Fmat%*%wmean
  # Rescale to sum to 1 (proportions):
  offspring<- offspring/sum(offspring)

  return(offspring)
}
