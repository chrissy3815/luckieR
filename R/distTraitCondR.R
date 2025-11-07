#' Distribution of lifespan conditional on LRO
#'
#' Calculates the distribution of lifespan conditional on LRO in the
#' presence of environmental variation.
#'
#' @param Plist A list of survival/growth transition matrices.
#'   Plist\[\[q\]\]\[i,j\] is the probability of transitioning from
#'   state j to state i in environment q.
#' @param Flist A list of fecundity matrices.  Flist\[\[q\]\]\[i,j\]
#'   is the expected number of state i offspring from a state j parent
#'   in environment q
#' @param Q The environment transition matrix.  Q\[i,j\] is the
#'   probability of transitioning from environment j to environment i.
#' @param c0 A vector specifying the offspring state distribution:
#'   c0\[j\] is the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param maxAge The maximum attainable age
#' @param percentileCutoff A value between 0 and 1.  Calculations are
#'   performed for values of LRO out to this percentile.  Optional.
#'   The default value is 0.99. 
#' @param Fdist The clutch size distribution.  The recognized options
#'   are "Poisson" and "Bernoulli".  Optional.  The default value is
#'   "Poisson".
#' @details The details of this calculation can be found in Robin
#'   E. Snyder and Stephen P. Ellner.
#'   2024. "To prosper, live long: Understanding the sources of
#'   reproductive skew and extreme reproductive success in structured
#'   populations."  The American Naturalist 204(2),
#'   https://doi.org/10.1086/730557 and its online supplement. 
#' @return Returns a list containing the following:
#' * probLifespanCondR: A matrix whose \[i,j\]th entry is the
#'   probability that an individual with LRO = i-1 will have a
#'   lifespan of j years (i.e. age j-1).
#' * distKidsAtDeath: A vector whose jth entry is the probability of
#'   having an LRO of j-1.
#' * sdLifespanCondR: A vector whose jth entry is the standard
#'   deviation of lifespan conditional on having an LRO of j-1.
#' * CVLifespanCondR: A vector whose jth entry is the coefficient of
#'   variation of lifespan conditional on having an LRO of j-1.
#' * maxKidsIndex: calculations were performed out to values of LRO
#'   equal to maxKidsIndex-1.
#' * normalSurvProb: A vector whose jth entry is the unconditional
#'   probability of surviving j years (i.e. until age j-1).
#' @examples
#' P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' P2 = matrix (c(0, 0.2, 0, 0, 0, 0.3, 0, 0, 0.4), 3, 3)
#' F1 = matrix (0, 3, 3); F1[1,] = 0:2
#' F2 = matrix (0, 3, 3); F1[1,] = 0.8*(0:2)
#' Plist = list (P1, P2)
#' Flist = list (F1, F2)
#' Q = matrix (1/2, 2, 2)
#' c0 = c(1,0,0)
#' maxClutchSize = 20
#' maxAge=20
#' out = distLifespanCondR2 (Plist, Flist, Q, c0, maxClutchSize, maxAge)
distLifespanCondR2 = function (Plist, Flist, Q,
                               c0, maxClutchSize, maxAge,
                               percentileCutoff = 0.99,
                               Fdist="Poisson") {
  require (Matrix)
  debugging = TRUE
  
  mT = maxClutchSize + 1
  mA = maxAge + 1
  
  mz = dim(Plist[[1]])[1]
  numEnv = dim(Q)[1]
  bigmz = numEnv*mz

  ## u0 is the stationary environmental distribution, given by the
  ## dominant eigenvector of Q
  u0 = eigen(Q)$vectors[,1]
  u0 = u0 / sum(u0)

  ## m0 is the stationary state cross-classified by size and
  ## environment
  m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)

  ## Define megamatrices M and F
  F = M = matrix (0, bigmz, bigmz)
  for (i in 1:numEnv) {
    for (j in 1:numEnv) {
      M[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Plist[[j]]*Q[i,j]
      F[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Flist[[j]]*Q[i,j]
    }
  }

  ## B[i,j] is the probability that a class-j individual has i-1 kids.
  ## We assume Poisson-distributed number of offspring.
  ## The columns of B should sum to 1 and they do.
  B = mk_B (maxClutchSize, F, Fdist)

  ## Construct A, the transition matrix for a (number of kids) x stage x
  ## env. model.   
  out = make_AxT (B, M, mT)
  A = out$A
  mzA = bigmz*mT
  K = out$K

  #####################################################################
  ## making the bullet matrices
  #####################################################################

  ## Make "bullet" matrices for A.
  Fbullet = Qbullet = Gbullet = matrix (0, mzA, mzA)

  ## Sbullet updates survival
  Sbullet = diag (colSums(A))

  ## Fbullet updates number of kids.
  for (j in 1:mT) {
    for (i in j:(mT-1)) {
      for (z in 1:bigmz) {
        if (Fdist == "Poisson") {
          Fbullet[(i-1)*bigmz + z, (j-1)*bigmz + z] =
            dpois (i-j, lambda=sum(F[,z]))
        } else if (Fdist == "Bernoulli") {
          Fbullet[(i-1)*bigmz + z, (j-1)*bigmz + z] =
            dbinom (i-j, prob=sum(F[,z]), size=1)
        } else {
          cat ("Did not recognize Fdist choice.\n")
        }
      }
    }
  }
  ## Make mT absorbing.  If you have #kids index < mT, the probability
  ## of going to mT while in stage z = 1 - all other kid transition
  ## probabilities with stage z.

  for (j in 1:(mT-1)) {
    for (z in 1:bigmz) {
      Fbullet[(mT-1)*bigmz + z, (j-1)*bigmz + z] =
        ##      1 - sum(Fbullet[1:((mT-1)*bigmz), (j-1)*bigmz + z])
        1 - sum(Fbullet[(0:(mT-2))*bigmz + z, (j-1)*bigmz + z])
    }
  }
  ## If you have #kids index = mT, you stay there.
  for (z in 1:bigmz) 
    Fbullet[(mT-1)*bigmz + z, (mT-1)*bigmz + z] = 1

  ## Qbullet updates env.
  smallQbullet = matrix (0, bigmz, bigmz)
  for (i in 1:numEnv) {
    for (j in 1:numEnv) {
      smallQbullet[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = diag(mz)*Q[i,j]
    }
  }
  ## Qbullet is the same for each number of kids
  for (k in 1:mT) 
    Qbullet[(k-1)*bigmz + 1:bigmz, (k-1)*bigmz + 1:bigmz] = smallQbullet

  ## Gbullet updates stage
  smallGbullet = matrix (0, bigmz, bigmz)
  for (q in 1:numEnv) {
    S = diag (colSums (Plist[[q]]))
    Sinv = diag (1 / diag(S))
    G = Plist[[q]] %*% Sinv
    smallGbullet[(q-1)*mz + 1:mz, (q-1)*mz + 1:mz] = G
  }
  for (k in 1:mT) 
    Gbullet[(k-1)*bigmz + 1:bigmz, (k-1)*bigmz + 1:bigmz] = smallGbullet

  ## sanity check
  if (debugging)
    cat ("Checking the bullet matrices for all #kids:",
         range((Qbullet %*% Gbullet %*% Sbullet %*% Fbullet - A)),
         "should = 0.\n")


  ####################################################################
  ## Now make the additionally extended space matrix that allows us to
  ## stop partway through a time step.  Our new prefix is es to
  ## indicate this additional extension.  ("es" = "extended space")
  ####################################################################
  
  esmzA = 4*mzA

  ##esbigF = esA = matrix (0, esmzA, esmzA)
  esA = matrix (0, esmzA, esmzA)
  esA[mzA + 1:mzA, 1:mzA] = Fbullet
  esA[2*mzA + 1:mzA, mzA + 1:mzA] = Sbullet
  esA[3*mzA + 1:mzA, 2*mzA + 1:mzA] = Gbullet
  esA[1:mzA, 3*mzA + 1:mzA] = Qbullet

  ##################################################################
  ## What is the distribution of #kids at death when reproduction
  ## comes before survival and growth?
  ##################################################################

  ## Create a survival probability vector 
  surv = apply (esA, 2, sum);  die = 1-surv; 

  ## And the fundamental matrix of esA
  fundesA <- solve(diag(ncol(esA))-esA)

  ## Make omega
  omega = matrix(0,nrow(fundesA),ncol(fundesA)); 

  for(j in 1:ncol(fundesA))
    omega[,j] = die*fundesA[,j]

  distAtBirth = matrix(0, bigmz, 4*mT)
  distAtBirth[,1] = m0  ## Everyone starts off with zero kids at the
  ## beginning of the time step

  ## Get distribution of states at death
  distAtDeath = omega %*% matrix(distAtBirth, ncol=1)  # state vector
  ## distAtDeath = matrix(distAtDeath, nrow=bigmz) # state matrix
  ## distKidsAtDeath = apply(distAtDeath, 2, sum)[1:mT]
  distAtDeath = array (distAtDeath, dim=c(bigmz, mT, 4)) # state
                                       # matrix
  distKidsAtDeath = apply (distAtDeath, 2, sum)

  ## What threshold number of kids represents the 99th %ile of the LRO
  ## distribution? 
  cumDistKidsAtDeath = cumsum (distKidsAtDeath)
  R90 = which (cumDistKidsAtDeath > 0.9)[1]
  R99 = which (cumDistKidsAtDeath > 0.99)[1]
  if (debugging) {
    cat ("The 90th %ile of the LRO distribution occurs at", R90,
         "kids.\n")
    cat ("The 99th %ile of the LRO distribution occurs at", R99,
         "kids.\n")
  }
  RCutoff = which (cumDistKidsAtDeath > percentileCutoff)[1]
  if (RCutoff > maxClutchSize)
    stop ("RCutoff = ", RCutoff, "but maxClutchSize = ", maxClutchSize, "\n")

  probLifespanCondR = matrix (0, RCutoff+2, mA)
  probThresholdOrMore = matrix (0, RCutoff+2, esmzA)
  probSurvAgeA = numeric (mA+1)
  probSuccessCondZ = matrix (0, RCutoff+1, esmzA)
  probSuccess = numeric (RCutoff+1)
  jointProbLifespanR = probLifespanCondR = matrix (0, RCutoff+1, mA)
  normalSurvProb = numeric (RCutoff+2)

  AaM0 = c(m0, rep(0, (mT-1)*bigmz))
  mzAZero = rep(0, mzA)
  extendedInit = c(AaM0, rep(mzAZero, 3))

  ## The case of R = 0 (kidsIndexThreshold = 1) #########
  probThresholdOrMore[1,] = 1
  esAa = diag (esmzA)
  esA4 = esA %*% esA %*% esA %*% esA
  for (a in 1:(mA+1)) {
    cat ("a = ", a, "\n")
    probSurvAgeA[a] = colSums (esAa) %*% extendedInit
    if (a >= 2)
      probLifespanCondR[1,a-1] =
        probSurvAgeA[a-1] - probSurvAgeA[a]
    ## Keep track of survival probabilities for unconditional kernel.
    ## we need to hit the last matrix 4 times to advance one time step
    esAa = esA4 %*% esAa
  }
  normalSurvProb = probSurvAgeA

  ## The case of R > 0
  ##for (kidsIndexThreshold in 2:(RCutoff+2)) {

  ## k is the actual number of kids for the threshold.  We start with
  ## the threshold for success at LRO = 1 kid.
  for (k in 1:(RCutoff+1)) {
    cat ("#kids threshold = ", k, "\n")

    ## Get esA conditional on success
    cutoff = k*bigmz
    transientStates = c(1:cutoff, mzA + 1:cutoff, 2*mzA + 1:cutoff,
                        3*mzA + 1:cutoff)
    out = makeCondKernel (esA, transientStates)
   ## condEsA = out$MCond
    probThresholdOrMore[k+1,] = out$q2Extended

    probSuccessCondZ[k,] =
      probThresholdOrMore[k,] - probThresholdOrMore[k+1,]
    probSuccess[k] = probSuccessCondZ[k,] %*% extendedInit

    if (debugging)
      cat ("Pr(success) = ", probSuccess[k], "\n")
    
    condEsA = matrix (0, esmzA, esmzA)
    ## For any starting state with #kids > threshold, probSuccessCondZ
    ## is zero and I think condEsA should be zero too.  Only do the
    ## following for z s.t. probSuccessCondZ > 0.
    for (j in transientStates)
      condEsA[,j] = esA[,j] * probSuccessCondZ[k,] /
        probSuccessCondZ[k,j]

    ## Useful later
    condEsA4 = condEsA %*% condEsA %*% condEsA %*% condEsA

    ## Initial distribution conditional on having LRO = exactly
    ## kidsIndexThreshold-1 
    initCondSuccess =
      probSuccessCondZ[k,] * extendedInit
    initCondSuccess = initCondSuccess / sum(initCondSuccess)

    ## What is the lifespan dist. given that you had at least
    ## numKidsThreshold kids?
    condEsAa = diag(esmzA)
    for (a in 1:(mA+1)) {
      cat ("a = ", a, "\n")
      probSurvAgeA[a] = colSums (condEsAa) %*% initCondSuccess
      ##if (a >= 2)
##        probLifespanCondR[k+1,a-1] =
##          probSurvAgeA[a-1] - probSurvAgeA[a]
      ## Keep track of survival probabilities for unconditional kernel.
      ## we need to hit the last matrix 4 times to advance one time step
      condEsAa = condEsA4 %*% condEsAa
    }

    probLifespanCondR[k,] = -diff (probSurvAgeA)
  }  ## end loop over k

  ## Sanity check: passes if maxAge is large enough
  if (debugging)
    cat ("sum_x Pr(lifespan = x | R) is ",
         apply (probLifespanCondR, 1, sum), "should = 1s.\n")

  ## Why are we calculating probLifespanCondR twice? #############

  if (FALSE) {
    ## Get P(L, R)
    for (kidsIndexThreshold in 1:(RCutoff+1))
      jointProbLifespanR[kidsIndexThreshold,] =
        probLifespanCondR[kidsIndexThreshold,] *
        probSuccess[kidsIndexThreshold]

    ## Sanity check
    cat ("sum_{L, R} Pr(L, R) = ", sum(jointProbLifespanR), " should = 1.\n")

    for (kidsIndexThreshold in 1:(RCutoff+1))
      probLifespanCondR[kidsIndexThreshold,] =
        jointProbLifespanR[kidsIndexThreshold,] /
        distKidsAtDeath[kidsIndexThreshold]
    
    cat ("sum_x Pr(lifespan = x | R) is ",
         apply (probLifespanCondR, 1, sum), "\n")
  }
  
  meanLifespanCondR = sdLifespanCondR = CVLifespanCondR =
    numeric (RCutoff+1)
  
  for (k in 1:(RCutoff+1)) {
    meanLifespanCondR[k] = sum((1:mA)*probLifespanCondR[k,])
    varLifespanCondR = sum((1:mA)^2*probLifespanCondR[k,]) -
      meanLifespanCondR[k]^2
    sdLifespanCondR[k] = sqrt(varLifespanCondR)
    CVLifespanCondR[k] = sdLifespanCondR[k] / meanLifespanCondR[k]
  }

  ## Now let's get Pr(R | Lifespan) (except that I've decided I don't
  ## want that after all, but let's keep it here, just in case I want
  ## it later.)

  ## Passed sanity check: I have compared the output below to
  ## calculating Pr(lifespan) from an age-structured model and they
  ## are identical.
  distLifespan = apply (jointProbLifespanR, 2, sum)

  probRCondLifespan = matrix (0, mA, RCutoff+1)
  for (a in 1:mA)
    probRCondLifespan[a,] = jointProbLifespanR[,a] / distLifespan[a]
  
  return (list(probLifespanCondR=probLifespanCondR,
               distKidsAtDeath=distKidsAtDeath,
               sdLifespanCondR=sdLifespanCondR,
               CVLifespanCondR=CVLifespanCondR,
               maxKidsIndex=RCutoff+1,
               normalSurvProb=normalSurvProb)) }

########################################################################
## calculates the distribution of LRO (#kids at death) in the presence
## of environmental variation, assuming that reproduction happens
## before survival and growth.  This is unlike the version in the IPM Book,
## which effectively assumes that survival (and growth?) happens
## first.  B[i,j] is the probability of a parent of size j having i-1
## offspring.  If B is not provided, then the code constructs B
## according to the clutch size distribution specified by Fdist,
## though currently, only the "Poisson" option is implemented.  (TO
## DO: implement Bernoulli distribution.)
########################################################################

#' Distribution of lifetime reproductive output
#'
#' Calculates the distribution of lifetime reproductive output (LRO)
#' the presence of environmental variation.  Assumes a pre-breeding
#' census (reproductiion happens before survival and growth).
#' @param Plist A list of survival/growth transition matrices.
#'   Plist\[\[q\]\]\[i,j\] is the probability of transitioning from
#'   state j to state i in environment q.
#' @param Flist A list of fecundity matrices.  Flist\[\[q\]\]\[i,j\]
#'   is the expected number of state i offspring from a state j parent
#'   in environment q
#' @param Q The environment transition matrix.  Q\[i,j\] is the
#'   probability of transitioning from environment j to environment i.
#' @param c0 A vector specifying the offspring state distribution:
#'   c0\[j\] is the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param maxLRO The maximum LRO to consider
#' @param Fdist The clutch size distribution.  The recognized options
#'   are "Poisson" and "Bernoulli".  Optional.  The default value is
#'   "Poisson".
#' @details Called by probTraitCondLRO.
#'
#' It is also possible to calculate the distribution of LRO
#' by cross-classifying states by stage and number of offspring
#' produced so far and calculating the state distribution at death.
#' If s(z') is the survival probability for cross-classified state z',
#' M(z', z) is the probability of transitioning from cross-classified
#' state z to cross-classified state z', and N is the fundamental
#' matrix for M: N = (I - M)^\{-1\}, then the state distribution at
#' death, conditional on starting in state z, is (1 - s(z')) * N(z',
#' z), where * denotes a Hadamard product (element-by-element
#' multiplication, not matrix multiplication).  (See, e.g., eq. 3.2.8
#' in Data-driven Modeling of Structured Populations: A Practical
#' Guide to the Integral Projection Model, by Stephen P. Ellner, Dylan
#' Z. Childs and Mark Rees, ed. 1, 2015.)  However, this implicitly assumes
#' that reproduction  happens after survival and growth, i.e. a
#' post-breeding census. 
#' @examples
#' P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' P2 = matrix (c(0, 0.2, 0, 0, 0, 0.3, 0, 0, 0.4), 3, 3)
#' F1 = matrix (0, 3, 3); F1[1,] = 0:2
#' F2 = matrix (0, 3, 3); F1[1,] = 0.8*(0:2)
#' Plist = list (P1, P2)
#' Flist = list (F1, F2)
#' Q = matrix (1/2, 2, 2)
#' c0 = c(1,0,0)
#' maxClutchSize = 10
#' maxLRO=20
#' out = calcDistLRO (Plist, Flist, Q, c0, maxClutchSize, maxLRO)
calcDistLRO = function (Plist, Flist, Q,
                        c0, maxClutchSize, maxLRO,
                        Fdist="Poisson") {
##  require (Matrix)
  debugging = TRUE

  mT = maxLRO + 1
  
  mz = dim(Plist[[1]])[1]
  numEnv = dim(Q)[1]
  bigmz = numEnv*mz

  ## u0 is the stationary environmental distribution, given by the
  ## dominant eigenvector of Q
  u0 = eigen(Q)$vectors[,1]
  u0 = u0 / sum(u0)

  ## m0 is the stationary state cross-classified by size and
  ## environment
  m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)

  ## Define megamatrices M and F
  F = M = matrix (0, bigmz, bigmz)
  for (i in 1:numEnv) {
    for (j in 1:numEnv) {
      M[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Plist[[j]]*Q[i,j]
      F[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Flist[[j]]*Q[i,j]
    }
  }

  ## B[i,j] is the probability that a class-j individual has i-1 kids.
  ## The columns of B should sum to 1.
  B = mk_B (maxClutchSize, F, Fdist)

  ## Sanity check: is maxLRO large enough?
##  cat ("calcDistLRO: x = ", x, ": min(colSums(B)) = ", min(colSums(B)),
##       "\n")

  ## Construct A, the transition matrix for a (number of kids) x stage x
  ## env. model.   
  out = make_AxT (B, M, mT)
  cat ("calcDistLRO: making A...\n")
  A = out$A
  mzA = bigmz*mT
  K = out$K

  #####################################################################
  ## making the bullet matrices
  #####################################################################

  ## Make "bullet" matrices for A.
  Fbullet = matrix (0, mzA, mzA)

  ## BUG: We ought to be able to make Fbullet with just the info. in
  ## B, but you need to respect the fact that i-j can't be bigger than
  ## Fbullet[(i-1)*bigmz + z, (j-1)*bigmz + z] =
  ## B[i-j + 1, z]


  ## Fbullet updates number of kids.
  for (j in 1:mT) {
    for (i in j:(mT-1)) {
      for (z in 1:bigmz) {
        if (Fdist == "Poisson") {
          Fbullet[(i-1)*bigmz + z, (j-1)*bigmz + z] =
            dpois (i-j, lambda=sum(F[,z]))
        } else if (Fdist == "Bernoulli") {
          Fbullet[(i-1)*bigmz + z, (j-1)*bigmz + z] =
            dbinom (i-j, prob=sum(F[,z]), size=1)
        } else {
          cat ("Did not recognize Fdist choice.\n")
        }
      }
    }
  }
  ## Make mT absorbing.  If you have #kids index < mT, the probability
  ## of going to mT while in stage z = 1 - all other kid transition
  ## probabilities with stage z.

  for (j in 1:(mT-1)) {
    for (z in 1:bigmz) {
      Fbullet[(mT-1)*bigmz + z, (j-1)*bigmz + z] =
        ##      1 - sum(Fbullet[1:((mT-1)*bigmz), (j-1)*bigmz + z])
        1 - sum(Fbullet[(0:(mT-2))*bigmz + z, (j-1)*bigmz + z])
    }
  }
  ## If you have #kids index = mT, you stay there.
  for (z in 1:bigmz) 
    Fbullet[(mT-1)*bigmz + z, (mT-1)*bigmz + z] = 1

  ## Mbullet updates everything else
  ## Why doesn't the calculation work when I use bdiag?  No idea.
  if (FALSE) {
    foo = list ()
    for (k in 1:mT)
      foo[[k]] = M
    Mbullet = bdiag (foo)
  } else {
    Mbullet = matrix (0, mzA, mzA)
    for (k in 1:mT)
      Mbullet[(k-1)*bigmz + 1:bigmz, (k-1)*bigmz + 1:bigmz] = M
  }

  
  ## sanity check
  if (debugging)
    cat ("Checking the bullet matrices for all #kids:",
         range((Mbullet %*% Fbullet - A)),
         "should = 0.\n")


  ####################################################################
  ## Now make the additionally extended space matrix that allows us to
  ## stop partway through a time step.  Our new prefix is es to
  ## indicate this additional extension.  ("es" = "extended space")
  ####################################################################

  ## We don't need a space that separately updates survival, growth,
  ## and environment.  We just need to do reproduction first, then
  ## the rest.
  esmzA = 2*mzA
  esA = matrix (0, esmzA, esmzA)
  esA[mzA + 1:mzA, 1:mzA] = Fbullet
  esA[1:mzA, mzA + 1:mzA] = Mbullet

  
  ##################################################################
  ## What is the distribution of #kids at death when reproduction
  ## comes before survival and growth?
  ##################################################################

  ## Create a survival probability vector 
  surv = apply (esA, 2, sum);  die = 1-surv; 

  ## And the fundamental matrix of esA
  cat ("calcDistLRO: calculating fundamental matrix...\n")
  fundesA <- solve(diag(ncol(esA))-esA)

  ## Make omega
  omega = matrix(0,nrow(fundesA),ncol(fundesA)); 

  for(j in 1:ncol(fundesA))
    omega[,j] = die*fundesA[,j]

  ##  distAtBirth = matrix(0, bigmz, 4*mT)
  distAtBirth = matrix(0, bigmz, 2*mT)
  distAtBirth[,1] = m0  ## Everyone starts off with zero kids at the
  ## beginning of the time step

  ## Get distribution of states at death
  distAtDeath = omega %*% matrix(distAtBirth, ncol=1)  # state vector
  distAtDeath = array (distAtDeath, dim=c(bigmz, mT, 2)) # state matrix
  distKidsAtDeath = apply (distAtDeath, 2, sum)

  return (distKidsAtDeath) 
}

#' Distribution of a trait conditional on LRO
#'
#' Calculates Prob(X | R), where X is a trait value and R is lifetime
#' reproductive output (LRO).
#' @param PlistAllTraits A list of lists of survival/growth transition
#'   matrices. Plist\[\[x\]\]\[\[q\]\]\[i,j\] is the probability of
#'   transitioning from state j to state i in environment q when an
#'   individual has trait x. 
#' @param FlistAllTraits A list of lists of fecundity  matrices.
#'   Flist\[\[x\]\]\[\[q\]\]\[i,j\] is the expected number of state i
#'   offspring from a state j, trait x parent in environment q.
#' @param Q The environment transition matrix.  Q\[i,j\] is the
#'   probability of transitioning from environment j to environment
#'   i.
#' @param c0 A vector specifying the offspring state distribution:
#'   c0\[j\] is the probability that an individual is born in state j
#' @param maxClutchSize The maximum clutch size to consider
#' @param maxLRO The maximum LRO to consider
#' @param traitDist A vector whose jth entry is the unconditional
#' probability that the trait takes on the jth value
#' @param Fdist The clutch size distribution.  The recognized options
#'   are "Poisson" and "Bernoulli".  Optional.  The default value is
#'   "Poisson".
#' @details The details of this calculation can be found in Robin
#' E. Snyder and Stephen P. Ellner.  2018.  "Pluck or Luck: Does Trait
#' Variation or Chance Drive Variation in Lifetime Reproductive
#' Success?"  The American Naturalist 191(4).  DOI:
#' 10.1086/696125  
#' @return A list containing the following:
#' * probXCondR: a matrix whose i,jth component is the probability
#' that the trait has the jth value given that LRO = i-1
#' * probR: A vector whose jth entry is the unconditional probability
#' that LRO = j-1
#' * probRCondX: A matrix whose i,jth entry is the probability that
#' LRO = i-1 given that the trait takes the jth value.
#' @examples
#' PlistAllTraits = FlistAllTraits = list ()
#' P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' P2 = 0.9*P1
#' PlistAllTraits[[1]] = list (P1, P2)
#' P1 = matrix (c(0, 0.5, 0, 0, 0, 0.2, 0, 0, 0.7), 3, 3)
#' P2 = 0.9*P1
#' PlistAllTraits[[2]] = list (P1, P2)
#' F1 = matrix (0, 3, 3); F1[1,] = 0:2
#' F2 = 0.8*F1
#' FlistAllTraits[[1]] = list (F1, F2)
#' F1 = matrix (0, 3, 3); F1[1,] = 0.9*(0:2)
#' F2 = 1.1*F1
#' FlistAllTraits[[2]] = list (F1, F2)
#' Q = matrix (1/2, 2, 2)
#' c0 = c(1,0,0)
#' traitDist = rep(0.5, 2)
#' out = probTraitCondLRO (PlistAllTraits, FlistAllTraits, Q,
#'       c0, maxClutchSize=10, maxLRO=15, traitDist)
probTraitCondLRO = function (PlistAllTraits, FlistAllTraits, Q,
                             c0, maxClutchSize, maxLRO,
                             traitDist,
                             Fdist="Poisson") {
  numTraits = length (PlistAllTraits)
  mT = maxLRO + 1

  ## u0 is the stationary environmental distribution, given by the
  ## dominant eigenvector of Q
  u0 = eigen(Q)$vectors[,1]
  u0 = u0 / sum(u0)

  ## m0 is the stationary state cross-classified by size and
  ## environment
  m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)

  ## P(LRO | trait value) for each trait value
  probRCondX = matrix (0, mT, numTraits)
  ## P(trait value | LRO)
  probXCondR = matrix (0, mT, numTraits)
  
  for (x in 1:numTraits) {
    cat ("Trait", x, "...\n")
    Plist = PlistAllTraits[[x]]
    Flist = FlistAllTraits[[x]]
    probRCondX[,x] = calcDistLRO (Plist, Flist, Q,
                                  c0, maxClutchSize, maxLRO,
                                  Fdist)
  }
  ## Sanity check
  cat ("Testing probRCondX:", colSums(probRCondX), "should = 1.\n")

  ## Now use Bayes thm to get P(X | R). ################

  ## Get joint probability P(R, X)
  probRAndX = matrix (0, mT, numTraits)
  for (i in 1:numTraits)
    ## joint probability P(R, X) = P(R | X) P(X).  
    probRAndX[,i] = probRCondX[,i] * traitDist[i]

  ## Marginalize over X to get P(R)
  probR = apply (probRAndX, 1, sum)
  cat ("Testing probR:", sum(probR), "should = 1.\n")

  ## Condition on R
  for (x in 1:numTraits)
    probXCondR[,x] = probRAndX[,x]/probR

  return (list(probXCondR=probXCondR,
               probR=probR, probRCondX=probRCondX))
}
