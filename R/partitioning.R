##########################################################################
##
## partitioning: functions for partitioning variance and skewness
##
##########################################################################

#' Partitions variance and skewness of LRO
#' 
#' Partitions variance and skewness of lifetime reproductive output
#' (LRO) for models without trait or environmental variation. 
#'
#' @param P The survival/growth transition matrix: P\[i,j\] is the
#'   probability of transitioning from state j to state i.
#' @param F The fecundity matrix: F\[i,j\] is the expected number of
#'   state i offspring from a state j parent
#' @param c0 A vector specifying the offspring state distribution:
#'   c0\[j\] is the probability that an individual is born in state j
#' @param maxAge The maximum age an individual can attain.  Optional.
#'   The default value is 100.
#' @param esR1 The 1st order reward matrix for the extended state
#'   space (see Details).  Optional.  The default assumes
#'   Poisson-distributed clutch sizes.
#' @param esR2 The 2nd order reward matrix for the extended state
#'   space (see Details).  Optional.  The default assumes
#'   Poisson-distributed clutch sizes.
#' @param esR3 The 3rd order reward matrix for the extended state
#'   space (see Details).  Optional.  The default assumes
#'   Poisson-distributed clutch sizes.
#' @param bsR1 The 1st order reward matrix for the extended state
#'   space that includes the pre-birth "ur-state" (see Details).
#'   Optional.  The default assumes Poisson-distributed clutch sizes.
#' @param bsR2 The 2nd order reward matrix for the extended state
#'   space that includes the pre-birth "ur-state" (see Details).
#'   Optional.  The default assumes Poisson-distributed clutch sizes.
#' @param bsR3 The 3rd order reward matrix for the extended state
#'   space that includes the pre-birth "ur-state" (see Details).
#'   Optional.  The default assumes Poisson-distributed clutch sizes.
#' @param survThreshold The threshold to use in determining lifespan
#'   (see return values).  Optional.  Default value is 0.05.
#' @param Fdist The clutch size distribution.  The recognized options
#'   are "Poisson" and "Bernoulli".  Optional.  The default value is
#'   "Poisson".
#' @details A pre-breeding census is assumed: reproduction happens
#'   before survival and growth.  The details of this calculation,
#'   including the definitions of the extended states, can be found in
#'   Robin E. Snyder and Stephen P. Ellner.  2024.
#'   "To prosper, live long: Understanding the sources of reproductive skew and extreme reproductive success in structured populations."
#'   The American Naturalist 204(2) and its online supplemnt.
#' @return A list containing the following: * birthStateVar: the
#'   contribution to Var(LRO) from birth state luck * survTrajecVar: a
#'   vector whose jth entry contains the contribution to Var(LRO) from
#'   survival trajectory luck at age j-1.  * growthTrajecVar: a vector
#'   whose jth entry contains the contribution to Var(LRO) from growth
#'   trajectory luck at age j-1.  * fecVar: a vector whose jth entry
#'   contains the contribution to Var(LRO) from fecundity luck at age
#'   j-1.  * totVar: the total variance in LRO * birthStateSkewness:
#'   the contribution to LRO skewness from birth state luck *
#'   survTrajecSkewness: a vector whose jth entry contains the
#'   contribution to LRO skewness from survival trajectory luck at age
#'   j-1.  * growthTrajecSkewness: a vector whose jth entry contains
#'   the contribution to LRO skewness from growth trajectory luck at
#'   age j-1.  * fecSkewness: a vector whose jth entry contains the
#'   contribution to LRO skewness from fecundity luck at age j-1.  *
#'   totSkewness: the total skewness in LRO * lifespan: the age by
#'   which a proportion survThreshold of a cohort will be dead
#' @examples
#' P = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' F = matrix (0, 3, 3); F[1,] = 0:2
#' c0 = c(1,0,0)
#' out = partitionVarSkewnessNoEnvVar (P, F, c0)
partitionVarSkewnessNoEnvVar = function (P, F, c0, maxAge=100,
                                             survThreshold=0.05,
                                             Fdist="Poisson",
                                             esR1=NULL, esR2=NULL,
                                             esR3=NULL,
                                             bsR1=NULL, bsR2=NULL,
                                             bsR3=NULL) {
  
  ## tolerance for error checking.  0.001 = 0.1% tolerance
  percentTol = 0.001

  mz = dim(P)[1]
  
  ## Sanity check input
  if (Fdist == "Bernoulli") {
    if (sum(colSums(F) > 1))
        stop("Probability of having an offspring > 1!  Columns of fecundity matrix sum to > 1 but clutch size is Bernoulli-distributed.")
  }

  ## Sanity check input
  if (length(c0) != mz)
    stop ("Length of c0 does not match dimension of P")

  ## Sanity check input
  if (dim(F)[1] != mz)
    stop ("P and F should have the same dimensions.")

  ## Make survival, growth, and fecundity bullet matrices
  Sbullet = diag (colSums (P))
  
  ## Does not work if any stage has zero probability of survival. 
  ## Replaced by code below that seems to deal with that issue. 
  #Sinv = diag (1 / diag(Sbullet))
  #Gbullet = P %*% Sinv
  
  Gbullet = P;  
  for(ic in 1:ncol(Gbullet)){
        if(Sbullet[ic,ic]>0) Gbullet[,ic]=Gbullet[,ic]/Sbullet[ic,ic]
  }

  ## Sanity check
  epsilon = 0.00001
  if (sum(abs(range (Gbullet%*%Sbullet-P))) > epsilon)
    warning ("Gbullet %*% Sbullet does not equal P.")
  
  ## #kids isn't part of the state definition, so Fbullet = I
  Fbullet = diag (mz)

  ## create matrices for extended state space, which I will denote with
  ## the prefix "es".

  esmz = 3*mz
  esF = esP = matrix (0, esmz, esmz)

  esP[1:mz, 2*mz + 1:mz] = Gbullet
  esP[mz + 1:mz, 1:mz] = Fbullet
  esP[2*mz + 1:mz, mz + 1:mz] = Sbullet
  
  esF[mz + 1:mz, 1:mz] = F

  ################################################################
  ## Extended state model that includes ur-stage alpha_z.  Used for
  ## calculating contributions from birth stage.  The "bs" before
  ## variable names stands for "birth state."
  ################################################################

  bsmz = 1 + mz
  ## initial state --- everybody starts off in the ur-state.
  bsC0 = rep(0, bsmz)
  bsC0[1] = 1

  bsF = bsP = matrix (0, bsmz, bsmz)
  ## From ur-state to birth state
  bsP[1 + (1:mz),1] = c0
  ## From normal state to normal state
  bsP[1 + 1:mz, 1 + 1:mz] = P

  ## You only start reproducing once you're past the ur-state.  Note
  ## that all we need are the column sums of bsF.
  bsF[1, 1 + 1:mz] = colSums (F)

  ## set up reward matrices ###############################

  if (is.null (esR1) & is.null (esR2) & is.null (esR3) &
        is.null (bsR1) & is.null (bsR2) & is.null (bsR3)) {
    esR1 = esR2 = esR3 = matrix (0, esmz+1, esmz+1)
    bsR1 = bsR2 = bsR3 = matrix (0, bsmz+1, bsmz+1)

    ## First moment of clutch size (Poisson)
    for (j in 1:esmz) 
      esR1[,j] = sum(esF[,j])
    for (j in 1:bsmz) 
      bsR1[,j] = sum(bsF[,j])

    if (Fdist == "Poisson") {
      ## Second moment of clutch size  
      esR2 = esR1 + esR1^2
      bsR2 = bsR1 + bsR1^2
      ## Third moment of clutch size
      esR3 = esR1+ 3*esR1*esR1 + esR1*esR1*esR1
      bsR3 = bsR1+ 3*bsR1*bsR1 + bsR1*bsR1*bsR1
    } else if (Fdist == "Bernoulli") {
      esR2 = esR3 = esR1
      bsR2 = bsR3 = bsR1
    } else {
      stop("Currently only supports Poisson- and Bernoulli-distributed clutch sizes.")
    }
  } else if (!is.null (esR1) & !is.null (esR2) & !is.null (esR3) &
               !is.null (bsR1) & !is.null (bsR2) & !is.null (bsR3)) {
    ## User-supplied reward matrices are fine too.
  } else {
    stop("Values need to be specified for all or none of esR1, esR2, esR3, bs$1, bsR2, bsR3.")
  }
  
  ## Get moments conditional on init. state #########################

  out = calcMoments (esP, esR1, esR2, esR3)
  esRho1Vec = out$rho1Vec
  esMu2Vec = out$mu2Vec
  esSkewnessVec = out$skewnessVec

  out = calcMoments (bsP, bsR1, bsR2, bsR3)
  bsRho1Vec = out$rho1Vec
  bsMu2Vec = out$mu2Vec
  bsMu3Vec = out$mu3Vec
  bsSkewnessVec = out$skewnessVec

  ## Start calculating luck terms ###################################

  survTrajecSkewness = growthTrajecSkewness = 
    fecSkewness = numeric (maxAge)
  fecUpdateSkewness = survUpdateSkewness =
    growthUpdateSkewness =  numeric (maxAge)
  survUpdateVar = growthUpdateVar = fecUpdateVar =
    numeric (maxAge) 
  survTrajecVar = growthTrajecVar = fecVar =
    numeric (maxAge)
  surv = numeric (maxAge)

  ## for ur-state var.
  ## expectation wrt z of Var(success | z)
  expZVarCondZ = bsMu2Vec %*% bsC0
  ## variance wrt z of Exp(success | z)
  varZExpCondZ = sum(bsC0*bsRho1Vec^2) - (sum(bsC0*bsRho1Vec)^2)
  ## Law of total variance
  urVar = expZVarCondZ + varZExpCondZ
  
  justBirthStateVar = bsMu2Vec %*% bsP %*% bsC0
  birthStateVar = urVar - justBirthStateVar

  ## for ur-state skewness
  expZMu3CondZ = bsMu3Vec %*% bsC0
  ## 3rd moment wrt z of Exp(success | z)
  rho3ZExpCondZ = sum(bsC0*(bsRho1Vec)^3)
  ## 2nd moment wrt z of Exp(success | z)
  rho2ZExpCondZ = sum(bsC0*(bsRho1Vec)^2)
  ## 1st moment wrt z of Exp(success | z)
  expZExpCondZ = sum(bsC0*bsRho1Vec)
  ## 3rd central moment wrt z of Exp(success | z)
  mu3ZExpCondZ = rho3ZExpCondZ -
    3*expZExpCondZ*rho2ZExpCondZ +
    2*expZExpCondZ^3
  ## covariance wrt z of Exp(success | z) and Var(success | z)
  covZExpVarCondZ = sum(bsC0*bsRho1Vec*bsMu2Vec) -
    sum(bsC0*bsRho1Vec)*sum(bsC0*bsMu2Vec)
  ## Law of total cumulance
  urMu3 = expZMu3CondZ + mu3ZExpCondZ +
    3*covZExpVarCondZ
  urSkewness = urMu3 / urVar^1.5
  
  justBirthStateSkewness = bsSkewnessVec %*% bsP %*% bsC0
  birthStateSkewness = urSkewness - justBirthStateSkewness
  
  ## Starting the age loop
  PaC0 = c0  ## P^a c_0
  mzZero = rep(0, mz)

  ## Sanity check
  N = solve(diag(mz) - P)
  rho1Vec = colSums (F %*% N)
  foo1 = esRho1Vec %*% c(c0, mzZero, mzZero)
  foo2 = rho1Vec %*% c0
  if ((foo1 < (1 - percentTol)*foo2) | (foo1 > (1 + percentTol)*foo2))
    warning("Calculation of mean via esRho1Vec isn't the same as calculation of mean via rho1Vec.") 

  fecUpdateSkewness[1] = esSkewnessVec %*% c(mzZero, PaC0, mzZero)
  survUpdateSkewness[1] = esSkewnessVec %*% 
    c(mzZero, mzZero, Sbullet %*% PaC0)
  growthUpdateSkewness[1] = esSkewnessVec %*% 
    c(P %*% PaC0, mzZero, mzZero)

  fecUpdateVar[1] = esMu2Vec %*% c(mzZero, PaC0, mzZero)
  survUpdateVar[1] = esMu2Vec %*% c(mzZero, mzZero, Sbullet %*% PaC0)
  growthUpdateVar[1] = esMu2Vec %*% c(P %*% PaC0, mzZero, mzZero)

  surv[1] = sum (PaC0)
  ## a is actually age + 1, since the first year of life is age 0
  for (a in 2:maxAge) {
    PaC0 = P %*% PaC0

    fecUpdateSkewness[a] = esSkewnessVec %*% c(mzZero, PaC0, mzZero)
    survUpdateSkewness[a] = esSkewnessVec %*% 
      c(mzZero, mzZero, Sbullet %*% PaC0)
    growthUpdateSkewness[a] = esSkewnessVec %*% 
      c(P %*% PaC0, mzZero, mzZero)

    fecUpdateVar[a] = esMu2Vec %*% c(mzZero, PaC0, mzZero)
    survUpdateVar[a] = esMu2Vec %*% c(mzZero, mzZero, Sbullet %*% PaC0)
    growthUpdateVar[a] = esMu2Vec %*% c(P %*% PaC0, mzZero, mzZero)

    ## fecundity skewness and var.
    fecSkewness[a] = growthUpdateSkewness[a-1] - fecUpdateSkewness[a]
    fecVar[a] = growthUpdateVar[a-1] - fecUpdateVar[a]
    
    ## survival trajectory skewness and var.
    survTrajecSkewness[a-1] = fecUpdateSkewness[a-1] -
      survUpdateSkewness[a-1]
    survTrajecVar[a-1] = fecUpdateVar[a-1] - survUpdateVar[a-1]
    
    ## growth trajectory skewness and var.
    growthTrajecSkewness[a-1] = survUpdateSkewness[a-1] -
      growthUpdateSkewness[a-1] 
    growthTrajecVar[a-1] = survUpdateVar[a-1] - growthUpdateVar[a-1]

    surv[a] = sum(PaC0)
  }

  ## By what age is a proportion survThreshold of individuals dead?
  lifespan = which(surv < survThreshold)[1]

  ## And get the first value of fecSkewness and fecVar
  fecSkewness[1] = justBirthStateSkewness - fecUpdateSkewness[1]
  fecVar[1] = justBirthStateVar - fecUpdateVar[1]

  ## sanity checks on variance contributions
  if (min(survTrajecVar) < 0)
    stop("Survival trajectory luck contributions to variance are negative.")
  if (min(growthTrajecVar) < 0)
    stop("Growth trajectory luck contributions to variance are negative.")
  if (min(fecVar) < 0)
    stop("Fecundity luck contributions to variance are negative.")
  if (birthStateVar < 0)
    stop("Birth state contribution to variance is negative.")  

  totSurvTrajecSkewness = sum(survTrajecSkewness)
  totGrowthTrajecSkewness = sum(growthTrajecSkewness)
  totFecSkewness = sum(fecSkewness)
  totSkewness = bsSkewnessVec %*% bsC0
  ## Sanity check
  foo = totSurvTrajecSkewness + totGrowthTrajecSkewness +
    totFecSkewness + birthStateSkewness
  if ((foo < (1 - percentTol)*totSkewness) |
      (foo > (1 + percentTol)*totSkewness))
    warning("Skewness components do not sum to total skewness as calculated by bsSkewnessVec.")

  totSurvTrajecVar = sum(survTrajecVar)
  totGrowthTrajecVar = sum(growthTrajecVar)
  totFecVar = sum(fecVar)
  totVar = bsMu2Vec %*% bsC0
  ## Sanity check
  foo = totSurvTrajecVar + totGrowthTrajecVar +
      totFecVar + birthStateVar
  if ((foo < (1 - percentTol)*totVar) |
      (foo > (1 + percentTol)*totVar))
    warning("Variance components do not sum to total variance as calculated by bsMu2Vec.")

  return (out=list(birthStateVar=birthStateVar,
                   survTrajecVar=survTrajecVar,
                   growthTrajecVar=growthTrajecVar,
                   fecVar=fecVar,
                   totVar=totVar,
                   birthStateSkewness=birthStateSkewness,
                   survTrajecSkewness=survTrajecSkewness,
                   growthTrajecSkewness=growthTrajecSkewness,
                   fecSkewness=fecSkewness,
                   totSkewness=totSkewness,
                   lifespan=lifespan))
}

#' Partitions variance and skewness of LRO
#' 
#' Partitions variance and skewness of lifetime reproductive output
#' (LRO) for models with environmental variation but without trait
#' variation.  
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
#' @param maxAge The maximum age an individual can attain.  Optional.
#'   The default value is 100.
#' @param esR1 The 1st order reward matrix for the extended state
#'   space (see Details).  Optional.  The default assumes
#'   Poisson-distributed clutch sizes.
#' @param esR2 The 2nd order reward matrix for the extended state
#'   space (see Details).  Optional.  The default assumes
#'   Poisson-distributed clutch sizes.
#' @param esR3 The 3rd order reward matrix for the extended state
#'   space (see Details).  Optional.  The default assumes
#'   Poisson-distributed clutch sizes.
#' @param bsR1 The 1st order reward matrix for the extended state
#'   space that includes the pre-birth "ur-state" (see Details).
#'   Optional.  The default assumes Poisson-distributed clutch sizes.
#' @param bsR2 The 2nd order reward matrix for the extended state
#'   space that includes the pre-birth "ur-state" (see Details).
#'   Optional.  The default assumes Poisson-distributed clutch sizes.
#' @param bsR3 The 3rd order reward matrix for the extended state
#'   space that includes the pre-birth "ur-state" (see Details).
#'   Optional.  The default assumes Poisson-distributed clutch sizes.
#' @param survThreshold The threshold to use in determining lifespan
#'   (see return values).  Optional.  Default value is 0.05.
#' @param Fdist The clutch size distribution.  The recognized options
#'   are "Poisson" and "Bernoulli".  Optional.  The default value is
#'   "Poisson".
#' @details A pre-breeding census is assumed: reproduction happens
#'   before survival and growth.  The details of this calculation,
#'   including the definitions of the extended states, can be found in
#'   Robin E. Snyder and Stephen P. Ellner.  2024.
#'   "To prosper, live long: Understanding the sources of reproductive
#'   skew and extreme reproductive success in structured populations." 
#'   The American Naturalist 204(2) and its online supplemnt.
#' @return A list containing the following: * birthStateVar: the
#'   contribution to Var(LRO) from birth state luck * birthEnvVar: the
#'   contribution to Var(LRO) from birth environment luck *
#'   survTrajecVar: a vector whose jth entry contains the contribution
#'   to Var(LRO) from survival trajectory luck at age j-1.  *
#'   growthTrajecVar: a vector whose jth entry contains the
#'   contribution to Var(LRO) from growth trajectory luck at age j-1.
#'   * envTrajecVar: a vector whose jth entry contains the
#'   contribution to Var(LRO) from environment trajectory luck at age
#'   j-1.  * fecVar: a vector whose jth entry contains the
#'   contribution to Var(LRO) from fecundity luck at age j-1.  *
#'   totVar: the total variance in LRO * birthStateSkewness: the
#'   contribution to LRO skewness from birth state luck *
#'   birthEnvSkewness: the contribution to LRO skewness from birth
#'   environment luck * survTrajecSkewness: a vector whose jth entry
#'   contains the contribution to LRO skewness from survival
#'   trajectory luck at age j-1.  * growthTrajecSkewness: a vector
#'   whose jth entry contains the contribution to LRO skewness from
#'   growth trajectory luck at age j-1.  * envTrajecSkewness: a vector
#'   whose jth entry contains the contribution to LRO skewness from
#'   environment trajectory luck at age j-1.  * fecSkewness: a vector
#'   whose jth entry contains the contribution to LRO skewness from
#'   fecundity luck at age j-1.  * totSkewness: the total skewness in
#'   LRO * lifespan: the age by which a proportion survThreshold of a
#'   cohort will be dead
#' @examples
#' P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' P2 = matrix (c(0, 0.2, 0, 0, 0, 0.3, 0, 0, 0.4), 3, 3)
#' F1 = matrix (0, 3, 3); F1[1,] = 0:2
#' F2 = matrix (0, 3, 3); F1[1,] = 0.8*(0:2)
#' Plist = list (P1, P2)
#' Flist = list (F1, F2)
#' Q = matrix (1/2, 2, 2)
#' c0 = c(1,0,0)
#' out = partitionVarSkewnessEnvVar (Plist, Flist, Q, c0)
partitionVarSkewnessEnvVar = function (Plist, Flist, Q, c0,
                                       maxAge=100, 
                                       survThreshold=0.05,
                                       Fdist="Poisson",
                                       esR1=NULL, esR2=NULL, esR3=NULL,
                                       bsR1=NULL, bsR2=NULL, bsR3=NULL)
{
  ##  require (Matrix)
  
  ## tolerance for error checking.  0.001 = 0.1% tolerance
  percentTol = 0.001
  
  mz = dim(Plist[[1]])[1]
  numEnv = dim(Q)[1]
  bigmz = numEnv*mz

  ## Sanity check input
  if (Fdist == "Bernoulli") {
    for (q in 1:numEnv)
      if (sum(colSums(Flist[[q]]) > 1))
        stop("Probability of having an offspring > 1!  Columns of fecundity matrix in environment ",
             q, " sum to > 1 but clutch size is Bernoulli-distributed.")
  }

  ## Sanity check input
  if (length(c0) != mz)
    stop ("Length of c0 does not match dimension of P")

  ## Sanity check input
  if (dim(Flist[[1]])[1] != mz)
    stop ("P and F should have the same dimensions.")

  ## u0 is the stationary environmental distribution, given by the
  ## dominant eigenvector of Q
  u0 = eigen(Q)$vectors[,1]
  u0 = u0 / sum(u0)

  m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)

  ## Define megamatrices M and F
  F = M = matrix (0, bigmz, bigmz)
  for (i in 1:numEnv) {
    for (j in 1:numEnv) {
      M[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Plist[[j]]*Q[i,j]
      F[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Flist[[j]]*Q[i,j]
    }
  }

  ## Make survival and growth matrices
  Slist = Glist = list (numEnv)
  for (q in 1:numEnv) {
    Slist[[q]] = diag (colSums (Plist[[q]]))
    Sinv = diag (1 / diag(Slist[[q]]))
    Glist[[q]] = Plist[[q]] %*% Sinv
  }

  ## Make "bullet" matrices.  bdiag makes sparse matrices, which can't
  ## be pasted into regular matrices, so wrap these in as.matrix().
  Sbullet = as.matrix(bdiag (Slist))
  Gbullet = as.matrix(bdiag (Glist))
  Pbullet = Gbullet %*% Sbullet
  Qbullet = matrix (0, bigmz, bigmz)
  for (i in 1:numEnv) {
    for (j in 1:numEnv) {
      Qbullet[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = diag(mz)*Q[i,j]
    }
  }
  ## This code assumes that #kids is not part of the state definition,
  ## so Fbullet = I
  Fbullet = diag(bigmz)

  ## Sanity check
  epsilon = 0.00001
  if (sum(abs(range(Qbullet %*% Gbullet %*% Sbullet - M))) > epsilon)
    warning("Qbullet %*% Gbullet %*% Sbullet does not equal M.")

  ## create matrices for extended state space, which I will denote with
  ## the prefix "es".

  esbigmz = 4*bigmz
  esF = esM = matrix (0, esbigmz, esbigmz)

  esM[bigmz + 1:bigmz, 1:bigmz] = Fbullet
  esM[2*bigmz + 1:bigmz, bigmz + 1:bigmz] = Sbullet
  esM[3*bigmz + 1:bigmz, 2*bigmz + 1:bigmz] = Gbullet
  esM[1:bigmz, 3*bigmz + 1:bigmz] = Qbullet
  
  esF[bigmz + 1:bigmz, 1:bigmz] = F

  ################################################################
  ## Extended state model that includes ur-stage alpha_z and
  ## ur-environment alpha_q.  Used for calculating contributions from
  ## birth stage and birth environment.  The "bs" before variable names
  ## stands for "birth state."
  ################################################################

  bsbigmz = 1 + mz + bigmz
  ## initial state
  bsM0 = rep(0, bsbigmz)
  bsM0[1] = 1

  bsF = bsM = matrix (0, bsbigmz, bsbigmz)
  ## From ur-stage to birth stage and ur-environment
  ## BUG.  Should c0 be m0?  No, I don't think so.  Should c0 be
  ## passed in?
  bsM[1 + (1:mz),1] = c0
  ## From birth stage and ur-environment to m0
  for (q in 1:numEnv) {
    bsM[1 + q*mz + 1:mz, 1 + 1:mz] = diag (u0[q], mz)
  }
  ## From normal stage x env. state to normal stage x env. state
  bsM[1 + mz + 1:bigmz, 1 + mz + 1:bigmz] = M
  ## You only start reproducing once you're past the ur-stages.  I guess
  ## reproduction should go to the ur-stage x ur-environment state.  Not
  ## that it matters much, since all we need are the column sums of bsF.
  bsF[1, 1 + mz + 1:bigmz] = colSums (F)

  ################################################################
  ## set up reward matrices
  ################################################################
  
  if (is.null (esR1) & is.null (esR2) & is.null (esR3) &
        is.null (bsR1) & is.null (bsR2) & is.null (bsR3)) {
    esR1 = esR2 = esR3 = matrix (0, esbigmz+1, esbigmz+1)
    bsR1 = bsR2 = bsR3 = matrix (0, bsbigmz+1, bsbigmz+1)

    ## First moment of clutch size
    for (j in 1:esbigmz) 
      esR1[,j] = sum(esF[,j])
    for (j in 1:bsbigmz) 
      bsR1[,j] = sum(bsF[,j])

    if (Fdist == "Poisson") {
      ## Second moment of clutch size  
      esR2 = esR1 + esR1^2
      bsR2 = bsR1 + bsR1^2
      ## Third moment of clutch size
      esR3 = esR1+ 3*esR1*esR1 + esR1*esR1*esR1
      bsR3 = bsR1+ 3*bsR1*bsR1 + bsR1*bsR1*bsR1
    } else if (Fdist == "Bernoulli") {
      esR2 = esR3 = esR1
      bsR2 = bsR3 = bsR1
    } else {
      stop("Currently only supports Poisson- and Bernoulli-distributed clutch sizes.")
    }
  } else if (!is.null (esR1) & !is.null (esR2) & !is.null (esR3) &
               !is.null (bsR1) & !is.null (bsR2) & !is.null (bsR3)) {
    ## User-supplied reward matrices are fine too.
  } else {
    stop("Values need to be specified for all or none of esR1, esR2, esR3, bs$1, bsR2, bsR3.")
  }

  message ("Calculating es moments...\n")
  out = calcMoments (esM, esR1, esR2, esR3)
  esRho1Vec = out$rho1Vec
  esMu2Vec = out$mu2Vec
  esSkewnessVec = out$skewnessVec

  message ("Calculating bs moments...\n")
  out = calcMoments (bsM, bsR1, bsR2, bsR3)
  bsRho1Vec = out$rho1Vec
  bsMu2Vec = out$mu2Vec
  bsMu3Vec = out$mu3Vec
  bsSkewnessVec = out$skewnessVec

  ## Start calculating luck terms
  survTrajecSkewness = growthTrajecSkewness = envTrajecSkewness =
    fecSkewness = numeric (maxAge)
  fecUpdateSkewness = survUpdateSkewness =
    growthUpdateSkewness = envUpdateSkewness = numeric (maxAge)
  survUpdateVar = growthUpdateVar = envUpdateVar = fecUpdateVar =
    numeric (maxAge) 
  survTrajecVar = growthTrajecVar = envTrajecVar = fecVar =
    numeric (maxAge)
  surv = numeric (maxAge)

  ## for ur-state var.
  ## expectation wrt z of Var(success | z)
  expZVarCondZ = bsMu2Vec %*% bsM0
  ## variance wrt z of Exp(success | z)
  varZExpCondZ = sum(bsM0*bsRho1Vec^2) - (sum(bsM0*bsRho1Vec)^2)
  ## Law of total variance
  urVar = expZVarCondZ + varZExpCondZ
  
  justBirthStageVar = bsMu2Vec %*% bsM %*% bsM0
  birthStageEnvVar = bsMu2Vec %*% bsM %*% bsM %*% bsM0
  
  birthStateVar = urVar - justBirthStageVar
  birthEnvVar = justBirthStageVar - birthStageEnvVar

  ## for ur-state skewness
  expZMu3CondZ = bsMu3Vec %*% bsM0
  ## 3rd moment wrt z of Exp(success | z)
  rho3ZExpCondZ = sum(bsM0*(bsRho1Vec)^3)
  ## 2nd moment wrt z of Exp(success | z)
  rho2ZExpCondZ = sum(bsM0*(bsRho1Vec)^2)
  ## 1st moment wrt z of Exp(success | z)
  expZExpCondZ = sum(bsM0*bsRho1Vec)
  ## 3rd central moment wrt z of Exp(success | z)
  mu3ZExpCondZ = rho3ZExpCondZ -
    3*expZExpCondZ*rho2ZExpCondZ +
    2*expZExpCondZ^3
  ## covariance wrt z of Exp(success | z) and Var(success | z)
  covZExpVarCondZ = sum(bsM0*bsRho1Vec*bsMu2Vec) -
    sum(bsM0*bsRho1Vec)*sum(bsM0*bsMu2Vec)
  ## Law of total cumulance
  urMu3 = expZMu3CondZ + mu3ZExpCondZ +
    3*covZExpVarCondZ
  urSkewness = urMu3 / urVar^1.5
  
  justBirthStageSkewness = bsSkewnessVec %*% bsM %*% bsM0
  birthStageEnvSkewness = bsSkewnessVec %*% bsM %*% bsM %*% bsM0

  birthStateSkewness = urSkewness - justBirthStageSkewness
  birthEnvSkewness = justBirthStageSkewness - birthStageEnvSkewness

  ## Starting the age loop
  MaM0 = m0
  bigmzZero = rep(0, bigmz)

  ## Sanity check --- passes
  N = solve(diag(bigmz) - M)
  rho1Vec = colSums (F %*% N)
  foo1 = esRho1Vec %*% c(bigmzZero, bigmzZero, bigmzZero, m0)
  foo2 = rho1Vec %*% m0
  if ((foo1 < (1 - percentTol)*foo2) |
      (foo1 > (1 + percentTol)*foo2))
    warning("Calculating the mean via esRho1Vec doesn't get the same answer as calculating it via rho1Vec.")

  fecUpdateSkewness[1] = esSkewnessVec %*%
    c(bigmzZero, Fbullet %*% MaM0, bigmzZero, bigmzZero)
  survUpdateSkewness[1] = esSkewnessVec %*% 
    c(bigmzZero, bigmzZero, Sbullet %*% Fbullet %*% MaM0, bigmzZero)
  growthUpdateSkewness[1] = esSkewnessVec %*% 
    c(bigmzZero, bigmzZero, bigmzZero, Pbullet %*% MaM0)
  envUpdateSkewness[1] = esSkewnessVec %*%
    c(M %*% MaM0, bigmzZero, bigmzZero, bigmzZero)

  fecUpdateVar[1] = esMu2Vec %*%
    c(bigmzZero, Fbullet %*% MaM0, bigmzZero, bigmzZero)
  survUpdateVar[1] = esMu2Vec %*%
    c(bigmzZero, bigmzZero, Sbullet %*% MaM0, bigmzZero) 
  growthUpdateVar[1] = esMu2Vec %*%
    c(bigmzZero, bigmzZero, bigmzZero, Pbullet %*% MaM0)
  envUpdateVar[1] = esMu2Vec %*%
    c(M %*% MaM0, bigmzZero, bigmzZero, bigmzZero)

  surv[1] = sum (MaM0)

  for (a in 2:maxAge) {
    MaM0 = M %*% MaM0

    fecUpdateSkewness[a] = esSkewnessVec %*%
      c(bigmzZero, Fbullet %*% MaM0, bigmzZero, bigmzZero)
    survUpdateSkewness[a] = esSkewnessVec %*% 
      c(bigmzZero, bigmzZero, Sbullet %*% Fbullet %*% MaM0, bigmzZero)
    growthUpdateSkewness[a] = esSkewnessVec %*% 
      c(bigmzZero, bigmzZero, bigmzZero, Pbullet %*% MaM0)
    envUpdateSkewness[a] = esSkewnessVec %*%
      c(M %*% MaM0, bigmzZero, bigmzZero, bigmzZero)

    fecUpdateVar[a] = esMu2Vec %*%
      c(bigmzZero, Fbullet %*% MaM0, bigmzZero, bigmzZero)
    survUpdateVar[a] = esMu2Vec %*%
      c(bigmzZero, bigmzZero, Sbullet %*% MaM0, bigmzZero) 
    growthUpdateVar[a] = esMu2Vec %*%
      c(bigmzZero, bigmzZero, bigmzZero, Pbullet %*% MaM0)
    envUpdateVar[a] = esMu2Vec %*%
      c(M %*% MaM0, bigmzZero, bigmzZero, bigmzZero)

    ## fecundity skewness and var.
    fecSkewness[a] = envUpdateSkewness[a-1] - fecUpdateSkewness[a]
    fecVar[a] = envUpdateVar[a-1] - fecUpdateVar[a]

    ## survival trajectory skewness and var.
    survTrajecSkewness[a-1] = fecUpdateSkewness[a-1] -
      survUpdateSkewness[a-1]
    survTrajecVar[a-1] = fecUpdateVar[a-1] - survUpdateVar[a-1]

    ## growth trajectory skewness and var.
    growthTrajecSkewness[a-1] = survUpdateSkewness[a-1] -
      growthUpdateSkewness[a-1] 
    growthTrajecVar[a-1] = survUpdateVar[a-1] - growthUpdateVar[a-1]
    
    ## environment trajectory skewness
    envTrajecSkewness[a-1] = growthUpdateSkewness[a-1] -
      envUpdateSkewness[a-1]
    envTrajecVar[a-1] = growthUpdateVar[a-1] - envUpdateVar[a-1]

    surv[a] = sum (MaM0)
  }  

  ## And get the first value of fecSkewness and fecVar
  fecSkewness[1] = birthStageEnvSkewness - fecUpdateSkewness[1]
  fecVar[1] = birthStageEnvVar - fecUpdateVar[1]

  ## sanity checks on variance contributions
  if (min(survTrajecVar) < 0)
    stop("Survival trajectory luck contributions to variance are negative.")
  if (min(growthTrajecVar) < 0)
    stop("Growth trajectory luck contributions to variance are negative.")
  if (min(fecVar) < 0)
    stop("Fecundity luck contributions to variance are negative.")
  if (min(envTrajecVar) < 0)
    stop("Environment trajectory luck contributions to variance are negative.")
  if (birthStateVar < 0)
    stop("Birth state contribution to variance is negative.")
  if (birthEnvVar < 0)
    stop("Birth environment contribution to variance is negative.")

  ## By what age are survThreshold propor. of individuals dead?
  lifespan = which(surv < survThreshold)[1]

  totSurvTrajecSkewness = sum(survTrajecSkewness)
  totGrowthTrajecSkewness = sum(growthTrajecSkewness)
  totEnvTrajecSkewness = sum(envTrajecSkewness)
  totFecSkewness = sum(fecSkewness)
  totSkewness = bsSkewnessVec %*% bsM0

  ## sanity check
  foo = totSurvTrajecSkewness + totGrowthTrajecSkewness +
         totEnvTrajecSkewness + totFecSkewness +
         birthStateSkewness + birthEnvSkewness
  if ((foo < (1 - percentTol)*totSkewness) |
      (foo > (1 + percentTol)*totSkewness))
    warning("Skewness luck components do not sum to total skewness as calculated by bsSkewnessVec.")
  
  totSurvTrajecVar = sum(survTrajecVar)
  totGrowthTrajecVar = sum(growthTrajecVar)
  totEnvTrajecVar = sum(envTrajecVar)
  totFecVar = sum(fecVar)
  totVar = bsMu2Vec %*% bsM0
  ## Sanity check
  foo = totSurvTrajecVar +
         totGrowthTrajecVar +
         totEnvTrajecVar + totFecVar +
         birthStateVar + birthEnvVar
  if ((foo < (1 - percentTol)*totVar) |
      (foo > (1 + percentTol)*totVar))
    warning("Variance luck components do not sum to total variance as calculated by bsMu2Vec.")

  return (out=list(birthStateVar=birthStateVar,
                   birthEnvVar=birthEnvVar,
                   survTrajecVar=survTrajecVar,
                   growthTrajecVar=growthTrajecVar,
                   envTrajecVar=envTrajecVar,
                   fecVar=fecVar,
                   totVar=totVar,
                   birthStateSkewness=birthStateSkewness,
                   birthEnvSkewness=birthEnvSkewness,
                   survTrajecSkewness=survTrajecSkewness,
                   growthTrajecSkewness=growthTrajecSkewness,
                   envTrajecSkewness=envTrajecSkewness,
                   fecSkewness=fecSkewness,
                   totSkewness=totSkewness,
                   lifespan=lifespan))
}


## Partitions variance and skewness of LRO in the presence of both
## env. variation and trait variation.  Includes the contribution from
## trait variation.

#' Partitions variance and skewness of LRO
#' 
#' Partitions variance and skewness of lifetime reproductive output
#' (LRO) for models with environmental variation and trait
#' variation.  
#'
#' @param PlistAllTraits A list of lists of survival/growth transition
#'   matrices. Plist\[\[x\]\]\[\[q\]\]\[i,j\] is the probability of
#'   transitioning from state j to state i in environment q when an
#'   individual has trait x.
#' @param FlistAllTraits A list of lists of fecundity matrices.
#'   Flist\[\[x\]\]\[\[q\]\]\[i,j\] is the expected number of state i
#'   offspring from a state j, trait x parent in environment q.
#' @param Q The environment transition matrix.  Q\[i,j\] is the
#'   probability of transitioning from environment j to environment i.
#' @param c0 A vector specifying the offspring state distribution:
#'   c0\[j\] is the probability that an individual is born in state j
#' @param traitDist A vector whose jth entry is the probability of
#'   having trait j
#' @param maxAge The maximum age an individual can attain.  Optional.
#'   The default value is 100.
#' @param survThreshold The threshold to use in determining lifespan
#'   (see return values).  Optional.  Default value is 0.05.
#' @param Fdist The clutch size distribution.  The recognized options
#'   are "Poisson" and "Bernoulli".  Optional.  The default value is
#'   "Poisson".
#' @details A pre-breeding census is assumed: reproduction happens
#'   before survival and growth.  The details of this calculation,
#'   including the definitions of the extended states, can be found in
#'   Robin E. Snyder and Stephen P. Ellner.  2024.
#'   "To prosper, live long: Understanding the sources of reproductive
#'   skew and extreme reproductive success in structured populations." 
#'   The American Naturalist 204(2) and its online supplemnt.
#' @return A list containing the following: * varFromTraits = the
#'   contribution to Var(LRO) from trait variation *
#'   skewnessFromTraits = the contribution to LRO skewness from trait
#'   variation * totVar: the total variance in LRO * totSkewness: the
#'   total skewness in LRO * birthEnvVar: the contribution to Var(LRO)
#'   from birth environment luck * birthEnvSkewness: the contribution
#'   to LRO skewness from birth environment luck * birthStateVar: the
#'   contribution to Var(LRO) from birth state luck *
#'   birthStateSkewness: the contribution to LRO skewness from birth
#'   state luck * survTrajecVar: a vector whose jth entry contains the
#'   contribution to Var(LRO) from survival trajectory luck at age
#'   j-1.  * survTrajecSkewness: a vector whose jth entry contains the
#'   contribution to LRO skewness from survival trajectory luck at age
#'   j-1.  * growthTrajecVar: a vector whose jth entry contains the
#'   contribution to Var(LRO) from growth trajectory luck at age j-1.
#'   * growthTrajecSkewness: a vector whose jth entry contains the
#'   contribution to LRO skewness from growth trajectory luck at age
#'   j-1.  * envTrajecVar: a vector whose jth entry contains the
#'   contribution to Var(LRO) from environment trajectory luck at age
#'   j-1.  * envTrajecSkewness: a vector whose jth entry contains the
#'   contribution to LRO skewness from environment trajectory luck at
#'   age j-1.  * fecVar: a vector whose jth entry contains the
#'   contribution to Var(LRO) from fecundity luck at age j-1.  *
#'   fecSkewness: a vector whose jth entry contains the contribution
#'   to LRO skewness from fecundity luck at age j-1.  * totVarCondX: a
#'   vector whose jth entry is Var(LRO | trait = j) *
#'   totSkewnessCondX: a vector whose jth entry is the LRO skewness
#'   conditional on trait = j * birthEnvVarCondX: a vector whose jth
#'   entry is the contribution of birth environment luck to Var(LRO |
#'   trait = j) * birthEnvSkewnessCondX: a vector whose jth entry is
#'   the contribution of birth environment luck to LRO skewness
#'   conditional on trait = j * birthStateVarCondX: a vector whose jth
#'   entry is the contribution of birth state luck to Var(LRO | trait
#'   = j) * birthStateSkewnessCondX: a vector whose jth entry is the
#'   contribution of birth state luck to LRO skewness conditional on
#'   trait = j * survTrajecVarCondX: a matrix whose i,jth entry
#'   contains the contribution to Var(LRO | trait = i) from survival
#'   trajectory luck at age j-1.  * survTrajecSkewnessCondX: a matrix
#'   whose i,jth entry contains the contribution to LRO skewness
#'   conditional on trait = i from survival trajectory luck at age
#'   j-1.  * growthTrajecVarCondX: a matrix whose i,jth entry contains
#'   the contribution to Var(LRO | trait = i) from growth trajectory
#'   luck at age j-1.  * growthTrajecSkewnessCondX: a matrix whose
#'   i,jth entry contains the contribution to LRO skewness conditional
#'   on trait = i from growth trajectory luck at age j-1.  *
#'   envTrajecVarCondX: a matrix whose i,jth entry contains the
#'   contribution to Var(LRO | trait = i) from environment trajectory
#'   luck at age j-1.  * envTrajecSkewnessCondX: a matrix whose i,jth
#'   entry contains the contribution to LRO skewness conditional on
#'   trait = i from environment trajectory luck at age j-1.  *
#'   fecVarCondX: a matrix whose i,jth entry contains the contribution
#'   to Var(LRO | trait = i) from fecundity luck at age j-1.  *
#'   fecSkewnessCondX: a matrix whose i,jth entry contains the
#'   contribution to LRO skewness conditional on trait = i from
#'   fecundity luck at age j-1.  * lifespanCondX: a vector whose jth
#'   entry is the age by which survThreshold proportion of a cohort
#'   with trait j would be dead
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
#' out = partitionVarSkewnessEnvVarAndTraits (PlistAllTraits,
#'       FlistAllTraits, Q, c0, traitDist)
partitionVarSkewnessEnvVarAndTraits = function (PlistAllTraits,
                                                FlistAllTraits, Q,
                                                c0, traitDist,
                                                maxAge=100,
                                                survThreshold=0.05,
                                                Fdist="Poisson") {

  ## tolerance for error checking.  0.001 = 0.1% tolerance
  percentTol = 0.001

  mz = dim(PlistAllTraits[[1]][[1]])[1]

  ## u0 is the stationary environmental distribution, given by the
  ## dominant eigenvector of Q
  u0 = eigen(Q)$vectors[,1]
  u0 = u0 / sum(u0)

  ## number of traits
  numTraits = length (PlistAllTraits)

  ## number of environment states
  numEnv = dim(Q)[1]

  bigmz = mz*numEnv
  ## Initial cross-classified state
  m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)

  ## Sanity check input
  if (Fdist == "Bernoulli") {
    for (x in 1:numTraits) {
      for (q in 1:numEnv)
        if (sum(colSums(FlistAllTraits[[x]][[q]]) > 1))
          stop("Probability of having an offspring > 1!  Columns of fecundity matrix in environment ",
               q, " with trait ", x," sum to > 1 but clutch size is Bernoulli-distributed.")
    }
  }

  ## Sanity check input
  if (length(c0) != mz)
    stop ("Length of c0 does not match dimension of P")

  ## Sanity check input
  if (dim(FlistAllTraits[[1]][[1]])[1] != mz)
    stop ("P and F should have the same dimensions.")

  lifespanCondX = birthEnvSkewnessCondX =
    birthEnvVarCondX = birthStateSkewnessCondX =
      birthStateVarCondX = numeric (numTraits)
  expRCondX = totVarCondX = totSkewnessCondX = numeric (numTraits)

  survTrajecVarCondX = growthTrajecVarCondX = envTrajecVarCondX = fecVarCondX =
    matrix (0, numTraits, maxAge)
  survTrajecSkewnessCondX = growthTrajecSkewnessCondX = envTrajecSkewnessCondX =
    fecSkewnessCondX = matrix (0, numTraits, maxAge)

  ## For sanity checking
  expRCondXZ = matrix (0, numTraits, bigmz)

  for (x in 1:numTraits) {

    message ("About to define Plist for x =", x, "\n")
    Plist = PlistAllTraits[[x]]
    Flist = FlistAllTraits[[x]]

    out = partitionVarSkewnessEnvVar (Plist=Plist, Flist=Flist,
                                      Q=Q, c0=c0, maxAge=maxAge,
                                      survThreshold=survThreshold,
                                      Fdist=Fdist)

    birthEnvSkewnessCondX[x] = out$birthEnvSkewness
    survTrajecSkewnessCondX[x,] = out$survTrajecSkewness
    growthTrajecSkewnessCondX[x,] = out$growthTrajecSkewness
    envTrajecSkewnessCondX[x,] = out$envTrajecSkewness
    fecSkewnessCondX[x,] = out$fecSkewness
    totSkewnessCondX[x] = out$totSkewness
    surv = out$surv
    lifespanCondX[x] = out$lifespan
    
    birthEnvVarCondX[x] = out$birthEnvVar
    survTrajecVarCondX[x,] = out$survTrajecVar
    growthTrajecVarCondX[x,] = out$growthTrajecVar
    envTrajecVarCondX[x,] = out$envTrajecVar
    fecVarCondX[x,] = out$fecVar
    totVarCondX[x] = out$totVar
    
    ## Get expected LRO while you're at it.
    ## Define megamatrices M and F
    F = M = matrix (0, bigmz, bigmz)
    for (i in 1:numEnv) {
      for (j in 1:numEnv) {
        M[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Plist[[j]]*Q[i,j]
        F[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Flist[[j]]*Q[i,j]
      }
    }
    N = solve (diag(bigmz) - M)
    expRCondXZ[x,] = (colSums (F %*% N))
    expRCondX[x] = (colSums (F %*% N)) %*% m0
  } ## end loop over traits
  
  ## Average over traits
  birthStateSkewness = traitAve (birthStateSkewnessCondX, traitDist)
  birthEnvSkewness = traitAve (birthEnvSkewnessCondX, traitDist)
  fecSkewness = apply (fecSkewnessCondX, 2, traitAve, traitDist)
  survTrajecSkewness = apply (survTrajecSkewnessCondX, 2, traitAve, traitDist)
  growthTrajecSkewness = apply (growthTrajecSkewnessCondX, 2, traitAve, traitDist)
  envTrajecSkewness = apply (envTrajecSkewnessCondX, 2, traitAve, traitDist)

  birthStateVar = traitAve (birthStateVarCondX, traitDist)
  birthEnvVar = traitAve (birthEnvVarCondX, traitDist)
  fecVar = apply (fecVarCondX, 2, traitAve, traitDist)
  survTrajecVar = apply (survTrajecVarCondX, 2, traitAve, traitDist)
  growthTrajecVar = apply (growthTrajecVarCondX, 2, traitAve, traitDist)
  envTrajecVar = apply (envTrajecVarCondX, 2, traitAve, traitDist)

  aveLifespan = traitAve (lifespanCondX, traitDist)

  totSurvTrajecSkewness = sum (survTrajecSkewness)
  totGrowthTrajecSkewness = sum (growthTrajecSkewness)
  totEnvTrajecSkewness = sum (envTrajecSkewness)
  totFecSkewness = sum (fecSkewness)

  totSurvTrajecVar = sum (survTrajecVar)
  totGrowthTrajecVar = sum (growthTrajecVar)
  totEnvTrajecVar = sum (envTrajecVar)
  totFecVar = sum (fecVar)

  ## Use the law of total cumulance to get the full variance and
  ## skewness.

  ## Variance
  expXVarRCondX = traitAve (totVarCondX, traitDist)
  varXExpRCondX = traitVar (expRCondX, traitDist)
  totVar = expXVarRCondX + varXExpRCondX
  
  ## To get total skewness or variance, we need to create transition
  ## matrices and reward matrices for a model in which state is
  ## classified by trait x stage x env.  I'll use the prefix "ts" to
  ## suggest "triple state" as a way to distinguish these matrices.
  ##
  ## And since we need to start off in the ur-state, we need to include
  ## pre-birth states (the ur-trait, ur-stage, and
  ## ur-environment)), so the full prefix is bsts...
  ##
  ## A guide to index ranges:
  ##
  ## 1 = the ur state.
  ##
  ## 2--(numTraits+1): trait is assigned but nothing else.  
  ## 
  ## (numTraits+2)--(numTraits*mz + numTraits+1): trait and z are
  ## assigned, but not q. 
  ##
  ## (numTraits*mz + numTraits+2)--(numTraits*mz + numTraits+ 1 + mz*numEnv):
  ## trait 1, (z,q).  
  ##
  ## Dare I say etc.?  For trait 2, (z, q), trait 3, (z, q), and so on.
  ##
  ## See SI section S1 of Snyder and Ellner 2024 for a guide to these
  ## calculations.

  message ("Calculating the effect of trait variation...\n")

  tsbigmz = numTraits * bigmz
  bstsbigmz = 1 + numTraits + numTraits*mz + tsbigmz

  Mlist = bigFlist = list()
  for (x in 1:numTraits) {
    Mlist[[x]] = bigFlist[[x]] = matrix (NA, bigmz, bigmz)
    for (i in 1:numEnv) {
      for (j in 1:numEnv) {
        Mlist[[x]][(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] =
          PlistAllTraits[[x]][[j]]*Q[i,j]
        bigFlist[[x]][(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] =
          FlistAllTraits[[x]][[j]]*Q[i,j]
      }
    }
  }

  tsM = as.matrix(bdiag (Mlist))
  tsF = as.matrix(bdiag (bigFlist))
  rm (Mlist, bigFlist)
  
  bstsF = bstsM = matrix (0, bstsbigmz, bstsbigmz)
  ## From ur-stage to just traits (submatrix W_{21} in SI)
  bstsM[1 + (1:numTraits),1] = traitDist
  ## From just traits to traits and stage (submatrix W_{32} in SI)
  for (x in 1:numTraits) 
    bstsM[1 + numTraits + (x-1)*mz + 1:mz,1 + x] = c0
  ## From traits and stage to traits and stage and environment (m0)
  ## (submatrix W_{43} in SI)
  for (x in 1:numTraits) {
    for (q in 1:numEnv) 
      bstsM[1 + numTraits + numTraits*mz + (x-1)*bigmz + (q-1)*mz + 1:mz,
            1 + numTraits + (x-1)*mz + 1:mz] =
        diag (u0[q], mz)
  }
  ## From normal stage x env. state to normal stage x env. state
  ## (submatrix W_{44} in SI)
  bstsM[1 + numTraits + numTraits*mz + 1:tsbigmz,
        1 + numTraits + numTraits*mz + 1:tsbigmz] = tsM

  ## You only start reproducing once you're past the ur-stages.  
  ## Reproduction should go to the ur-stage x ur-environment state.  Not
  ## that it matters much, since all we need are the column sums of
  ## bsF.
  bstsF[1, 1 + numTraits + numTraits*mz + 1:tsbigmz] =
      colSums (tsF)

  ## The following reward matrix calculations come from van Daalen,
  ## S. F., and H. Caswell. 2017. Lifetime reproductive output:
  ## individual stochasticity, variance, and sensitivity
  ## analysis. Theoretical Ecology 10:355 – 374.

  ## set up bstsMplus, the transition matrix with absorbing state (dead)
  bstsMplus = matrix (0, bstsbigmz+1, bstsbigmz+1)
  bstsMplus[1:bstsbigmz, 1:bstsbigmz] = bstsM
  bstsMplus[bstsbigmz+1, 1:bstsbigmz] = 1 - colSums(bstsM)
  bstsMplus[bstsbigmz+1, bstsbigmz+1] = 1
  ## Set up bstsZ, which cleaves off dead states
  bstsZ = cbind (diag(bstsbigmz), rep(0, bstsbigmz))

  ## Set up the vector 1_s, where s = bigmz + 1
  bstsOneS = rep(0, bstsbigmz+1)
  bstsOneS[bstsbigmz+1] = 1

  ## set up reward matrices
  bstsR1 = bstsR2 = bstsR3 = matrix (0, bstsbigmz+1, bstsbigmz+1)
  bstsR1Stripped = bstsR2Stripped = matrix (0, bstsbigmz, bstsbigmz)

  ## First moment of clutch size.
  foo = colSums(bstsF)
  for (j in 1:bstsbigmz)
    bstsR1[,j] = foo[j]

  if (Fdist == "Poisson") {
    ## bstsR2 = bstsR1^2
    bstsR2 = bstsR1 + bstsR1^2
    ## bstsR3 = bstsR1^3
    bstsR3 = bstsR1 + 3*bstsR1^2 + bstsR1^3
  } else if (Fdist == "Bernoulli") {
    bstsR2 = bstsR3 = bstsR1
  } else {
    stop("Currently only supports Poisson- and Bernoulli-distributed clutch sizes.")
  }
  
  ## R1Stripped has the absorbed state cleaved off
  bstsR1Stripped = bstsR1[1:bstsbigmz,1:bstsbigmz]
  bstsR2Stripped = bstsR2[1:bstsbigmz,1:bstsbigmz]

  ## Ex(R)
  message ("Calculating bstsRho1Vec...\n")
  bstse = rep (1, bstsbigmz+1)
  bstsN = solve (diag(bstsbigmz) - bstsM)

  bstsRho1Vec = t(bstsN) %*% bstsZ %*% t(bstsMplus * bstsR1) %*% bstse
  ## Sanity check
  epsilon = 0.00001
  foo = bstsRho1Vec[(numTraits*mz + numTraits+2):
                    (numTraits*mz + numTraits+ 1 + mz*numEnv)]
  if (sum(abs(range(foo - expRCondXZ[1,]))) > epsilon)
    warning("bstsRho1Vec differs substantially from expRCondXZ[1,].")
        
  message ("Calculating bstsRho2Vec...\n")
  bstsRho2Vec = t(bstsN) %*% (bstsZ %*% t(bstsMplus * bstsR2) %*% bstse +
                              2*t(bstsM * bstsR1Stripped) %*%
                              bstsRho1Vec)
      
  ## Ex(R^3 | x)
  message ("Calculating bstsRho3Vec...\n")
  bstsRho3Vec = t(bstsN) %*% (bstsZ %*% t(bstsMplus * bstsR3) %*% bstse +
                              3*t(bstsM*bstsR2Stripped) %*% bstsRho1Vec +
                              3*t(bstsM*bstsR1Stripped) %*% bstsRho2Vec)

  bstsRho1Vec = t(bstsRho1Vec)
  bstsRho2Vec = t(bstsRho2Vec)
  bstsRho3Vec = t(bstsRho3Vec)

  ## Do I get the correct bstsMu2Vec for trait 3, all values of (z,q)?
  bstsMu2Vec = bstsRho2Vec - bstsRho1Vec^2
  bstsMu3Vec = bstsRho3Vec - 3*bstsRho1Vec*bstsRho2Vec + 2*bstsRho1Vec^3
  bstsSkewnessVec = bstsMu3Vec / bstsMu2Vec^1.5

  urState = rep(0, bstsbigmz)
  urState[1] = 1
  totSkewness = bstsSkewnessVec %*% urState
  totVar2 = bstsMu2Vec %*% urState

  if ((totVar < (1 - percentTol)*totVar2) |
      (totVar > (1 + percentTol)*totVar2))
    warning("Variance calculated from bstsM differs from variance calculated from variances for each trait value.")

  skewnessFromTraits = totSkewness - traitAve (totSkewnessCondX, traitDist)
  varFromTraits = totVar - traitAve (totVarCondX, traitDist)

  ## sanity check
  if (varFromTraits < 0)
    stop ("Trait contribution to variance is negative")

  ## Sanity check
  if ((varFromTraits < (1 - percentTol)*varXExpRCondX) |
      (varFromTraits > (1 + percentTol)*varXExpRCondX))
    warning("The portion of variance due to traits that is calculated by subtracting Ex_x Var(R | x) from total variance is not equal to Var_x Ex(R | x).")

  return (out=list(varFromTraits=varFromTraits,
                   skewnessFromTraits=skewnessFromTraits,
                   totVar=totVar, totSkewness=totSkewness,
                   birthEnvVar=birthEnvVar,
                   birthEnvSkewness=birthEnvSkewness,
                   birthStateVar=birthStateVar,
                   birthStateSkewness=birthStateSkewness,
                   survTrajecVar=survTrajecVar,
                   survTrajecSkewness=survTrajecSkewness,
                   growthTrajecVar=growthTrajecVar,
                   growthTrajecSkewness=growthTrajecSkewness,
                   envTrajecVar=envTrajecVar,
                   envTrajecSkewness=envTrajecSkewness,
                   fecVar=fecVar,
                   fecSkewness=fecSkewness,
                   totVarCondX=totVarCondX,
                   totSkewnessCondX=totSkewnessCondX,
                   birthEnvVarCondX=birthEnvVarCondX,
                   birthEnvSkewnessCondX=birthEnvSkewnessCondX,
                   birthStateVarCondX=birthStateVarCondX,
                   birthStateSkewnessCondX=birthStateSkewnessCondX,
                   survTrajecVarCondX=survTrajecVarCondX,
                   survTrajecSkewnessCondX=survTrajecSkewnessCondX,
                   growthTrajecVarCondX=growthTrajecVarCondX,
                   growthTrajecSkewnessCondX=growthTrajecSkewnessCondX,
                   envTrajecVarCondX=envTrajecVarCondX,
                   envTrajecSkewnessCondX=envTrajecSkewnessCondX,
                   fecVarCondX=fecVarCondX,
                   fecSkewnessCondX=fecSkewnessCondX,
                   lifespanCondX=lifespanCondX))
}

############### Post-breeding census versions ###################

#' Partitions variance and skewness of LRO for a post-breeding census
#' 
#' Partitions variance and skewness of lifetime reproductive output
#' (LRO) for models without trait or environmental variation, assuming
#' a post-breeding census.
#'
#' @param P The survival/growth transition matrix: P\[i,j\] is the
#'   probability of transitioning from state j to state i.
#' @param F The fecundity matrix: F\[i,j\] is the expected number of
#'   state i offspring from a state j parent
#' @param c0 A vector specifying the offspring state distribution:
#'   c0\[j\] is the probability that an individual is born in state j
#' @param maxAge The maximum age an individual can attain.  Optional.
#'   The default value is 100.
#' @param esR1 The 1st order reward matrix for the extended state
#'   space (see Details).  Optional.  The default assumes
#'   Poisson-distributed clutch sizes.
#' @param esR2 The 2nd order reward matrix for the extended state
#'   space (see Details).  Optional.  The default assumes
#'   Poisson-distributed clutch sizes.
#' @param esR3 The 3rd order reward matrix for the extended state
#'   space (see Details).  Optional.  The default assumes
#'   Poisson-distributed clutch sizes.
#' @param bsR1 The 1st order reward matrix for the extended state
#'   space that includes the pre-birth "ur-state" (see Details).
#'   Optional.  The default assumes Poisson-distributed clutch sizes.
#' @param bsR2 The 2nd order reward matrix for the extended state
#'   space that includes the pre-birth "ur-state" (see Details).
#'   Optional.  The default assumes Poisson-distributed clutch sizes.
#' @param bsR3 The 3rd order reward matrix for the extended state
#'   space that includes the pre-birth "ur-state" (see Details).
#'   Optional.  The default assumes Poisson-distributed clutch sizes.
#' @param survThreshold The threshold to use in determining lifespan
#'   (see return values).  Optional.  Default value is 0.05.
#' @param Fdist The clutch size distribution.  The recognized options
#'   are "Poisson" and "Bernoulli".  Optional.  The default value is
#'   "Poisson".
#' @details A post-breeding census is assumed: reproduction happens
#'   after survival and growth.  The details of this calculation,
#'   including the definitions of the extended states, can be found in
#'   Robin E. Snyder and Stephen P. Ellner.  2024.
#'   "To prosper, live long: Understanding the sources of reproductive skew and extreme reproductive success in structured populations."
#'   The American Naturalist 204(2) and its online supplemnt.
#' @return A list containing the following: * birthStateVar: the
#'   contribution to Var(LRO) from birth state luck * survTrajecVar: a
#'   vector whose jth entry contains the contribution to Var(LRO) from
#'   survival trajectory luck at age j-1.  * growthTrajecVar: a vector
#'   whose jth entry contains the contribution to Var(LRO) from growth
#'   trajectory luck at age j-1.  * fecVar: a vector whose jth entry
#'   contains the contribution to Var(LRO) from fecundity luck at age
#'   j-1.  * totVar: the total variance in LRO * birthStateSkewness:
#'   the contribution to LRO skewness from birth state luck *
#'   survTrajecSkewness: a vector whose jth entry contains the
#'   contribution to LRO skewness from survival trajectory luck at age
#'   j-1.  * growthTrajecSkewness: a vector whose jth entry contains
#'   the contribution to LRO skewness from growth trajectory luck at
#'   age j-1.  * fecSkewness: a vector whose jth entry contains the
#'   contribution to LRO skewness from fecundity luck at age j-1.  *
#'   totSkewness: the total skewness in LRO * lifespan: the age by
#'   which a proportion survThreshold of a cohort will be dead
#' @examples
#' P = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' F = matrix (0, 3, 3); F[1,] = 0:2
#' c0 = c(1,0,0)
#' out = partitionVarSkewnessNoEnvVarPostBreeding (P, F, c0)
partitionVarSkewnessNoEnvVarPostBreeding = function (P, F, c0, maxAge=100,
            survThreshold=0.05,
            Fdist="Poisson",
            esR1=NULL, esR2=NULL,
            esR3=NULL,
            bsR1=NULL, bsR2=NULL,
            bsR3=NULL) {
    
  ## tolerance for error checking.  0.001 = 0.1% tolerance
  percentTol = 0.001

  mz = dim(P)[1]
  
  ## Sanity check input
  if (Fdist == "Bernoulli") {
    if (sum(colSums(F) > 1))
        stop("Probability of having an offspring > 1!  Columns of fecundity matrix sum to > 1 but clutch size is Bernoulli-distributed.")
  }

  ## Sanity check input
  if (length(c0) != mz)
    stop ("Length of c0 does not match dimension of P")

  ## Sanity check input
  if (dim(F)[1] != mz)
    stop ("P and F should have the same dimensions.")

  ## Make survival, growth, and fecundity bullet matrices
  Sbullet = diag (colSums (P))
  
  ## Does not work if any stage has zero probability of survival. 
  ## Replaced by code below that seems to deal with that issue. 
  #Sinv = diag (1 / diag(Sbullet))
  #Gbullet = P %*% Sinv
  
  Gbullet = P;  
  for(ic in 1:ncol(Gbullet)){
        if(Sbullet[ic,ic]>0) Gbullet[,ic]=Gbullet[,ic]/Sbullet[ic,ic]
  }

  ## Sanity check
  epsilon = 0.00001
  if (sum(abs(range (Gbullet%*%Sbullet-P))) > epsilon)
    warning ("Gbullet %*% Sbullet does not equal P.")
  
  ## #kids isn't part of the state definition, so Fbullet = I
  Fbullet = diag (mz)

  ## create matrices for extended state space, which I will denote with
  ## the prefix "es".

  esmz = 3*mz
  esF = esP = matrix (0, esmz, esmz)

  esP[1:mz, 2*mz + 1:mz] = Fbullet
  esP[mz + 1:mz, 1:mz] = Sbullet
  esP[2*mz + 1:mz, mz + 1:mz] = Gbullet
  
  esF[1:mz, 2*mz + 1:mz] = F

  ################################################################
  ## Extended state model that includes ur-stage alpha_z.  Used for
  ## calculating contributions from birth stage.  The "bs" before
  ## variable names stands for "birth state."
  ################################################################

  bsmz = 1 + mz
  ## initial state --- everybody starts off in the ur-state.
  bsC0 = rep(0, bsmz)
  bsC0[1] = 1

  bsF = bsP = matrix (0, bsmz, bsmz)
  ## From ur-state to birth state
  bsP[1 + (1:mz),1] = c0
  ## From normal state to normal state
  bsP[1 + 1:mz, 1 + 1:mz] = P

  ## You only start reproducing once you're past the ur-state.  Note
  ## that all we need are the column sums of bsF.
  bsF[1, 1 + 1:mz] = colSums (F)

  ## set up reward matrices ###############################

  if (is.null (esR1) & is.null (esR2) & is.null (esR3) &
        is.null (bsR1) & is.null (bsR2) & is.null (bsR3)) {
    esR1 = esR2 = esR3 = matrix (0, esmz+1, esmz+1)
    bsR1 = bsR2 = bsR3 = matrix (0, bsmz+1, bsmz+1)

    ## First moment of clutch size (Poisson)
    for (j in 1:esmz) 
      esR1[,j] = sum(esF[,j])
    for (j in 1:bsmz) 
      bsR1[,j] = sum(bsF[,j])

    if (Fdist == "Poisson") {
      ## Second moment of clutch size  
      esR2 = esR1 + esR1^2
      bsR2 = bsR1 + bsR1^2
      ## Third moment of clutch size
      esR3 = esR1+ 3*esR1*esR1 + esR1*esR1*esR1
      bsR3 = bsR1+ 3*bsR1*bsR1 + bsR1*bsR1*bsR1
    } else if (Fdist == "Bernoulli") {
      esR2 = esR3 = esR1
      bsR2 = bsR3 = bsR1
    } else {
      stop("Currently only supports Poisson- and Bernoulli-distributed clutch sizes.")
    }
  } else if (!is.null (esR1) & !is.null (esR2) & !is.null (esR3) &
               !is.null (bsR1) & !is.null (bsR2) & !is.null (bsR3)) {
    ## User-supplied reward matrices are fine too.
  } else {
    stop("Values need to be specified for all or none of esR1, esR2, esR3, bs$1, bsR2, bsR3.")
  }
  
  ## Get moments conditional on init. state #########################

  out = calcMoments (esP, esR1, esR2, esR3)
  esRho1Vec = out$rho1Vec
  esMu2Vec = out$mu2Vec
  esSkewnessVec = out$skewnessVec

  out = calcMoments (bsP, bsR1, bsR2, bsR3)
  bsRho1Vec = out$rho1Vec
  bsMu2Vec = out$mu2Vec
  bsMu3Vec = out$mu3Vec
  bsSkewnessVec = out$skewnessVec

  ## Start calculating luck terms ###################################

  survTrajecSkewness = growthTrajecSkewness = 
    fecSkewness = numeric (maxAge)
  fecUpdateSkewness = survUpdateSkewness =
    growthUpdateSkewness =  numeric (maxAge)
  survUpdateVar = growthUpdateVar = fecUpdateVar =
    numeric (maxAge) 
  survTrajecVar = growthTrajecVar = fecVar =
    numeric (maxAge)
  surv = numeric (maxAge)

  ## for ur-state var.
  ## expectation wrt z of Var(success | z)
  expZVarCondZ = bsMu2Vec %*% bsC0
  ## variance wrt z of Exp(success | z)
  varZExpCondZ = sum(bsC0*bsRho1Vec^2) - (sum(bsC0*bsRho1Vec)^2)
  ## Law of total variance
  urVar = expZVarCondZ + varZExpCondZ
  
  justBirthStateVar = bsMu2Vec %*% bsP %*% bsC0
  birthStateVar = urVar - justBirthStateVar

  ## for ur-state skewness
  expZMu3CondZ = bsMu3Vec %*% bsC0
  ## 3rd moment wrt z of Exp(success | z)
  rho3ZExpCondZ = sum(bsC0*(bsRho1Vec)^3)
  ## 2nd moment wrt z of Exp(success | z)
  rho2ZExpCondZ = sum(bsC0*(bsRho1Vec)^2)
  ## 1st moment wrt z of Exp(success | z)
  expZExpCondZ = sum(bsC0*bsRho1Vec)
  ## 3rd central moment wrt z of Exp(success | z)
  mu3ZExpCondZ = rho3ZExpCondZ -
    3*expZExpCondZ*rho2ZExpCondZ +
    2*expZExpCondZ^3
  ## covariance wrt z of Exp(success | z) and Var(success | z)
  covZExpVarCondZ = sum(bsC0*bsRho1Vec*bsMu2Vec) -
    sum(bsC0*bsRho1Vec)*sum(bsC0*bsMu2Vec)
  ## Law of total cumulance
  urMu3 = expZMu3CondZ + mu3ZExpCondZ +
    3*covZExpVarCondZ
  urSkewness = urMu3 / urVar^1.5
  
  justBirthStateSkewness = bsSkewnessVec %*% bsP %*% bsC0
  birthStateSkewness = urSkewness - justBirthStateSkewness
  
  ## Starting the age loop
  PaC0 = c0  ## P^a c_0
  mzZero = rep(0, mz)

  ## Sanity check
  N = solve(diag(mz) - P)
  rho1Vec = colSums (F %*% N)
  foo1 = esRho1Vec %*% c(c0, mzZero, mzZero)
  foo2 = rho1Vec %*% c0
  if ((foo1 < (1 - percentTol)*foo2) | (foo1 > (1 + percentTol)*foo2))
    warning("Calculation of mean via esRho1Vec isn't the same as calculation of mean via rho1Vec.")

  survUpdateSkewness[1] = esSkewnessVec %*%
    c(mzZero, Sbullet %*% PaC0, mzZero)  
  growthUpdateSkewness[1] = esSkewnessVec %*%
    c(mzZero, mzZero, P %*% PaC0) 
  fecUpdateSkewness[1] = esSkewnessVec %*%
    c(P %*% PaC0, mzZero, mzZero)

  survUpdateVar[1] = esMu2Vec %*% c(mzZero, Sbullet %*% PaC0, mzZero)  
  growthUpdateVar[1] = esMu2Vec %*% c(mzZero, mzZero, P %*% PaC0) 
  fecUpdateVar[1] = esMu2Vec %*% c(P %*% PaC0, mzZero, mzZero)

  surv[1] = sum (PaC0)
  ## a is actually age + 1, since the first year of life is age 0
  for (a in 2:maxAge) {
    PaC0 = P %*% PaC0

    survUpdateSkewness[a] = esSkewnessVec %*% 
      c(mzZero, Sbullet %*% PaC0, mzZero)  
    growthUpdateSkewness[a] = esSkewnessVec %*% 
      c(mzZero, mzZero, P %*% PaC0) 
    fecUpdateSkewness[a] = esSkewnessVec %*% 
      c(P %*% PaC0, mzZero, mzZero)

    survUpdateVar[a] = esMu2Vec %*% c(mzZero, Sbullet %*% PaC0, mzZero)  
    growthUpdateVar[a] = esMu2Vec %*% c(mzZero, mzZero, P %*% PaC0) 
    fecUpdateVar[a] = esMu2Vec %*% c(P %*% PaC0, mzZero, mzZero)

    ## survival trajectory skewness and var.
    survTrajecSkewness[a] = fecUpdateSkewness[a-1] -
      survUpdateSkewness[a]
    survTrajecVar[a] = fecUpdateVar[a-1] - survUpdateVar[a]
    
    ## growth trajectory skewness and var.
    growthTrajecSkewness[a-1] = survUpdateSkewness[a-1] -
      growthUpdateSkewness[a-1] 
    growthTrajecVar[a-1] = survUpdateVar[a-1] - growthUpdateVar[a-1]

    ## fecundity skewness and var.
    fecSkewness[a-1] = growthUpdateSkewness[a-1] - fecUpdateSkewness[a-1]
    fecVar[a-1] = growthUpdateVar[a-1] - fecUpdateVar[a-1]
    
    surv[a] = sum(PaC0)
  }

  ## By what age is a proportion survThreshold of individuals dead?
  lifespan = which(surv < survThreshold)[1]

  ## And get the first value of survTrajecSkewness and survTrajecVar
  survSkewness[1] = birthStageEnvSkewness - survUpdateSkewness[1]
  survVar[1] = birthStageEnvVar - survUpdateVar[1]

  ## sanity checks on variance contributions
  if (min(survTrajecVar) < 0)
    stop("Survival trajectory luck contributions to variance are negative.")
  if (min(growthTrajecVar) < 0)
    stop("Growth trajectory luck contributions to variance are negative.")
  if (min(fecVar) < 0)
    stop("Fecundity luck contributions to variance are negative.")
  if (birthStateVar < 0)
    stop("Birth state contribution to variance is negative.")
  
  totSurvTrajecSkewness = sum(survTrajecSkewness)
  totGrowthTrajecSkewness = sum(growthTrajecSkewness)
  totFecSkewness = sum(fecSkewness)
  totSkewness = bsSkewnessVec %*% bsC0
  ## Sanity check
  foo = totSurvTrajecSkewness + totGrowthTrajecSkewness +
    totFecSkewness + birthStateSkewness
  if ((foo < (1 - percentTol)*totSkewness) |
      (foo > (1 + percentTol)*totSkewness))
    warning("Skewness components do not sum to total skewness as calculated by bsSkewnessVec.")

  totSurvTrajecVar = sum(survTrajecVar)
  totGrowthTrajecVar = sum(growthTrajecVar)
  totFecVar = sum(fecVar)
  totVar = bsMu2Vec %*% bsC0
  ## Sanity check
  foo = totSurvTrajecVar + totGrowthTrajecVar +
      totFecVar + birthStateVar
  if ((foo < (1 - percentTol)*totVar) |
      (foo > (1 + percentTol)*totVar))
    warning("Variance components do not sum to total variance as calculated by bsMu2Vec.")

  return (out=list(birthStateVar=birthStateVar,
                   survTrajecVar=survTrajecVar,
                   growthTrajecVar=growthTrajecVar,
                   fecVar=fecVar,
                   totVar=totVar,
                   birthStateSkewness=birthStateSkewness,
                   survTrajecSkewness=survTrajecSkewness,
                   growthTrajecSkewness=growthTrajecSkewness,
                   fecSkewness=fecSkewness,
                   totSkewness=totSkewness,
                   lifespan=lifespan))
}


#' Partitions variance and skewness of LRO
#' 
#' Partitions variance and skewness of lifetime reproductive output
#' (LRO) for models with environmental variation but without trait
#' variation, assuming a post-breeding census.
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
#' @param maxAge The maximum age an individual can attain.  Optional.
#'   The default value is 100.
#' @param esR1 The 1st order reward matrix for the extended state
#'   space (see Details).  Optional.  The default assumes
#'   Poisson-distributed clutch sizes.
#' @param esR2 The 2nd order reward matrix for the extended state
#'   space (see Details).  Optional.  The default assumes
#'   Poisson-distributed clutch sizes.
#' @param esR3 The 3rd order reward matrix for the extended state
#'   space (see Details).  Optional.  The default assumes
#'   Poisson-distributed clutch sizes.
#' @param bsR1 The 1st order reward matrix for the extended state
#'   space that includes the pre-birth "ur-state" (see Details).
#'   Optional.  The default assumes Poisson-distributed clutch sizes.
#' @param bsR2 The 2nd order reward matrix for the extended state
#'   space that includes the pre-birth "ur-state" (see Details).
#'   Optional.  The default assumes Poisson-distributed clutch sizes.
#' @param bsR3 The 3rd order reward matrix for the extended state
#'   space that includes the pre-birth "ur-state" (see Details).
#'   Optional.  The default assumes Poisson-distributed clutch sizes.
#' @param survThreshold The threshold to use in determining lifespan
#'   (see return values).  Optional.  Default value is 0.05.
#' @param Fdist The clutch size distribution.  The recognized options
#'   are "Poisson" and "Bernoulli".  Optional.  The default value is
#'   "Poisson".
#' @details A post-breeding census is assumed: reproduction happens
#'   after survival and growth.  The details of this calculation,
#'   including the definitions of the extended states, can be found in
#'   Robin E. Snyder and Stephen P. Ellner.  2024.
#'   "To prosper, live long: Understanding the sources of reproductive
#'   skew and extreme reproductive success in structured populations." 
#'   The American Naturalist 204(2) and its online supplemnt.
#' @return A list containing the following: * birthStateVar: the
#'   contribution to Var(LRO) from birth state luck * birthEnvVar: the
#'   contribution to Var(LRO) from birth environment luck *
#'   survTrajecVar: a vector whose jth entry contains the contribution
#'   to Var(LRO) from survival trajectory luck at age j-1.  *
#'   growthTrajecVar: a vector whose jth entry contains the
#'   contribution to Var(LRO) from growth trajectory luck at age j-1.
#'   * envTrajecVar: a vector whose jth entry contains the
#'   contribution to Var(LRO) from environment trajectory luck at age
#'   j-1.  * fecVar: a vector whose jth entry contains the
#'   contribution to Var(LRO) from fecundity luck at age j-1.  *
#'   totVar: the total variance in LRO * birthStateSkewness: the
#'   contribution to LRO skewness from birth state luck *
#'   birthEnvSkewness: the contribution to LRO skewness from birth
#'   environment luck * survTrajecSkewness: a vector whose jth entry
#'   contains the contribution to LRO skewness from survival
#'   trajectory luck at age j-1.  * growthTrajecSkewness: a vector
#'   whose jth entry contains the contribution to LRO skewness from
#'   growth trajectory luck at age j-1.  * envTrajecSkewness: a vector
#'   whose jth entry contains the contribution to LRO skewness from
#'   environment trajectory luck at age j-1.  * fecSkewness: a vector
#'   whose jth entry contains the contribution to LRO skewness from
#'   fecundity luck at age j-1.  * totSkewness: the total skewness in
#'   LRO * lifespan: the age by which a proportion survThreshold of a
#'   cohort will be dead
#' @examples
#' P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' P2 = matrix (c(0, 0.2, 0, 0, 0, 0.3, 0, 0, 0.4), 3, 3)
#' F1 = matrix (0, 3, 3); F1[1,] = 0:2
#' F2 = matrix (0, 3, 3); F1[1,] = 0.8*(0:2)
#' Plist = list (P1, P2)
#' Flist = list (F1, F2)
#' Q = matrix (1/2, 2, 2)
#' c0 = c(1,0,0)
#' out = partitionVarSkewnessEnvVarPostBreeding (Plist, Flist, Q, c0)
## BUG: Growth trajec. var. is negative, as is fec. var.
partitionVarSkewnessEnvVarPostBreeding = function (Plist, Flist, Q, c0,
                                       maxAge=100, 
                                       survThreshold=0.05,
                                       Fdist="Poisson",
                                       esR1=NULL, esR2=NULL, esR3=NULL,
                                       bsR1=NULL, bsR2=NULL, bsR3=NULL)
{
  ##  require (Matrix)
  
  ## tolerance for error checking.  0.001 = 0.1% tolerance
  percentTol = 0.001
  
  mz = dim(Plist[[1]])[1]
  numEnv = dim(Q)[1]
  bigmz = numEnv*mz

  ## Sanity check input
  if (Fdist == "Bernoulli") {
    for (q in 1:numEnv)
      if (sum(colSums(Flist[[q]]) > 1))
        stop("Probability of having an offspring > 1!  Columns of fecundity matrix in environment ",
             q, " sum to > 1 but clutch size is Bernoulli-distributed.")
  }

  ## Sanity check input
  if (length(c0) != mz)
    stop ("Length of c0 does not match dimension of P")

  ## Sanity check input
  if (dim(Flist[[1]])[1] != mz)
    stop ("P and F should have the same dimensions.")

  ## u0 is the stationary environmental distribution, given by the
  ## dominant eigenvector of Q
  u0 = eigen(Q)$vectors[,1]
  u0 = u0 / sum(u0)

  m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)

  ## Define megamatrices M and F
  F = M = matrix (0, bigmz, bigmz)
  for (i in 1:numEnv) {
    for (j in 1:numEnv) {
      M[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Plist[[j]]*Q[i,j]
      F[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Flist[[j]]*Q[i,j]
    }
  }

  ## Make survival and growth matrices
  Slist = Glist = list (numEnv)
  for (q in 1:numEnv) {
    Slist[[q]] = diag (colSums (Plist[[q]]))
    Sinv = diag (1 / diag(Slist[[q]]))
    Glist[[q]] = Plist[[q]] %*% Sinv
  }

  ## Make "bullet" matrices.  bdiag makes sparse matrices, which can't
  ## be pasted into regular matrices, so wrap these in as.matrix().
  Sbullet = as.matrix(bdiag (Slist))
  Gbullet = as.matrix(bdiag (Glist))
  Pbullet = Gbullet %*% Sbullet
  Qbullet = matrix (0, bigmz, bigmz)
  for (i in 1:numEnv) {
    for (j in 1:numEnv) {
      Qbullet[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = diag(mz)*Q[i,j]
    }
  }
  ## This code assumes that #kids is not part of the state definition,
  ## so Fbullet = I
  Fbullet = diag(bigmz)

  ## Sanity check
  epsilon = 0.00001
  if (sum(abs(range(Qbullet %*% Gbullet %*% Sbullet - M))) > epsilon)
    warning("Qbullet %*% Gbullet %*% Sbullet does not equal M.")

  ## create matrices for extended state space, which I will denote with
  ## the prefix "es".

  esbigmz = 4*bigmz
  esF = esM = matrix (0, esbigmz, esbigmz)

  esM[bigmz + 1:bigmz, 1:bigmz] = Fbullet
  esM[2*bigmz + 1:bigmz, bigmz + 1:bigmz] = Sbullet
  esM[3*bigmz + 1:bigmz, 2*bigmz + 1:bigmz] = Gbullet
  esM[1:bigmz, 3*bigmz + 1:bigmz] = Qbullet
  
  esF[bigmz + 1:bigmz, 1:bigmz] = F

  ################################################################
  ## Extended state model that includes ur-stage alpha_z and
  ## ur-environment alpha_q.  Used for calculating contributions from
  ## birth stage and birth environment.  The "bs" before variable names
  ## stands for "birth state."
  ################################################################

  bsbigmz = 1 + mz + bigmz
  ## initial state
  bsM0 = rep(0, bsbigmz)
  bsM0[1] = 1

  bsF = bsM = matrix (0, bsbigmz, bsbigmz)
  ## From ur-stage to birth stage and ur-environment
  ## BUG.  Should c0 be m0?  No, I don't think so.  Should c0 be
  ## passed in?
  bsM[1 + (1:mz),1] = c0
  ## From birth stage and ur-environment to m0
  for (q in 1:numEnv) {
    bsM[1 + q*mz + 1:mz, 1 + 1:mz] = diag (u0[q], mz)
  }
  ## From normal stage x env. state to normal stage x env. state
  bsM[1 + mz + 1:bigmz, 1 + mz + 1:bigmz] = M
  ## You only start reproducing once you're past the ur-stages.  I guess
  ## reproduction should go to the ur-stage x ur-environment state.  Not
  ## that it matters much, since all we need are the column sums of bsF.
  bsF[1, 1 + mz + 1:bigmz] = colSums (F)

  ################################################################
  ## set up reward matrices
  ################################################################
  
  if (is.null (esR1) & is.null (esR2) & is.null (esR3) &
        is.null (bsR1) & is.null (bsR2) & is.null (bsR3)) {
    esR1 = esR2 = esR3 = matrix (0, esbigmz+1, esbigmz+1)
    bsR1 = bsR2 = bsR3 = matrix (0, bsbigmz+1, bsbigmz+1)

    ## First moment of clutch size
    for (j in 1:esbigmz) 
      esR1[,j] = sum(esF[,j])
    for (j in 1:bsbigmz) 
      bsR1[,j] = sum(bsF[,j])

    if (Fdist == "Poisson") {
      ## Second moment of clutch size  
      esR2 = esR1 + esR1^2
      bsR2 = bsR1 + bsR1^2
      ## Third moment of clutch size
      esR3 = esR1+ 3*esR1*esR1 + esR1*esR1*esR1
      bsR3 = bsR1+ 3*bsR1*bsR1 + bsR1*bsR1*bsR1
    } else if (Fdist == "Bernoulli") {
      esR2 = esR3 = esR1
      bsR2 = bsR3 = bsR1
    } else {
      stop("Currently only supports Poisson- and Bernoulli-distributed clutch sizes.")
    }
  } else if (!is.null (esR1) & !is.null (esR2) & !is.null (esR3) &
               !is.null (bsR1) & !is.null (bsR2) & !is.null (bsR3)) {
    ## User-supplied reward matrices are fine too.
  } else {
    stop("Values need to be specified for all or none of esR1, esR2, esR3, bs$1, bsR2, bsR3.")
  }

  message ("Calculating es moments...\n")
  out = calcMoments (esM, esR1, esR2, esR3)
  esRho1Vec = out$rho1Vec
  esMu2Vec = out$mu2Vec
  esSkewnessVec = out$skewnessVec

  message ("Calculating bs moments...\n")
  out = calcMoments (bsM, bsR1, bsR2, bsR3)
  bsRho1Vec = out$rho1Vec
  bsMu2Vec = out$mu2Vec
  bsMu3Vec = out$mu3Vec
  bsSkewnessVec = out$skewnessVec

  ## Start calculating luck terms
  survTrajecSkewness = growthTrajecSkewness = envTrajecSkewness =
    fecSkewness = numeric (maxAge)
  fecUpdateSkewness = survUpdateSkewness =
    growthUpdateSkewness = envUpdateSkewness = numeric (maxAge)
  survUpdateVar = growthUpdateVar = envUpdateVar = fecUpdateVar =
    numeric (maxAge) 
  survTrajecVar = growthTrajecVar = envTrajecVar = fecVar =
    numeric (maxAge)
  surv = numeric (maxAge)

  ## for ur-state var.
  ## expectation wrt z of Var(success | z)
  expZVarCondZ = bsMu2Vec %*% bsM0
  ## variance wrt z of Exp(success | z)
  varZExpCondZ = sum(bsM0*bsRho1Vec^2) - (sum(bsM0*bsRho1Vec)^2)
  ## Law of total variance
  urVar = expZVarCondZ + varZExpCondZ
  
  justBirthStageVar = bsMu2Vec %*% bsM %*% bsM0
  birthStageEnvVar = bsMu2Vec %*% bsM %*% bsM %*% bsM0
  
  birthStateVar = urVar - justBirthStageVar
  birthEnvVar = justBirthStageVar - birthStageEnvVar

  ## for ur-state skewness
  expZMu3CondZ = bsMu3Vec %*% bsM0
  ## 3rd moment wrt z of Exp(success | z)
  rho3ZExpCondZ = sum(bsM0*(bsRho1Vec)^3)
  ## 2nd moment wrt z of Exp(success | z)
  rho2ZExpCondZ = sum(bsM0*(bsRho1Vec)^2)
  ## 1st moment wrt z of Exp(success | z)
  expZExpCondZ = sum(bsM0*bsRho1Vec)
  ## 3rd central moment wrt z of Exp(success | z)
  mu3ZExpCondZ = rho3ZExpCondZ -
    3*expZExpCondZ*rho2ZExpCondZ +
    2*expZExpCondZ^3
  ## covariance wrt z of Exp(success | z) and Var(success | z)
  covZExpVarCondZ = sum(bsM0*bsRho1Vec*bsMu2Vec) -
    sum(bsM0*bsRho1Vec)*sum(bsM0*bsMu2Vec)
  ## Law of total cumulance
  urMu3 = expZMu3CondZ + mu3ZExpCondZ +
    3*covZExpVarCondZ
  urSkewness = urMu3 / urVar^1.5
  
  justBirthStageSkewness = bsSkewnessVec %*% bsM %*% bsM0
  birthStageEnvSkewness = bsSkewnessVec %*% bsM %*% bsM %*% bsM0

  birthStateSkewness = urSkewness - justBirthStageSkewness
  birthEnvSkewness = justBirthStageSkewness - birthStageEnvSkewness

  ## Starting the age loop
  MaM0 = m0
  bigmzZero = rep(0, bigmz)

  ## Sanity check --- passes
  N = solve(diag(bigmz) - M)
  rho1Vec = colSums (F %*% N)
  foo1 = esRho1Vec %*% c(bigmzZero, bigmzZero, bigmzZero, m0)
  foo2 = rho1Vec %*% m0
  if ((foo1 < (1 - percentTol)*foo2) |
      (foo1 > (1 + percentTol)*foo2))
    warning("Calculating the mean via esRho1Vec doesn't get the same answer as calculating it via rho1Vec.")

  survUpdateSkewness[1] = esSkewnessVec %*%
    c(bigmzZero, Sbullet %*% MaM0, bigmzZero, bigmzZero) 
  growthUpdateSkewness[1] = esSkewnessVec %*%
    c(bigmzZero, bigmzZero, Pbullet %*% MaM0, bigmzZero)
  fecUpdateSkewness[1] = esSkewnessVec %*%
    c(bigmzZero, bigmzZero, bigmzZero, Fbullet %*% Pbullet %*% MaM0)
  envUpdateSkewness[1] = esSkewnessVec %*%
    c(M %*% MaM0, bigmzZero, bigmzZero, bigmzZero)
  
  survUpdateVar[1] = esMu2Vec %*%
    c(bigmzZero, Sbullet %*% MaM0, bigmzZero, bigmzZero) 
  growthUpdateVar[1] = esMu2Vec %*%
    c(bigmzZero, bigmzZero, Pbullet %*% MaM0, bigmzZero)
  fecUpdateVar[1] = esMu2Vec %*%
    c(bigmzZero, bigmzZero, bigmzZero, Fbullet %*% Pbullet %*% MaM0)
  envUpdateVar[1] = esMu2Vec %*%
    c(M %*% MaM0, bigmzZero, bigmzZero, bigmzZero)

  surv[1] = sum (MaM0)

  for (a in 2:maxAge) {
    MaM0 = M %*% MaM0

    survUpdateSkewness[a] = esSkewnessVec %*%
      c(bigmzZero, Sbullet %*% MaM0, bigmzZero, bigmzZero) 
    growthUpdateSkewness[a] = esSkewnessVec %*%
      c(bigmzZero, bigmzZero, Pbullet %*% MaM0, bigmzZero)
    fecUpdateSkewness[a] = esSkewnessVec %*%
      c(bigmzZero, bigmzZero, bigmzZero, Fbullet %*% Pbullet %*% MaM0)
    envUpdateSkewness[a] = esSkewnessVec %*%
      c(M %*% MaM0, bigmzZero, bigmzZero, bigmzZero)
   
    survUpdateVar[a] = esMu2Vec %*%
      c(bigmzZero, Sbullet %*% MaM0, bigmzZero, bigmzZero) 
    growthUpdateVar[a] = esMu2Vec %*%
      c(bigmzZero, bigmzZero, Pbullet %*% MaM0, bigmzZero)
    fecUpdateVar[a] = esMu2Vec %*%
      c(bigmzZero, bigmzZero, bigmzZero, Fbullet %*% Pbullet %*% MaM0)
    envUpdateVar[a] = esMu2Vec %*%
      c(M %*% MaM0, bigmzZero, bigmzZero, bigmzZero)

    ## survival trajectory skewness and var.
    survTrajecSkewness[a] = fecUpdateSkewness[a-1] -
      survUpdateSkewness[a]
    survTrajecVar[a] = fecUpdateVar[a-1] - survUpdateVar[a]

    ## growth trajectory skewness and var.
    growthTrajecSkewness[a-1] = survUpdateSkewness[a-1] -
      growthUpdateSkewness[a-1] 
    growthTrajecVar[a-1] = survUpdateVar[a-1] - growthUpdateVar[a-1]

    ## fecundity skewness and var.
    fecSkewness[a-1] = envUpdateSkewness[a-1] - fecUpdateSkewness[a-1]
    fecVar[a-1] = envUpdateVar[a-1] - fecUpdateVar[a-1]

    ## environment trajectory skewness
    envTrajecSkewness[a-1] = growthUpdateSkewness[a-1] -
      envUpdateSkewness[a-1]
    envTrajecVar[a-1] = growthUpdateVar[a-1] - envUpdateVar[a-1]

    surv[a] = sum (MaM0)
  }

  ## And get the first value of survTrajecSkewness and survTrajecVar
  survTrajecSkewness[1] = birthStageEnvSkewness - survUpdateSkewness[1]
  survTrajecVar[1] = birthStageEnvVar - survUpdateVar[1]

  ## sanity checks on variance contributions
  if (min(survTrajecVar) < 0)
    stop("Survival trajectory luck contributions to variance are negative.")
  if (min(growthTrajecVar) < 0)
    stop("Growth trajectory luck contributions to variance are negative.")
  if (min(fecVar) < 0)
    stop("Fecundity luck contributions to variance are negative.")
  if (min(envTrajecVar) < 0)
    stop("Environment trajectory luck contributions to variance are negative.")
  if (birthStateVar < 0)
    stop("Birth state contribution to variance is negative.")
  if (birthEnvVar < 0)
    stop("Birth environment contribution to variance is negative.")
  
  ## By what age are survThreshold propor. of individuals dead?
  lifespan = which(surv < survThreshold)[1]

  totSurvTrajecSkewness = sum(survTrajecSkewness)
  totGrowthTrajecSkewness = sum(growthTrajecSkewness)
  totEnvTrajecSkewness = sum(envTrajecSkewness)
  totFecSkewness = sum(fecSkewness)
  totSkewness = bsSkewnessVec %*% bsM0

  ## sanity check
  foo = totSurvTrajecSkewness + totGrowthTrajecSkewness +
         totEnvTrajecSkewness + totFecSkewness +
         birthStateSkewness + birthEnvSkewness
  if ((foo < (1 - percentTol)*totSkewness) |
      (foo > (1 + percentTol)*totSkewness))
    warning("Skewness luck components do not sum to total skewness as calculated by bsSkewnessVec.")
  
  totSurvTrajecVar = sum(survTrajecVar)
  totGrowthTrajecVar = sum(growthTrajecVar)
  totEnvTrajecVar = sum(envTrajecVar)
  totFecVar = sum(fecVar)
  totVar = bsMu2Vec %*% bsM0
  ## Sanity check
  foo = totSurvTrajecVar +
         totGrowthTrajecVar +
         totEnvTrajecVar + totFecVar +
         birthStateVar + birthEnvVar
  if ((foo < (1 - percentTol)*totVar) |
      (foo > (1 + percentTol)*totVar))
    warning("Variance luck components do not sum to total variance as calculated by bsMu2Vec.")

  return (out=list(birthStateVar=birthStateVar,
                   birthEnvVar=birthEnvVar,
                   survTrajecVar=survTrajecVar,
                   growthTrajecVar=growthTrajecVar,
                   envTrajecVar=envTrajecVar,
                   fecVar=fecVar,
                   totVar=totVar,
                   birthStateSkewness=birthStateSkewness,
                   birthEnvSkewness=birthEnvSkewness,
                   survTrajecSkewness=survTrajecSkewness,
                   growthTrajecSkewness=growthTrajecSkewness,
                   envTrajecSkewness=envTrajecSkewness,
                   fecSkewness=fecSkewness,
                   totSkewness=totSkewness,
                   lifespan=lifespan))
}
