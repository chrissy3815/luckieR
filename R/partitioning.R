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
#' @param debugging Will print the results of various sanity checks if
#'   debugging is set to TRUE.  Optional.  Default value is FALSE.
#' @details The details of this calculation, including the definitions
#'   of the extended states can be found in Robin E. Snyder and
#'   Stephen P. Ellner.  2024.  "To prosper, live long: Understanding
#'   the sources of reproductive skew and extreme reproductive success
#'   in structured populations."  The American Naturalist 204(2) and
#'   its online supplemnt.
#' @return A list containing the following:
#' * birthStateVar: the contribution to Var(LRO) from birth state luck
#' * survTrajecVar: a vector whose jth entry contains the contribution
#'   to Var(LRO) from survival trajectory luck at age j-1.
#' * growthTrajecVar: a vector whose jth entry contains the contribution
#'   to Var(LRO) from growth trajectory luck at age j-1.
#' * fecVar: a vector whose jth entry contains the contribution
#'   to Var(LRO) from fecundity luck at age j-1.
#' * totVar: the total variance in LRO
#' * survTrajecSkewness: a vector whose jth entry contains the contribution
#'   to LRO skewness from survival trajectory luck at age j-1.
#' * growthTrajecSkewness: a vector whose jth entry contains the contribution
#'   to LRO skewness from growth trajectory luck at age j-1.
#' * fecSkewness: a vector whose jth entry contains the contribution
#'   to LRO skewness from fecundity luck at age j-1.
#' * totSkewness: the total skewness in LRO
#' * lifespan: the age by which 99% of a cohort is expected to be dead
#' @examples
#' P = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
#' F = matrix (0, 3, 3); F[1,] = 0:2
#' c0 = c(1,0,0)
#' out = getVarSkewnessPartitionsNoEnvVar (P, F, c0)
getVarSkewnessPartitionsNoEnvVar = function (P, F, c0, maxAge=100,
                                             esR1=NULL, esR2=NULL,
                                             esR3=NULL,
                                             bsR1=NULL, bsR2=NULL,
                                             bsR3=NULL, debugging=FALSE) {
  # Input parameters:
  # P is the survival/growth transition matrix (called U in COMPADRE)
  # F is the matrix of fertility transitions
  # c0 is the birth state distribution (also called mixing distribution)

  mz = dim(P)[1]

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
  cat(range(Gbullet%*%Sbullet-P), "should equal 0", "\n"); 
  
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

  ## If the reward matrices aren't specified, we assume
  ## Poisson-distributed clutch sizes
  if (is.null(esR1)) {
    esR1 = esR2 = esR3 = matrix (0, esmz+1, esmz+1)
    bsR1 = bsR2 = bsR3 = matrix (0, bsmz+1, bsmz+1)

    ## First moment of clutch size (Poisson)
    for (j in 1:esmz) 
      esR1[,j] = sum(esF[,j])
    for (j in 1:bsmz) 
      bsR1[,j] = sum(bsF[,j])

    ## Second moment of clutch size  
    esR2 = esR1 + esR1^2
    bsR2 = bsR1 + bsR1^2
    ## Third moment of clutch size
    esR3 = esR1+ 3*esR1*esR1 + esR1*esR1*esR1
    bsR3 = bsR1+ 3*bsR1*bsR1 + bsR1*bsR1*bsR1
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

  ## Sanity check: passes
  if (debugging) {
    N = solve(diag(mz) - P)
    rho1Vec = colSums (F %*% N)
    cat (esRho1Vec %*% c(c0, mzZero, mzZero),
         "should = ", rho1Vec %*% c0, "\n")
  }

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

  ## By what age are 99% of individuals dead?
  lifespan = which(surv < 0.01)[1]

  ## And get the first value of fecSkewness and fecVar
  fecSkewness[1] = justBirthStateSkewness - fecUpdateSkewness[1]
  fecVar[1] = justBirthStateVar - fecUpdateVar[1]

  ## By what age are 99% of individuals dead?
  lifespan = which(surv < 0.01)[1]

  totSurvTrajecSkewness = sum(survTrajecSkewness)
  totGrowthTrajecSkewness = sum(growthTrajecSkewness)
  totFecSkewness = sum(fecSkewness)
  totSkewness = bsSkewnessVec %*% bsC0
  ## Sanity check: passes
  if (debugging) 
    cat (totSurvTrajecSkewness + totGrowthTrajecSkewness +
         totFecSkewness + birthStateSkewness, "should = ",
         totSkewness, "\n")

  totSurvTrajecVar = sum(survTrajecVar)
  totGrowthTrajecVar = sum(growthTrajecVar)
  totFecVar = sum(fecVar)
  totVar = bsMu2Vec %*% bsC0
  ## Sanity check: passes
  if (debugging)
    cat (totSurvTrajecVar + totGrowthTrajecVar +
         totFecVar + birthStateVar, "should = ",
         totVar, "\n")

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


getVarSkewnessPartitionsEnvVar = function (Plist, Flist, Q,
                                           c0, u0, 
                                           maxAge=100, esR1=NULL,
                                           esR2=NULL, esR3=NULL,
                                           bsR1=NULL, bsR2=NULL,
                                           bsR3=NULL)
{
  require (Matrix)
  debugging = TRUE
  
  mz = dim(Plist[[1]])[1]
  numEnv = dim(Q)[1]
  bigmz = numEnv*mz

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
  if (debugging)
    cat (sum(Qbullet %*% Gbullet %*% Sbullet - M), "should = 0.\n")

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
  
  ## If the reward matrices aren't specified, we assume
  ## Poisson-distributed clutch sizes
  if (is.null(esR1)) {
    esR1 = esR2 = esR3 = matrix (0, esbigmz+1, esbigmz+1)
    bsR1 = bsR2 = bsR3 = matrix (0, bsbigmz+1, bsbigmz+1)

    ## First moment of clutch size (Poisson)
    for (j in 1:esbigmz) 
      esR1[,j] = sum(esF[,j])
    for (j in 1:bsbigmz) 
      bsR1[,j] = sum(bsF[,j])

    ## Second moment of clutch size  
    esR2 = esR1 + esR1^2
    bsR2 = bsR1 + bsR1^2
    ## Third moment of clutch size
    esR3 = esR1+ 3*esR1*esR1 + esR1*esR1*esR1
    bsR3 = bsR1+ 3*bsR1*bsR1 + bsR1*bsR1*bsR1
  }

  cat ("Calculating es moments...\n")
  out = calcMoments (esM, esR1, esR2, esR3)
  esRho1Vec = out$rho1Vec
  esMu2Vec = out$mu2Vec
  esSkewnessVec = out$skewnessVec

  cat ("Calculating bs moments...\n")
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
  if (debugging) {
    N = solve(diag(bigmz) - M)
    rho1Vec = colSums (F %*% N)
    cat ("Checking Exp(R):", esRho1Vec %*%
                             c(bigmzZero, bigmzZero, bigmzZero, m0),
         "should = ", rho1Vec %*% m0, "\n")
  }

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

  ## By what age are (1 - survThreshold) propor. of individuals dead?
##  lifespan = which(surv < survThreshold)[1]

  totSurvTrajecSkewness = sum(survTrajecSkewness)
  totGrowthTrajecSkewness = sum(growthTrajecSkewness)
  totEnvTrajecSkewness = sum(envTrajecSkewness)
  totFecSkewness = sum(fecSkewness)
  totSkewness = bsSkewnessVec %*% bsM0

  ## sanity check
  if (debugging) 
    cat ("Checking sum of skewness luck:",
         totSurvTrajecSkewness + totGrowthTrajecSkewness +
         totEnvTrajecSkewness + totFecSkewness +
         birthStateSkewness + birthEnvSkewness, "should = ",
         totSkewness, "\n")
  
  totSurvTrajecVar = sum(survTrajecVar)
  totGrowthTrajecVar = sum(growthTrajecVar)
  totEnvTrajecVar = sum(envTrajecVar)
  totFecVar = sum(fecVar)
  totVar = bsMu2Vec %*% bsM0
  ## Sanity check: passes
  if (debugging)
    cat ("Checking sum of variance luck:",
         totSurvTrajecVar +
         totGrowthTrajecVar +
         totEnvTrajecVar + totFecVar +
         birthStateVar + birthEnvVar, "should = ",
         totVar, "\n")


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
                   surv=surv))
}


## Partitions variance and skewness of LRO in the presence of both
## env. variation and trait variation.  Includes the contribution from
## trait variation.
partitionVarSkewnessWithEnvAndTraits = function (numTraits, numEnv,
                                                 PlistAllTraits,
                                                 FlistAllTraits,
                                                 Q, c0, u0, maxAge,
                                                 survThreshold=0.05) {

  mz = length (c0)
  numEnv = length(u0)
  bigmz = mz*numEnv
  
  ## Initial cross-classified state
  m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)

  lifespanCondX = birthEnvSkewnessCondX =
    birthEnvVarCondX = birthStateSkewnessCondX =
      birthStateVarCondX = numeric (numTraits)
  expRCondX = totVarCondX = totSkewnessCondX = numeric (numTraits)

  survTrajecVarCondX = growthTrajecVarCondX = envTrajecVarCondX = fecVarCondX =
    matrix (0, numTraits, maxAge)
  survTrajecSkewnessCondX = growthTrajecSkewnessCondX = envTrajecSkewnessCondX =
    fecSkewnessCondX = matrix (0, numTraits, maxAge)

  if (debugging) 
    expRCondXZ = r2CondX = skewnessVecCondX = matrix (0, numTraits, bigmz)

  if (debugging)
    totVarCondX2 = numeric (numTraits)

  for (x in 1:numTraits) {

    cat ("About to define Plist for x =", x, "\n")
    Plist = PlistAllTraits[[x]]
    Flist = FlistAllTraits[[x]]

    ## No need to specify reward matrices below because offspring
    ## distrib. is assumed to be Poisson.

    ## bsMu2Vec is all zeros except a stonking big negative value in
    ## the first element, which is supposed to be variance. 
    out = getVarSkewnessPartitionsEnvVar (Plist, Flist, Q, 
                                          c0, u0, maxAge)
    birthEnvSkewnessCondX[x] = out$birthEnvSkewness
    survTrajecSkewnessCondX[x,] = out$survTrajecSkewness
    growthTrajecSkewnessCondX[x,] = out$growthTrajecSkewness
    envTrajecSkewnessCondX[x,] = out$envTrajecSkewness
    fecSkewnessCondX[x,] = out$fecSkewness
    totSkewnessCondX[x] = out$totSkewness
    surv = out$surv
    lifespanCondX[x] = which(surv < survThreshold)[1]
    
    birthEnvVarCondX[x] = out$birthEnvVar
    survTrajecVarCondX[x,] = out$survTrajecVar
    growthTrajecVarCondX[x,] = out$growthTrajecVar
    envTrajecVarCondX[x,] = out$envTrajecVar
    fecVarCondX[x,] = out$fecVar
    totVarCondX[x] = out$totVar
    ## Sanity check: passes
    if (debugging) {
      totVarCondX2[x] = birthEnvVarCondX[x] +
        sum(survTrajecVarCondX[x,] + growthTrajecVarCondX[x,] +
            envTrajecVarCondX[x,] + fecVarCondX[x,])
      cat ("Do the contributions to variance add up to tot. var.?",
           totVarCondX[x], "should = ",
           totVarCondX2[x], "\n")
    }
    
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

    if (debugging) {
      ## All but the seedling class reproduce
      foo = rep (1, mz)
      foo[1] = 0
      pb = rep (foo, numEnv)
      b = colSums (F)
      sigbsq = b
      rbarPib = expRCondXZ[x,] %*% M     # \bar{r} \pi_b
      r2CondX[x,] = (sigbsq + (pb*b)^2 + 2*pb*b*rbarPib) %*% N
      varRCondXZ = r2CondX[x,] - expRCondXZ[x,]^2
      expZVarRCondXZ = varRCondXZ %*% m0
      varZExpRCondXZ = sum(m0*expRCondXZ[x,]^2) -
        sum(m0*expRCondXZ[x,])^2
            ## Passes
      cat ("Var. check:", totVarCondX[x], "should = ",
           expZVarRCondXZ + varZExpRCondXZ, "\n")
      
      R1 = R2 = R3 = matrix (0, bigmz+1, bigmz+1)
      for (j in 1:bigmz)
        R1[,j] = sum(F[,j])
      R2 = R1 + R1^2 ## Poisson
      R3 = R1 + 3*R1^2 + R1^3  ## Poisson
      out = calcMoments (M, R1, R2, R3)

      skewnessVecCondX[x,] = out$skewness
      ## rho1 matches up perfectly
      ## plot (expRCondXZ[x,], main="Ex[R | z]")
      ## lines (out$rho1Vec[1,])
      
      ## plot (r2CondX[x,], main="Ex[R^2 | z]")
      ## lines (out$mu2Vec[1,] + out$rho1Vec[1,]^2)
      
      ## Passes
      cat ("Ex (R^2) check:", r2CondX[x,] %*% m0, "should = ",
      (out$mu2Vec[1,] + out$rho1Vec[1,]^2) %*% m0, "\n")
    }

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

  cat ("Calculating the effect of trait variation...\n")

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
  ## Is this right?  fulmar code has a mix of F and M.  But kittiwake
  ## code has all F.  I think the fulmar code has a bug. :-(
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

  ## You only start reproducing once you're past the ur-stages.  I guess
  ## reproduction should go to the ur-stage x ur-environment state.  Not
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

  ## bstsR2 = bstsR1^2
  bstsR2 = bstsR1 + bstsR1^2
  ## bstsR3 = bstsR1^3
  bstsR3 = bstsR1 + 3*bstsR1^2 + bstsR1^3

  ## R1Stripped has the absorbed state cleaved off
  bstsR1Stripped = bstsR1[1:bstsbigmz,1:bstsbigmz]
  bstsR2Stripped = bstsR2[1:bstsbigmz,1:bstsbigmz]

  ## Ex(R)
  cat ("Calculating bstsRho1Vec...\n")
  bstse = rep (1, bstsbigmz+1)
  bstsN = solve (diag(bstsbigmz) - bstsM)

  bstsRho1Vec = t(bstsN) %*% bstsZ %*% t(bstsMplus * bstsR1) %*% bstse
  ## Sanity check: passes
  if (debugging)  
    cat ("Checking bstsRho1Vec:",
         range (bstsRho1Vec[(numTraits*mz + numTraits+2):
                            (numTraits*mz + numTraits+ 1 + mz*numEnv)]
                -
                expRCondXZ[1,]), "should = 0 0.\n")

  ## Ex(R^2 | x) Do I get the correct bstsRho2Vec for traits 1 and 3 (and
  ## presumably 2), all values of (z,q)?  (E.g. compare
  ## bstsRho2Vec[248--1833] with r2 for x = 1.)

  cat ("Calculating bstsRho2Vec...\n")
  bstsRho2Vec = t(bstsN) %*% (bstsZ %*% t(bstsMplus * bstsR2) %*% bstse +
                              2*t(bstsM * bstsR1Stripped) %*%
                              bstsRho1Vec)
  ## Sanity check: passes
  if (debugging) 
    cat ("Checking bstsRho2Vec:",
         range (bstsRho2Vec[(numTraits*mz + numTraits+2):
                            (numTraits*mz + numTraits+ 1 + mz*numEnv)]
                -
                r2CondX[1,]),
         "should = 0 0.\n")
      
  ## Ex(R^3 | x)
  cat ("Calculating bstsRho3Vec...\n")
  bstsRho3Vec = t(bstsN) %*% (bstsZ %*% t(bstsMplus * bstsR3) %*% bstse +
                              3*t(bstsM*bstsR2Stripped) %*% bstsRho1Vec +
                              3*t(bstsM*bstsR1Stripped) %*% bstsRho2Vec)

  bstsRho1Vec = t(bstsRho1Vec)
  bstsRho2Vec = t(bstsRho2Vec)
  bstsRho3Vec = t(bstsRho3Vec)

  ## Do I get the correct bstsMu2Vec for trait 3, all values of (z,q)?
  ## (Compare bstsMu2Vec[41:52] with varRCondXZQ for x = 3 in
  ## fulmarSurvGrowthPluckPartition.R.)
  bstsMu2Vec = bstsRho2Vec - bstsRho1Vec^2
  bstsMu3Vec = bstsRho3Vec - 3*bstsRho1Vec*bstsRho2Vec + 2*bstsRho1Vec^3
  bstsSkewnessVec = bstsMu3Vec / bstsMu2Vec^1.5

  ## sanity check
  if (debugging)
    cat (range(
        bstsSkewnessVec[(numTraits*mz + numTraits+2):
                        (numTraits*mz + numTraits+ 1 + mz*numEnv)] -
        skewnessVecCondX[1,]), "should = 0 0\n")

  urState = rep(0, bstsbigmz)
  urState[1] = 1
  totSkewness = bstsSkewnessVec %*% urState
  totVar2 = bstsMu2Vec %*% urState

  ## Sanity check: passes
  if (debugging)
    cat ("Checking totVar:", totVar, "should = ", totVar2, "\n")

  skewnessFromTraits = totSkewness - traitAve (totSkewnessCondX, traitDist)
  varFromTraits = totVar - traitAve (totVarCondX, traitDist)

  ## Sanity check: passes
  if (debugging)
    cat ("Checking variance from traits:", varFromTraits, "should = ",
         varXExpRCondX, "\n")

  return (out=list(varFromTraits=varFromTraits,
                   skewnessFromTraits=skewnessFromTraits,
                   totVar=totVar, totSkewness=totSkewness,
                   birthEnvVar=birthEnvVar,
                   birthEnvSkewness=birthEnvSkewness,
                   birthStateVar=birthStateVar,
                   birthStateSkewness=birthStateSkewness,
                   survTrajecVar=survTrajecVar,
                   survTrajecSkewness=survTrajecSkewness,
                   totSurvTrajecVar=totSurvTrajecVar,
                   totSurvTrajecSkewness=totSurvTrajecSkewness,
                   growthTrajecVar=growthTrajecVar,
                   growthTrajecSkewness=growthTrajecSkewness,
                   totGrowthTrajecVar=totGrowthTrajecVar,
                   totGrowthTrajecSkewness=totGrowthTrajecSkewness,
                   envTrajecVar=envTrajecVar,
                   envTrajecSkewness=envTrajecSkewness,
                   totEnvTrajecVar=totEnvTrajecVar,
                   totEnvTrajecSkewness=totEnvTrajecSkewness,
                   fecVar=fecVar,
                   fecSkewness=fecSkewness,
                   totFecVar=totFecVar,
                   totFecSkewness=totFecSkewness,
                   birthEnvVarCondX=birthEnvVarCondX,
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
