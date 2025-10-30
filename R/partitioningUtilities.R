#############################################################################
##
## partitioningUtilities: Utility functions for partitioning variance
## and skewness, with the full kernel and with a conditional kernel.
##
#############################################################################

source ("megamatrixFunctions.R")

traitAve = function (x, traitDist) {
  return (sum(traitDist*x))
}

traitVar = function (x, traitDist) {
  return (sum(traitDist*x^2) - traitAve(x, traitDist)^2)
}

traitMu3 = function (x, traitDist) {
  return (sum(traitDist*(x - traitAve(x, traitDist))^3))
}

## Calculate probability of success conditional on initial state and
## kernel conditional becoming successful before death
makeCondKernel = function (M, transientStates) {

  ## How many states are there all told?
  numStates = dim(M)[1]
  ## Create U, the transition matrix for the transient states
  U = M[transientStates, transientStates]
  ## How many transient states are there?
  numTransient = length(transientStates)
  ## calculate a2, the probability of hitting an absorbing state in
  ## the next time step
  a2 = apply (M, 2, sum)[transientStates] - apply (U, 2, sum)
  ## calculate q2, the probability of being in Z2 before you die given
  ## that you are now in state i.
  NU = solve (diag(numTransient) - U)
  q2 = a2 %*% NU
  q2Extended = rep (1, numStates); q2Extended[transientStates] = q2

  ## Transition matrix conditional on reproducing before death.
  MCond = matrix (0, numStates, numStates)
  for (j in 1:numStates)
    MCond[,j] = M[,j] * q2Extended / q2Extended[j]
    
  return (out=list (q2Extended=q2Extended, MCond=MCond))
}

## Calculate probability of not being successful conditional on
## initial state and kernel conditional on not becoming successful
## before death
makeCondFailureKernel = function (M, transientStates) {

  ## How many states are there all told?
  numStates = dim(M)[1]
  ## Create U, the transition matrix for the transient states
  U = M[transientStates, transientStates]
  ## How many transient states are there?
  numTransient = length(transientStates)
  if (FALSE) {
    ## calculate a1, the probability of dying in the next time step
    a1 = (1 - apply (M, 2, sum))[transientStates]
    ## calculate q1, the probability of dying before hitting the
    ## "successful" absorbing state given that you are now in state i.
    NU = solve (diag(numTransient) - U)
    q1 = a1 %*% NU
  }
  ## calculate a2, the probability of hitting an absorbing state in
  ## the next time step
  a2 = apply (M, 2, sum)[transientStates] - apply (U, 2, sum)
  ## calculate q2, the probability of being in Z2 before you die given
  ## that you are now in state i.
  NU = solve (diag(numTransient) - U)
  q2 = a2 %*% NU
  q1 = 1 - q2
  
  ## Transition matrix conditional on not reproducing before death.
  MCond = matrix (0, numTransient, numTransient)
  for (j in 1:numTransient)
    MCond[,j] = M[transientStates,transientStates[j]] * q1 / q1[j]

  return (out=list (q1=q1, MCond=MCond))
}

## Calculate 1st (mean) and 2nd and 3rd (non-central) moments (mu2 and mu3) and
## skewness of given transition matrix and associated fecundity
## matrix.  Inputs: reward matrices for first, second, third moments.
calcMoments = function (M, R1, R2, R3) {

  numStates = dim(M)[1]

  ## set up Mplus, the transition matrix with absorbing state (dead)
  Mplus = matrix (0, numStates+1, numStates+1)
  Mplus[1:numStates, 1:numStates] = M
  Mplus[numStates+1, 1:numStates] = 1 - colSums(M)
  Mplus[numStates+1, numStates+1] = 1

  ## Set up Z, which cleaves off dead states
  Z = cbind (diag(numStates), rep(0, numStates))

  ## Set up the vector 1_s, where s = bigmz + 1
  oneS = rep(0, numStates+1)
  oneS[numStates+1] = 1

  ## R1Stripped has the absorbed state cleaved off
  R1Stripped = R1[1:numStates, 1:numStates]
  R2Stripped = R2[1:numStates, 1:numStates]

  ## Get zero-centered moments rho1, rho2, rho3

  ## Ex(R)
  e = rep (1, numStates+1)
  N = solve (diag(numStates) - M)
  cat ("Calculating rho1Vec...\n")
  rho1Vec = t(N) %*% Z %*% t(Mplus * R1) %*% e

  ## Ex(R^2 | x)
  cat ("Calculating rho2Vec...\n")
  rho2Vec = t(N) %*% (Z %*% t(Mplus * R2) %*% e +
                        2*t(M * R1Stripped) %*% rho1Vec)

  ## Ex(R^3 | x)
  cat ("Calculating rho3Vec...\n")
  rho3Vec = t(N) %*% (Z %*% t(Mplus * R3) %*% e +
                    3*t(M*R2Stripped) %*% rho1Vec +
                    3*t(M*R1Stripped) %*% rho2Vec)

  ## non-central moments
  rho1Vec = t(rho1Vec)
  rho2Vec = t(rho2Vec)
  rho3Vec = t(rho3Vec)

  ## central moments
  mu2Vec = rho2Vec - rho1Vec^2
  mu3Vec = rho3Vec - 3*rho1Vec*rho2Vec + 2*rho1Vec^3

  ## skewness. NOTE, skewness=0 whenever the variance=0. 
  skewnessVec = rep(0,length(mu3Vec)); 
  e = which(mu2Vec>0); 
  skewnessVec[e] = mu3Vec[e]/mu2Vec[e]^1.5

  return (out=list(rho1Vec=rho1Vec, mu2Vec=mu2Vec, mu3Vec=mu3Vec,
                   skewnessVec=skewnessVec))
}



######################################################### 
# Useful for making a kernel conditional on breeding or not breeding.
# Function to make big P, Definition 3
#########################################################
makePDef3 = function (P, F) {
  ## The number of original states
  mz = dim(P)[1]
  ## The number of transient states (the nonbreeding states) 
##  mzS = length(nonbreedingStates)
  
  esmz = 3*mz
  bigF = bigP = P0 = matrix (0, esmz, esmz)

  ## calculate pd, the probability of producing at least one offspring
  ## in the current year.
  b = colSums(F)
  pd = 1 - exp(-b)
  
  ## Make P for the extended state space Z1 \cup Z2 \cup Z3
  # Z1 to Z1 
  for (j in 1:mz)  bigP[j, 1:mz] = P[j,]*(1 - pd[j])
  # Z1 to Z2
  for (j in 1:mz)  bigP[mz+j, 1:mz] = P[j,]*pd[j]
  # Z2 to Z3
  bigP[(2*mz+1):esmz, (mz+1):(2*mz)] = P
  # Z3 to Z3
  bigP[(2*mz+1):esmz, (2*mz+1):esmz] = P
  
  return(bigP)
}  

########################################################################## 
# Function to make big P conditional on breeding, Definition 3
##########################################################################  
makePCondBreedDef3 = function (P, F, c0) {
  ## The number of original states
  mz = dim(P)[1]

  esmz = 3*mz
  bigF = matrix (0, esmz, esmz)
  
  ##bigP=NULL; 
  bigP = makePDef3(P, F)

  ## expanded state space init. distribution
  bigc0 = c(c0, rep(0, 2*mz))

  ## make F for the extended state space Z1 \cup Z2 \cup Z3.  Note that the
  ## only non-zero entries are from Z2 and Z3 to Z1.
  b = colSums (F)
  bigF[1, (mz+1):esmz] = b
  
  ## Calculate aM, the probability of going from Z1 to Z2 in one step
  aM = apply (bigP[(mz+1):(2*mz), 1:mz], 2, sum)
  
  ## Q2 = P restricted to Z1, N2 = (I - Q2)^{-1}
  Q2 = bigP[1:mz,1:mz];

  N2 = solve(diag(mz) - Q2)

  ## Calculate probEverBreedCondZ (called B(z-hat) in the ms), the
  ## probability of reaching Z2 before death.
  probEverBreedCondZ = c(aM %*% N2, rep(1, 2*mz))
  probNeverBreedCondZ = 1 - probEverBreedCondZ[1:mz]

  ## Calculate the transition kernel conditional on eventually breeding
  PCondBreed = matrix (NA, esmz, esmz)
  for (j in 1:esmz)
    PCondBreed[,j] = bigP[,j]*probEverBreedCondZ/probEverBreedCondZ[j]
  ## P conditional on never breeding is only defined for the "has not
  ## yet bred" states in 1:mz
  PCondNeverBreed = matrix (NA, mz, mz)
  for (j in 1:mz)
    PCondNeverBreed[,j] = bigP[1:mz,j]*probNeverBreedCondZ/probNeverBreedCondZ[j]


  
  ###################################################################
  ## Find the extended initial state bigc0 conditional on breeding
  ###################################################################
  
  ## unconditional probability of breeding
  probBreed = probEverBreedCondZ %*% bigc0
  ## unconditional probability of never breeding
  probNeverBreed = 1 - probBreed
  ## Joint probability of starting in state z and eventually breeding
  jointProbBreedAndZ = bigc0 * probEverBreedCondZ
  ## Joint probability of starting in state z and never breeding:
  ## note, only makes sense for the "has not yet bred" states in 1:mz.
  jointProbNeverBreedAndZ = c0 * probNeverBreedCondZ
  ## P(start in state z | eventually breed)
  bigc0CondBreed = jointProbBreedAndZ / rep(probBreed, esmz)
  ## P(start in state z | never breed)
  bigc0CondNeverBreed = jointProbNeverBreedAndZ / rep(probNeverBreed, mz)

  out = list (PCondBreed=PCondBreed, probEverBreedCondZ=probEverBreedCondZ,
              PCondNeverBreed=PCondNeverBreed,
              probNeverBreedCondZ=probNeverBreedCondZ,
              bigc0CondBreed = bigc0CondBreed,
              bigc0CondNeverBreed = bigc0CondNeverBreed,
              pM=sum(bigc0*probEverBreedCondZ), bigF=bigF)
  return (out)
}

## Useful for getting modal time to first reproduction when states are
## defined by reproductive status.
getModalTimeToHitState = function (absorbingStates,
                                   Plist, Flist, Q=NULL,
                                   m0, maxAge=100) {
  mz = dim(Plist[[1]])[1]
  if (is.null(Q)) {  ## no env. var.
    bigmz = mz
    M = P
  } else {
    numEnv = dim(Q)[1]
    bigmz = numEnv*mz
    M = matrix (0, bigmz, bigmz)
    for (i in 1:numEnv) {
      for (j in 1:numEnv) {
        M[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Plist[[j]]*Q[i,j]
        F[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Flist[[j]]*Q[i,j]
      }
    }
  }

  ## Create kernel conditional on hitting an absorbing state
  transient = (1:bigmz)[-absorbingStates]
  numTransient = length(transient)

  out = makeCondKernel (M, transient)
  MCond = out$MCond 

  ## Create kernel with all absorbing states going directly to a
  ## special absorbing state.
  MAbsorbing = matrix (0, bigmz+1, bigmz+1)
  MAbsorbing[transient, transient] = MCond[transient, transient]
  for (j in transient)
    MAbsorbing[bigmz+1, j] = sum(MCond[absorbingStates, j])
  MAbsorbing[bigmz+1, absorbingStates] = 1
  MAbsorbing[bigmz+1, bigmz+1] = 1

  ## Get distribution of times to absorption
  m0Absorbing = c(m0, 0)
  MaAbsorbingM0 = matrix (0, maxAge, bigmz+1)
  MaAbsorbingM0[1,] = m0Absorbing

  for (a in 2:maxAge)
    MaAbsorbingM0[a,] = MAbsorbing %*% MaAbsorbingM0[a-1,]

  hitTimes = diff (MaAbsorbingM0[,bigmz+1])
  hitTimeMode = which.max (hitTimes)
  hitTimeDist = hitTimes / sum(hitTimes)

  return (out = list(hitTimeMode=hitTimeMode,
                     hitTimeDist=hitTimeDist))
}

getModalTimeToFirstRepro = function (Plist, Flist, Q, m0, 
                                     maxAge=100) {
  mz = dim(Plist[[1]])[1]
  numEnv = dim(Q)[1]
  bigmz = numEnv*mz

  ## Define megamatrices M and F
  F = M = matrix (0, bigmz, bigmz)
  for (i in 1:numEnv) {
    for (j in 1:numEnv) {
      M[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Plist[[j]]*Q[i,j]
      F[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Flist[[j]]*Q[i,j]
    }
  }

  ## pb = prob. of breeding as a function of state
  ## Assume entries of F are means of Poisson distributions, so
  ## prob. of breeding is pb = 1 - exp(-mean)
  pb = 1 - exp(-colSums (F))

  big3mz = 3*bigmz
  bigF = matrix (0, big3mz, big3mz)
  bigM = matrix (0, big3mz, big3mz);
  
  ## Make M for the extended state space Z1 \cup Z2 \cup Z3.  Z1 = has
  ## never reproduced, Z2 = reproduced for the first time this year,
  ## Z3 = reproduced for the first time sometime in the past.  (See
  ## Snyder and Ellner Am. Nat. 2021, SI section "Conditioning on
  ## reaching maturity, attempting breeding, or producing an
  ## offspring", "Breeder definition 3".
  Z1 = 1:bigmz
  Z2 = (bigmz+1):(2*bigmz)
  Z3 = (2*bigmz+1):(3*bigmz)
  Z23 = c(Z2, Z3)
  
  # Z1 to Z1
  for (j in Z1) bigM[j, Z1] = M[j,]*(1 - pb[j])
  # Z1 to Z2 
  for (j in Z1) bigM[bigmz+j, Z1] = M[j,]*pb[j]
  # Z2 to Z3
  bigM[Z3, Z2] = M
  # Z3 to Z3
  bigM[Z3, Z3] = M

  ## expanded state space init. distribution
  bigm0 = c(m0, rep(0, 2*bigmz))

  ## make F for the extended state space Z1 \cup Z2 \cup Z3.  Note that the
  ## only non-zero entries are from Z2 and Z3 to Z1.
  bigF[1, Z23] = F[1,]

  ## Calculate aM, the probability of going from Z1 to Z2 in one step
  aM = apply (bigM[(bigmz+1):(2*bigmz), 1:bigmz], 2, sum)
  
  ## Q2 = M restricted to Z1, N2 = (I - Q2)^{-1}
  Q2 = bigM[1:bigmz,1:bigmz];

  N2 = solve(diag(bigmz) - Q2)

  ## Calculate q2Extended, the probability of reaching Z2 before
  ## death.
  q2Extended = c(aM %*% N2, rep(1, 2*bigmz))

  ## Transition matrix conditional on becoming big before death
  bigMCond = matrix (0, big3mz, big3mz)
  for (j in 1:big3mz)
    bigMCond[,j] = bigM[,j] * q2Extended / q2Extended[j]

  ## For calculating repro. states hitting time dist.  All repro. states
  ## transition to a special absorbing state.
  bigMReproAbsorbing = matrix (0, big3mz+1, big3mz+1)
  ## Transitions from never reproduced to never reproduced
  bigMReproAbsorbing[Z1, Z1] = bigMCond[Z1, Z1]
  ## Transitions from never reproduced to the absorbing state
  for (j in Z1)
    bigMReproAbsorbing[big3mz+1, j] = sum(bigMCond[Z23, j])
  ## Transitions from reproduced this year or sometime in the past to
  ## the absorbing state
  bigMReproAbsorbing[big3mz+1, Z23] = 1
  ## Transitions from the absorbing state to the absorbing state
  bigMReproAbsorbing[big3mz+1, big3mz+1] = 1

  bigMaAbsorbing = diag(1, big3mz+1)
  bigm0Absorbing = c(bigm0, 0)
  bigMaAbsorbingM0 = matrix (0, maxAge, big3mz+1)
  bigMaAbsorbingM0[1,] = bigm0Absorbing

  for (a in 2:maxAge) {
    bigMaAbsorbing = bigMReproAbsorbing %*% bigMaAbsorbing
    bigMaAbsorbingM0[a,] = bigMaAbsorbing %*% bigm0Absorbing
  }

  ## calculate big size hit time distribution
  reproArrivalAges = diff(bigMaAbsorbingM0[,big3mz+1])
  reproArrivalMode = which.max(reproArrivalAges)
  reproAgeDist = reproArrivalAges / sum(reproArrivalAges)

  return (out=list(reproArrivalMode=reproArrivalMode,
                   reproAgeDist=reproAgeDist,
                   q2Extended=q2Extended))
}
