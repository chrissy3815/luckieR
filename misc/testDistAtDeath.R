## Finally works!!  12/8/25
rm(list=ls(all=TRUE)); graphics.off();
## Newborns mature to juveniles with prob. 0.5 and die otherwise.
## Juveniles mature to adults with prob. 0.5 and die otherwise.
## There is no adult survival.
## In the language of Kendall et al.'s "Persistent problems in the
## construction of matrix population models," Fig. 2, sigma_N = G_J =
## 0.5.  P_j = sigma_A = 0.
P1 = matrix (c(0,0.5,0, 0,0,0.5, 0,0,0), 3, 3)
## If you lived to be an adult 1, you have a mean of 1 offspring. 
## In the language of Kendall et al.'s "Persistent problems in the
## construction of matrix population models," Fig. 2, b_A = 0.8.
F1 = matrix (c(0,0,0, 0.4,0,0, 0,0,0), 3, 3)
## LRO dist.: You live to become an adult (and therefore attempt
## reproduction) with probability 0.5^2 = 0.25 and die before
## attempting reproduction with probability 1 - 0.25 = 0.75.  The LRO
## distribution is therefore a weight of 0.75 at R = 0 plus
## 0.25*Pois(lambda=1) = (1 - 0.25)*c(1, rep(0, 5)) + 0.25*dpois(0:5,
## 0.8)
Plist = list (P1, P1)
Flist = list (F1, F1)
Q = matrix (1/2, 2, 2)
##c0 = c(1,0)
c0 = c(1,0,0)
maxClutchSize = 5
maxAge=3
maxLRO = 5
Fdist = "Poisson"
source ("../R/distTraitCondR.R")
source ("../R/megamatrixFunctions.R")
source ("../R/partitioningUtilities.R")
source ("~/luckCalculations/Utilities.R")

mT = maxLRO + 1
mA = maxAge + 1

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

## Making a 2nd environment with the same demographics --- leave this
## complication aside for now.
if (TRUE) {  
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
  ## The columns of B should sum to 1 and they do.
  ## B = mk_B (F, maxClutchSize, Fdist)
  B = mk_BPostBreeding (M, F, maxClutchSize, Fdist)

  ## Construct A, the transition matrix for a (number of kids) x stage x
  ## env. model.
  out = make_AxTPostBreeding (B, M, mT)
  A = out$A
  mzA = bigmz*mT
  ##  K = out$K

  ## Create a survival probability vector 
  surv = apply (A, 2, sum);  die = 1-surv; 

  ## And the fundamental matrix of esA
  fundA <- solve(diag(ncol(A)) - A)

  ## Make omega
  omega = matrix(0, nrow(A), ncol(A)); 

  for (j in 1:ncol(fundA))
    omega[,j] = die*fundA[,j]

  distAtBirth = matrix(0, bigmz, mT)
  distAtBirth[,1] = m0  ## Everyone starts off with zero kids at the
  ## beginning of the time step

  ## Get distribution of states at death
  distAtDeath = omega %*% matrix(distAtBirth, ncol=1)  # state vector
  distAtDeath = matrix (distAtDeath, nrow=bigmz)
  distKidsAtDeath = apply (distAtDeath, 2, sum)
  distKidsAtDeath = distKidsAtDeath / sum(distKidsAtDeath)

  distKidsAtDeathLuckieR = calcDistLROPostBreeding (Plist, Flist, Q,
                                                    c0, maxClutchSize,
                                                    maxLRO, "Poisson")
}

  ## Without the fake 2nd environment #######
if (TRUE) {
  ## B[i,j] is the probability that a class-j individual has i-1 kids.
  ## We assume Poisson-distributed number of offspring.
  ## The columns of B should sum to 1 and they do.
  if (TRUE) {
    B = mk_BPostBreeding (P1, F1, maxClutchSize, Fdist)
  } else {
    ## Do this by hand, because the way we make A, B should be the
    ## clutch size distribution conditional on surviving to
    ## reproduction.  Yeah, I suck.  Still comes up with everyone
    ## having no kids.
    B = matrix (0, maxClutchSize+1, mz)
    B[1,1] = 1 ## stage 1 has no kids
    B[1,3] = 1 ## Neither does stage 3
    B[,2] = dpois(0:maxClutchSize, lambda=0.8)
  }

  ## Construct A, the transition matrix for a (number of kids) x stage x
  ## env. model.
  out = make_AxTPostBreeding (B, P1, mT)
  A = out$A
  mzA = mz*mT

  ## Create a survival probability vector
  surv = apply (A, 2, sum);  die = 1-surv;

  ## And the fundamental matrix of esA
  fundA <- solve(diag(ncol(A)) - A)

  ## Make omega
  omega = matrix(0, nrow(A), ncol(A));

  for (j in 1:ncol(fundA))
    omega[,j] = die*fundA[,j]

  distAtBirth = matrix(0, mz, mT)
  distAtBirth[,1] = c0  ## Everyone starts off with zero kids at the
  ## beginning of the time step

  ## Get distribution of states at death
  distAtDeath = omega %*% matrix(distAtBirth, ncol=1)  # state vector
  distAtDeath = matrix (distAtDeath, nrow=mz)
  ## We get the same answer that we do with 2 env.
  distKidsAtDeath2 = apply (distAtDeath, 2, sum)
  ## my sense of what distAtDeath should be.  Mine sums to 1 and so
  ## does distAtDeath, but distAtDeath has more mass at stage 3, zero
  ## kids.
  myDistAtDeath = matrix (0, 3, maxClutchSize+1)
  myDistAtDeath[1,] = 0.5*c(1, rep(0,5))
  myDistAtDeath[2,] = 0.25*c(1, rep(0,5))
    if (Fdist == "Poisson") {
      myDistAtDeath[3,] = 0.25*dpois(0:5, 0.8)
    } else if (Fdist == "Bernoulli") {
      myDistAtDeath[3,] = 0.25*dbinom(0:5, prob=0.8, size=1)
    }
  myDistKidsAtDeath = apply (myDistAtDeath, 2, sum)
}
