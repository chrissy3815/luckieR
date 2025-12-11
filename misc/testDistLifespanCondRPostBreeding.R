## Finally works!!  12/8/25
rm(list=ls(all=TRUE)); graphics.off();
## Newborns mature to juveniles with prob. 0.5 and die otherwise.
## Juveniles mature to adults with prob. 0.5 and die otherwise.
## There is no adult survival.
## In the language of Kendall et al.'s "Persistent problems in the
## construction of matrix population models," Fig. 2, sigma_N = G_J =
## 0.5.  P_j = sigma_A = 0.
P = matrix (c(0,0.5,0, 0,0,0.5, 0,0,0), 3, 3)
## If you lived to be an adult 1, you have a mean of 1 offspring. 
## In the language of Kendall et al.'s "Persistent problems in the
## construction of matrix population models," Fig. 2, b_A = 0.8.
F = matrix (c(0,0,0, 0.4,0,0, 0,0,0), 3, 3)
## LRO dist.: You live to become an adult (and therefore attempt
## reproduction) with probability 0.5^2 = 0.25 and die before
## attempting reproduction with probability 1 - 0.25 = 0.75.  The LRO
## distribution is therefore a weight of 0.75 at R = 0 plus
## 0.25*Pois(lambda=1) = (1 - 0.25)*c(1, rep(0, 5)) + 0.25*dpois(0:5,
## 0.8)
Plist = list (P)
Flist = list (F)
##Q = matrix (1/2, 2, 2)
Q = matrix (1, 1, 1)
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

out = distLifespanCondR2PostBreeding (Plist, Flist, Q, c0,
                                      maxClutchSize, maxLRO, maxAge)


