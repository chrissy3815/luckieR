## makeCondKernel #####################################

## All (stage 1) mature into stage 2 which matures into stage 3 with
## no losses.  Stage 3 is absorbing.  Verify that the kernel
## conditional on reaching stage 3 is simply the unconditional kernel.
M = matrix (c(0,1,0, 0,0,1, 0,0,0.5), 3, 3)
transientStates = c(1,2)
out = makeCondKernel (M, transientStates)
test_that("makeCondKernel returns M when reaching the absorbing states is guaranteed", {
  expect_equal (out$MCond, M)
})

## makeCondFailureKernel ##################################

## All stage 1 mature into stage 2 which never mature into stage 3.
## Stage 3 is absorbing.  Verify that the kernel conditional on never
## reaching stage 3 is simply the unconditional transient kernel.
M = matrix (c(0,1,0, 0,0.8,0, 0,0,0.5), 3, 3)
transientStates = c(1,2)
out = makeCondFailureKernel (M, transientStates)
test_that("makeCondFailureKernel returns M when reaching the absorbing states is guaranteed not to occur", {
  expect_equal (out$MCond, M[transientStates, transientStates])
})

## calcMoments #################################################

## All juveniles become adults, who have exactly f kids and then die.
f = 2
M = matrix (c(0,1,0,0), 2, 2)
F = matrix (c(0,0,f,0), 2, 2)
## Clutch size is Bernoulli distributed
R1 = matrix (0, 3, 3)
for (j in 1:2)
  R1[,j] = sum(F[,j])
R2 = R1^2
R3 = R1^3
out = calcMoments (M, R1, R2, R3)
test_that ("Mean LRO is 2", {
  expect_equal (out$rho1Vec[1,], c(f,f))
})

## All juveniles become adults, who have exactly f kids and then die.
M = matrix (c(0,1,0,0), 2, 2)
F = matrix (c(0,0,f,0), 2, 2)
## Clutch size is Bernoulli distributed
R1 = matrix (0, 3, 3)
for (j in 1:2)
  R1[,j] = sum(F[,j])
R2 = R1^2
R3 = R1^3
out = calcMoments (M, R1, R2, R3)
test_that ("Variance of LRO is 0", {
  expect_equal (out$mu2Vec[1,], c(0,0))
})

## All juveniles become adults, who have exactly f kids and then die.
M = matrix (c(0,1,0,0), 2, 2)
F = matrix (c(0,0,f,0), 2, 2)
## Clutch size is Bernoulli distributed
R1 = matrix (0, 3, 3)
for (j in 1:2)
  R1[,j] = sum(F[,j])
R2 = R1^2
R3 = R1^3
out = calcMoments (M, R1, R2, R3)
test_that ("Skewness of LRO is 0", {
  expect_equal (out$skewnessVec, c(0,0))
})

## makeMCondLROThreshold #################################################
## All juveniles become adults, who have exactly 1 kid and then die.
M = matrix (c(0,1,0,0), 2, 2)
F = matrix (c(0,0,1,0), 2, 2)
m0 = c(1,0)
Fdist = "Bernoulli"
maxLRO = 1
maxClutchSize = 1

mT = maxLRO+1
mzA = 2*mT
bigmz = 2
## Fbullet updates number of kids.
Fbullet = matrix (0, mzA, mzA)

for (j in 1:mT) {
  for (i in j:(mT-1)) {
    for (z in 1:bigmz) {
      if (Fdist == "Poisson") {
        Fbullet[(i-1)*bigmz + z, (j-1)*bigmz + z] =
          dpois (i-j, lambda=sum(Fmat[,z]))
      } else if (Fdist == "Bernoulli") {
        Fbullet[(i-1)*bigmz + z, (j-1)*bigmz + z] =
          dbinom (i-j, prob=sum(Fmat[,z]), size=1)
      } else {
        stop ("Supported options for Fdist are 'Poisson' and 'Bernoulli.'\n")
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
      1 - sum(Fbullet[(0:(mT-2))*bigmz + z, (j-1)*bigmz + z])
  }
}
## If you have #kids index = mT, you stay there.
for (z in 1:bigmz)
  Fbullet[(mT-1)*bigmz + z, (mT-1)*bigmz + z] = 1

## Mbullet updates everything else
Mbullet = matrix (0, mzA, mzA)
for (k in 1:mT)
  Mbullet[(k-1)*bigmz + 1:bigmz, (k-1)*bigmz + 1:bigmz] = M
esmzA = 2*mzA
esA = matrix (0, esmzA, esmzA)
esA[mzA + 1:mzA, 1:mzA] = Fbullet
esA[1:mzA, mzA + 1:mzA] = Mbullet

out = makeMCondLROThreshold (M, F, threshold, m0, maxLRO,
                             maxClutchSize, Fdist)
test_that("makeMCondLROThreshold returns the unconditional kernel (in the extended state space) when everyone is guaranteed to make the LRO threshold",{
  expect_equal (out$ACondSucceed, esA)
})

## Less trivial model, set the threshold to 1 and make sure the result
## agrees with makePCondBreedDef3
M = matrix (c(0.25, 0.3, 0, 0.65), 2, 2)
F = matrix (c(0,0,0.8,0), 2, 2)
m0 = c(1,0)
Fdist = "Poisson"
maxLRO = 8
maxClutchSize = 8
threshold = 1
out1 = makeMCondLROThreshold (M, F, threshold, m0, maxLRO,
                              maxClutchSize, Fdist)
probSucceed = out1$probSucceed
out2 = makePCondBreedDef3 (M, F, m0)
probBreed = out2$pM

test_that("makeMCondLROThreshold with a threshold of LRO = 1 (do you breed at least once?) gives the same probability of success as makePCondBreedDef3",{
  expect_equal (probSucceed, probBreed)
})

## makeMCondLROThresholdPostBreeding #######################
## All juveniles become adults, who have 1 kid and then die
P = matrix (c(0,1,0, 0,0,1, 0,0,0), 3, 3)
F = matrix (c(0,0,0, 1,0,0, 0,0,0), 3, 3)
m0 = c(1,0,0)
threshold = 1
maxLRO = 1
maxClutchSize = 1

out = makeMCondLROThresholdPostBreeding (P, F, threshold, m0, maxLRO=maxLRO,
                                         maxClutchSize=maxClutchSize,
                                         Fdist="Bernoulli")
probSucceed = out$probSucceed
ACondSucceed = out$ACondSucceed

B = mk_B (F, maxClutchSize, Fdist="Bernoulli")
## Construct A, the transition matrix for a (number of kids) x stage x
## env. model.
mT = maxLRO + 1
out = make_AxT (B, P, mT)
A = out$A

test_that("makeMCondLROThresholdPostBreeding returns the unconditional kernel if success is guaranteed",{
  expect_equal (ACondSucceed, A)
})

