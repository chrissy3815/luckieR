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
## Why doesn't the calculation work when I use bdiag?  No idea.
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

