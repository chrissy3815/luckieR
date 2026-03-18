## All juveniles become adults, who have exactly 1 kid and then die.
## There are two traits, which have identical parameters, and two
## environments, which are identical.
P = matrix (c(0,1,0,0), 2, 2)
F = matrix (c(0,0,1,0), 2, 2)
c0 = c(1,0)
Fdist = "Bernoulli"
maxLRO = 1
maxClutchSize = 1
Plist = list (P, P)
Flist = list (F, F)
PlistAllTraits = list (Plist, Plist)
FlistAllTraits = list (Flist, Flist)
Q = matrix (1/2, 2, 2)
traitDist = rep(0.5, 2)

out = partitionVarSkewnessEnvVarAndTraits (PlistAllTraits,
                                           FlistAllTraits,
                                           Q, c0, traitDist)

test_that("Traits contribute no variance", {
  expect_equal (out$varFromTraits, 0)
})

test_that("Traits contribute no skewness", {
  expect_equal (out$skewnessFromTraits[1,1], 0)
})

test_that("Environment contributes no variance", {
  expect_equal (sum(out$envTrajecVar), 0)
})

test_that("Traits contribute no skewness", {
  expect_equal (sum(out$envTrajecSkewness), 0)
})


## Sanity check comparisons #######################################

P = matrix (c(0,0.3,0, 0,0.25,0.25, 0,0,0.5), 3, 3)
F = matrix (0, 3, 3); F[1,] = 0:2
c0 = c(1,0,0)
Plist = list (P, P)
Flist = list (F, F)
PlistAllTraits = list (Plist, Plist)
FlistAllTraits = list (Flist, Flist)
Q = matrix (1/2, 2, 2)
traitDist = rep(0.5, 2)
out1 = partitionVarSkewnessNoEnvVar (P, F, c0)
out2 = partitionVarSkewnessEnvVar(Plist, Flist, Q, c0)
out3 = partitionVarSkewnessEnvVarAndTraits (PlistAllTraits,
                                            FlistAllTraits,
                                            Q, c0, traitDist)
test_that("Having the same matrices for two environments gives the same survival trajectory luck for variance as having a single environment", {
  expect_equal (out1$survTrajecVar, out2$survTrajecVar)
})

test_that("Having the same matrices for two environments gives the same growth trajectory luck for skewness as having a single environment", {
  expect_equal (out1$survTrajecSkewness, out2$survTrajecSkewness)
})

test_that("Having the same matrices for both traits give the same environment trajectory luck for variance as having a single trait", {
  expect_equal (out2$envTrajecVar, out3$envTrajecVar)
})

test_that("Having the same matrices for both traits give the same fecundity luck for skewness as having a single trait", {
  expect_equal (out2$fecSkewness, out3$fecSkewness)
})

## And the post-breeding versions

P = matrix (c(0,0.3,0, 0,0.25,0.25, 0,0,0.5), 3, 3)
F = matrix (0, 3, 3); F[1,] = 0:2
c0 = c(1,0,0)
Plist = list (P, P)
Flist = list (F, F)
PlistAllTraits = list (Plist, Plist)
FlistAllTraits = list (Flist, Flist)
Q = matrix (1/2, 2, 2)
out1 = partitionVarSkewnessNoEnvVarPostBreeding (P, F, c0)
out2 = partitionVarSkewnessEnvVarPostBreeding(Plist, Flist, Q, c0)
out3 = partitionVarSkewnessEnvVarAndTraits (PlistAllTraits,
                                            FlistAllTraits,
                                            Q, c0, traitDist,
                                            postBreedingCensus=TRUE)

test_that("Having the same matrices for two environments gives the same survival trajectory luck for variance as having a single environment for a post-breeding census", {
  expect_equal (out1$survTrajecVar, out2$survTrajecVar)
})

test_that("Having the same matrices for two environments gives the same growth trajectory luck for skewness as having a single environment for a post-breeding census", {
  expect_equal (out1$survTrajecSkewness, out2$survTrajecSkewness)
})

test_that("Having the same matrices for both traits give the same environment trajectory luck for variance as having a single trait for a post-breeding census", {
  expect_equal (out2$envTrajecVar, out3$envTrajecVar)
})

test_that("Having the same matrices for both traits give the same fecundity luck for skewness as having a single trait for a post-breeding census", {
  expect_equal (out2$fecSkewness, out3$fecSkewness)
})
