## partitionVarSkewnessNoEnvVar ####################
P = matrix (c(0,0.3,0, 0,0.25,0.25, 0,0,0.5), 3, 3)
F = matrix (0, 3, 3); F[1,] = 0:2
c0 = c(1,0,0)
test_that("partitionVarSkewnessNoEnvVar works with Poisson dist. clutch sizes", {  
##  withr::local_options(width=20)
  expect_snapshot (
      partitionVarSkewnessNoEnvVar (P, F, c0)
##      waldo::partitionVarSkewnessNoEnvVar (P, F, c0)
  )
})

P = matrix (c(0,0.3,0, 0,0.25,0.25, 0,0,0.5), 3, 3)
F = matrix (0, 3, 3); F[1,] = 0.1*(0:2)
c0 = c(1,0,0)
test_that("partitionVarSkewnessNoEnvVar works with Bernoulli dist. clutch sizes", {  
  expect_snapshot (
      partitionVarSkewnessNoEnvVar (P, F, c0, maxAge=40,
                                    survThreshold=0.01,
                                    Fdist="Bernoulli")
  )
})

## partitionVarSkewnessEnvVar #######################

P1 = matrix (c(0,0.3,0, 0,0.25,0.25, 0,0,0.5), 3, 3)
P2 = matrix (c(0,0.4,0, 0,0.2,0.2, 0,0,0.6), 3, 3)
F1 = matrix (0, 3, 3); F1[1,] = 0:2
F2 = matrix (0, 3, 3); F1[1,] = 0.8*(0:2)
Plist = list (P1, P2)
Flist = list (F1, F2)
Q = matrix (1/2, 2, 2)
c0 = c(1,0,0)
test_that("partitionVarSkewnessEnvVar works with Poisson dist. clutch sizes", {
  expect_snapshot (
      partitionVarSkewnessEnvVar(Plist, Flist, Q, c0)
  )
})

P1 = matrix (c(0,0.3,0, 0,0.25,0.25, 0,0,0.5), 3, 3)
P2 = matrix (c(0,0.4,0, 0,0.2,0.2, 0,0,0.6), 3, 3)
F1 = matrix (0, 3, 3); F1[1,] = 0.1*(0:2)
F2 = matrix (0, 3, 3); F1[1,] = 0.2*(0:2)
Plist = list (P1, P2)
Flist = list (F1, F2)
Q = matrix (1/2, 2, 2)
c0 = c(1,0,0)
test_that("partitionVarSkewnessEnvVar works with Bernoulli dist. clutch sizes", {
  expect_snapshot (
      partitionVarSkewnessEnvVar (Plist, Flist, Q, c0, maxAge=60,
                                      survThreshold=0.01,
                                      Fdist="Bernoulli")
  )
})

## partitionVarSkewnessEnvVarAndTraits #######################
PlistAllTraits = FlistAllTraits = list ()
P1 = matrix (c(0,0.3,0, 0,0.25,0.25, 0,0,0.5), 3, 3)
P2 = 0.9*P1
PlistAllTraits[[1]] = list (P1, P2)
P1 = matrix (c(0,0.4,0, 0,0.2,0.2, 0,0,0.6), 3, 3)
P2 = 0.9*P1
PlistAllTraits[[2]] = list (P1, P2)
F1 = matrix (0, 3, 3); F1[1,] = 0:2
F2 = 0.8*F1
FlistAllTraits[[1]] = list (F1, F2)
F1 = matrix (0, 3, 3); F1[1,] = 0.2*(0:2)
F2 = 1.1*F1
FlistAllTraits[[2]] = list (F1, F2)
Q = matrix (1/2, 2, 2)
c0 = c(1,0,0)
traitDist = rep(0.5, 2)
test_that("partitionVarSkewnessEnvVarAndTraits works with Poisson-dist. clutch sizes", {
  expect_snapshot (
      partitionVarSkewnessEnvVarAndTraits (PlistAllTraits,
                                            FlistAllTraits, Q, c0,
                                           traitDist)
  )
})

PlistAllTraits = FlistAllTraits = list ()
P1 = matrix (c(0,0.3,0, 0,0.25,0.25, 0,0,0.5), 3, 3)
P2 = 0.9*P1
PlistAllTraits[[1]] = list (P1, P2)
P1 = matrix (c(0,0.4,0, 0,0.2,0.2, 0,0,0.6), 3, 3)
P2 = 0.9*P1
PlistAllTraits[[2]] = list (P1, P2)
F1 = matrix (0, 3, 3); F1[1,] = 0.1*(0:2)
F2 = 0.8*F1
FlistAllTraits[[1]] = list (F1, F2)
F1 = matrix (0, 3, 3); F1[1,] = 0.09*(0:2)
F2 = 1.1*F1
FlistAllTraits[[2]] = list (F1, F2)
Q = matrix (1/2, 2, 2)
c0 = c(1,0,0)
traitDist = rep(0.5, 2)
test_that("partitionVarSkewnessEnvVarAndTraits works with Bernoulli-dist. clutch sizes", {
  expect_snapshot (
      partitionVarSkewnessEnvVarAndTraits (PlistAllTraits,
                                           FlistAllTraits, Q, c0,
                                           traitDist,
                                           maxAge=60,
                                           survThreshold=0.05,
                                           Fdist="Bernoulli")
  )
})

## partitionVarSkewnessNoEnvVarPostBreeding ####################
P = matrix (c(0,0.3,0, 0,0.25,0.25, 0,0,0.5), 3, 3)
F = matrix (0, 3, 3); F[1,] = 0:2
c0 = c(1,0,0)
test_that("partitionVarSkewnessNoEnvVarPostBreeding works with Poisson dist. clutch sizes", {  
  expect_snapshot (
      partitionVarSkewnessNoEnvVar (P, F, c0)
  )
})

P = matrix (c(0,0.3,0, 0,0.25,0.25, 0,0,0.5), 3, 3)
F = matrix (0, 3, 3); F[1,] = 0.1*(0:2)
c0 = c(1,0,0)
test_that("partitionVarSkewnessNoEnvVarPostBreeding works with Bernoulli dist. clutch sizes", {  
  expect_snapshot (
      partitionVarSkewnessNoEnvVar (P, F, c0, maxAge=40,
                                    survThreshold=0.01,
                                    Fdist="Bernoulli")
  )
})


## Sanity check comparisons
## ###########################################

P = matrix (c(0,0.3,0, 0,0.25,0.25, 0,0,0.5), 3, 3)
F = matrix (0, 3, 3); F[1,] = 0:2
c0 = c(1,0,0)
Plist = list (P, P)
Flist = list (F, F)
PlistAllTraits = list (Plist, Plist)
FlistAllTraits = list (Flist, Flist)
Q = matrix (1/2, 2, 2)
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
