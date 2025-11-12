## partitionVarSkewnessNoEnvVar ####################
P = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
F = matrix (0, 3, 3); F[1,] = 0:2
c0 = c(1,0,0)
test_that("partitionVarSkewnessNoEnvVar works with Poisson dist. clutch sizes", {  
##  withr::local_options(width=20)
  expect_snapshot (
      partitionVarSkewnessNoEnvVar (P, F, c0)
##      waldo::partitionVarSkewnessNoEnvVar (P, F, c0)
  )
})


P = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
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

P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
P2 = matrix (c(0, 0.2, 0, 0, 0, 0.3, 0, 0, 0.4), 3, 3)
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

P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
P2 = matrix (c(0, 0.2, 0, 0, 0, 0.3, 0, 0, 0.4), 3, 3)
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
P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
P2 = 0.9*P1
PlistAllTraits[[1]] = list (P1, P2)
P1 = matrix (c(0, 0.5, 0, 0, 0, 0.2, 0, 0, 0.7), 3, 3)
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
P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
P2 = 0.9*P1
PlistAllTraits[[1]] = list (P1, P2)
P1 = matrix (c(0, 0.5, 0, 0, 0, 0.2, 0, 0, 0.7), 3, 3)
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
