## distLifespanCondR2 ####################
P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
P2 = matrix (c(0, 0.2, 0, 0, 0, 0.3, 0, 0, 0.4), 3, 3)
F1 = matrix (0, 3, 3); F1[1,] = 0:2
F2 = matrix (0, 3, 3); F1[1,] = 0.8*(0:2)
Plist = list (P1, P2)
Flist = list (F1, F2)
Q = matrix (1/2, 2, 2)
c0 = c(1,0,0)
test_that("distLifespanCondR2 works with Poisson dist. clutch sizes", {  
  expect_snapshot (
      distLifespanCondR2 (Plist, Flist, Q, c0, maxClutchSize=12,
                          maxLRO=20,
                          maxAge=20) 
  )
})

P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
P2 = matrix (c(0, 0.2, 0, 0, 0, 0.3, 0, 0, 0.4), 3, 3)
Plist = list (P1, P2)
F1 = matrix (0, 3, 3); F1[1,] = 0.1*(0:2)
F2 = matrix (0, 3, 3); F1[1,] = 0.08*(0:2)
Flist = list (F1, F2)
Q = matrix (1/2, 2, 2)
c0 = c(1,0,0)
test_that("distLifespanCondR2 works with Bernoulli dist. clutch sizes", {  
  expect_snapshot (
      distLifespanCondR2 (Plist, Flist, Q, c0, maxClutchSize=1,
                          maxAge=12,
                          maxLRO=20,
                          percentileCutoff=0.95,
                          Fdist="Bernoulli")
  )
})

## distLifespanCondR2PostBreeding ##################
## If you live to stage 3, you produce a Poisson-dist. number of
## offspring in one clutch, then you die.  If you don't live to stage
## 3, you do not reproduce.  There is no survival in stages 1 or 2 ---
## you mature or you die.  Therefore, if you had any offspring at all,
## you must have lived 3 years, precisely.
P = matrix (c(0,0.5,0, 0,0,0.5, 0,0,0), 3, 3)
F = matrix (c(0,0,0, 0.4,0,0, 0,0,0), 3, 3)
Plist = list (P)
Flist = list (F)
Q = matrix (1, 1, 1)
c0 = c(1,0,0)
maxClutchSize = 5
maxAge=3
maxLRO = 5
Fdist = "Poisson"

out = distLifespanCondR2PostBreeding (Plist, Flist, Q, c0,
                                      maxClutchSize, maxLRO, maxAge)
test_that("distLifespandCondR2PostBreeding gives the correct conditional lifespan distribution for a simple 3-stage model", {
  expect_equal (out$probLifespanCondR[2:5, c(1,2,4)],
                matrix (0, nrow=4, ncol=3))
  expect_equal (out$probLifespanCondR[2:5, 3],
                rep(1, 4))
})

## calcDistLROPostBreeding ####################################
P1 = matrix (c(0,0.5,0, 0,0,0.5, 0,0,0), 3, 3)
F1 = matrix (c(0,0,0, 0.4,0,0, 0,0,0), 3, 3)
Plist = list (P1, P1)
Flist = list (F1, F1)
Q = matrix (1/2, 2, 2)
##c0 = c(1,0)
c0 = c(1,0,0)
maxClutchSize = 5
maxLRO=5
maxAge=3
correctLRODist = (1 - 0.25)*c(1, rep(0, maxLRO)) +
  0.25*dpois(0:maxLRO,0.8)
## Correct for truncating the clutch size
correctLRODist[maxLRO+1] = 1 - sum(correctLRODist[1:maxLRO])

test_that("calcDistLRO works with Poisson-dist. clutch sizes", {  
  expect_equal (
      calcDistLROPostBreeding (Plist, Flist, Q, c0, maxClutchSize,maxLRO),
      correctLRODist
  )
})

correctLRODist = (1 - 0.25)*c(1, rep(0, maxClutchSize)) +
  0.25*dbinom(0:maxClutchSize, prob=0.8, size=1)
## Correct for truncating the clutch size
correctLRODist[maxClutchSize+1] = 1 - sum(correctLRODist[1:maxClutchSize])

test_that("calcDistLRO works with Poisson-dist. clutch sizes", {  
  expect_equal (
      calcDistLROPostBreeding (Plist, Flist, Q, c0, maxClutchSize,
                               maxLRO, "Bernoulli"),
      correctLRODist
  )
})


## probTraitCondLRO ####################################
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
F1 = matrix (0, 3, 3); F1[1,] = 0.9*(0:2)
F2 = 1.1*F1
FlistAllTraits[[2]] = list (F1, F2)
Q = matrix (1/2, 2, 2)
c0 = c(1,0,0)
traitDist = rep(0.5, 2)
test_that("probTraitCondLRO works with Poisson-dist. clutch sizes", {  
  expect_snapshot (
      probTraitCondLRO (PlistAllTraits, FlistAllTraits, Q,
                        c0, maxClutchSize=10, maxLRO=15, traitDist)
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
F1 = matrix (0, 3, 3); F1[1,] = 0.2*(0:2)
F2 = 0.9*F1
FlistAllTraits[[2]] = list (F1, F2)
Q = matrix (1/2, 2, 2)
c0 = c(1,0,0)
traitDist = rep(0.5, 2)
test_that("probTraitCondLRO works with Bernoulli-dist. clutch sizes", {  
  expect_snapshot (
      probTraitCondLRO (PlistAllTraits, FlistAllTraits, Q,
      c0, maxClutchSize=1, maxLRO=15, traitDist, Fdist="Bernoulli")
  )
})

## distLifespanCondR2NoEnv ##################

P = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
F = matrix (0, 3, 3); F[1,] = 0:2
c0 = c(1,0,0)
maxClutchSize = 12
maxLRO = 30
maxAge=20

Plist = list(P, P)
Flist = list(F, F)
Q = matrix (1/2, 2, 2)

correct = distLifespanCondR2 (Plist, Flist, Q, c0, maxClutchSize,
                              maxLRO, maxAge)
test_that("distLifespanCondR2NoEnv gives the same answer as distLifespanCondR2", {
  expect_equal (distLifespanCondR2NoEnv (P, F, c0, maxClutchSize,
                                         maxLRO, maxAge),
                correct
                )
})

## distLifespanCondR2PostBreedingNoEnv ##################
## BUG: distLifespanCondR2PostBreeding is throwing warnings.  Re-write.
## You don't need those bullet matrices.
P = matrix (c(0,0.3,0, 0,0,0.5, 0,0,0.5), 3, 3)
F = matrix (0, 3, 3); F[1,] = c(0, 0.4, 0.3)
c0 = c(1,0,0)
maxClutchSize = 12
maxLRO = 30
maxAge=20

Plist = list(P, P)
Flist = list(F, F)
Q = matrix (1/2, 2, 2)

correct = distLifespanCondR2PostBreeding (Plist, Flist, Q, c0, maxClutchSize,
                                          maxLRO, maxAge)
test_that("distLifespanCondR2PostBreedingNoEnv gives the same answer as distLifespanCondR2PostBreeding", {
  expect_equal (distLifespanCondR2PostBreedingNoEnv (P, F, c0, maxClutchSize,
                                         maxLRO, maxAge),
                correct
                )
})

## calcDistLRONoEnv ####################################
P = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
F = matrix (0, 3, 3); F1[1,] = 0:2
c0 = c(1,0,0)
maxClutchSize = 10
maxLRO=20
Plist = list(P, P)
Flist = list(F, F)
Q = matrix (1/2, 2, 2)

correct = calcDistLRO (Plist, Flist, Q, c0, maxClutchSize, maxLRO)
test_that ("calcDistLRONoEnv gives the same answer as calcDistLRO", {
  expect_equal (calcDistLRONoEnv (P, F, c0, maxClutchSize, maxLRO),
                correct)
})

## calcDistLROPostBreedingNoEnv ####################################
P = matrix (c(0,0.5,0, 0,0,0.5, 0,0,0), 3, 3)
F = matrix (c(0,0,0, 0.4,0,0, 0,0,0), 3, 3)
c0 = c(1,0,0)
maxClutchSize = 10
maxLRO=20
Plist = list(P, P)
Flist = list(F, F)
Q = matrix (1/2, 2, 2)

correct = calcDistLROPostBreeding (Plist, Flist, Q, c0, maxClutchSize, maxLRO)
test_that ("calcDistLROPostBreedingNoEnv gives the same answer as calcDistLROPostBreeding", {
  expect_equal (calcDistLROPostBreedingNoEnv (P, F, c0, maxClutchSize, maxLRO),
                correct)
})

## probTraitCondLRONoEnv ####################################333
P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
P2 = matrix (c(0, 0.5, 0, 0, 0, 0.2, 0, 0, 0.7), 3, 3)
Plist = list(P1, P2)
F1 = matrix (0, 3, 3); F1[1,] = 0:2
F2 = matrix (0, 3, 3); F1[1,] = 0.9*(0:2)
Flist = list (F1, F2)
c0 = c(1,0,0)
traitDist = rep(0.5, 2)

P1list = list(Plist[[1]], Plist[[1]])
P2list = list(Plist[[2]], Plist[[2]])
PlistAllTraits = list(P1list, P2list)

F1list = list(Flist[[1]], Flist[[1]])
F2list = list(Flist[[2]], Flist[[2]])
FlistAllTraits = list (F1list, F2list)

Q = matrix (1/2, 2, 2)

correct = probTraitCondLRO (PlistAllTraits, FlistAllTraits, Q,
                            c0, maxClutchSize=10, maxLRO=15, traitDist)

test_that("probTraitCondLRONoEnv gives the same answer as probTraitCondLRO", {
  expect_equal (probTraitCondLRONoEnv (Plist, Flist, c0,
                                       maxClutchSize=10, maxLRO=15, traitDist), correct)
})
