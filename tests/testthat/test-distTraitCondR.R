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
