## calcDistLRO ######################################
P1 = matrix (c(0,0.5, 0,0), 2, 2)
F1 = matrix (c(0,0, 0.8,0), 2, 2)
Plist = list (P1, P1)
Flist = list (F1, F1)
c0 = c(1,0)
Q = matrix (1/2, 2, 2)
maxClutchSize = 5
maxLRO=5
maxAge=3
correctLRODist = (1 - 0.5)*c(1, rep(0, maxLRO)) +
  0.5*dpois(0:maxLRO,0.8)
## Correct for truncating the clutch size
correctLRODist[maxLRO+1] = 1 - sum(correctLRODist[1:maxLRO])

test_that("calcDistLRO works with Poisson-dist. clutch sizes", {  
  expect_equal (
      calcDistLRO (Plist, Flist, Q, c0, maxClutchSize,maxLRO),
      correctLRODist
  )
})

correctLRODist = (1 - 0.5)*c(1, rep(0, maxClutchSize)) +
  0.5*dbinom(0:maxClutchSize, prob=0.8, size=1)
## Correct for truncating the clutch size
correctLRODist[maxClutchSize+1] = 1 - sum(correctLRODist[1:maxClutchSize])

test_that("calcDistLRO works with Bernoulli-dist. clutch sizes", {  
  expect_equal (
      calcDistLRO (Plist, Flist, Q, c0, maxClutchSize,
                               maxLRO, "Bernoulli"),
      correctLRODist
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

test_that("calcDistLRO works with Bernoulli-dist. clutch sizes", {  
  expect_equal (
      calcDistLROPostBreeding (Plist, Flist, Q, c0, maxClutchSize,
                               maxLRO, "Bernoulli"),
      correctLRODist
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
