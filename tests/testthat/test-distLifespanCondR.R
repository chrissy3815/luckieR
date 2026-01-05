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
