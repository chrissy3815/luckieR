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

## probTraitCondLROPostBreeding #############################

P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
P2 = matrix (c(0, 0.5, 0, 0, 0, 0.2, 0, 0, 0.7), 3, 3)
Plist = list(P1, P2)
F1 = matrix (0, 3, 3); F1[1,] = 0:2
F2 = matrix (0, 3, 3); F2[1,] = 0.9*(0:2)
Flist = list (F1, F2)
c0 = c(1,0,0)
traitDist = rep(0.5, 2)

test_that("probTraitCondLROPostBreeding works with Poisson-dist. clutch sizes", {
  expect_snapshot (
      probTraitCondLROPostBreeding (PlistAllTraits, FlistAllTraits, Q,
                                    c0, maxClutchSize=10,
                                    maxLRO=15, traitDist)
  )
})

## probTraitCondLRONoEnv ####################################
P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
P2 = matrix (c(0, 0.5, 0, 0, 0, 0.2, 0, 0, 0.7), 3, 3)
Plist = list(P1, P2)
F1 = matrix (0, 3, 3); F1[1,] = 0:2
F2 = matrix (0, 3, 3); F2[1,] = 0.9*(0:2)
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

## probTraitCondLROPostBreedingNoEnv #################################
## BUG: This test passes if I run just this section in the console but
## does not pass if I run test_file("test-probTraitCondLRO.R").

P1 = matrix (c(0, 0.3, 0, 0, 0, 0.5, 0, 0, 0.5), 3, 3)
P2 = matrix (c(0, 0.5, 0, 0, 0, 0.2, 0, 0, 0.7), 3, 3)
Plist = list(P1, P2)
F1 = matrix (0, 3, 3); F1[1,] = 0:2
F2 = matrix (0, 3, 3); F2[1,] = 0.9*(0:2)
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

correct = probTraitCondLROPostBreeding (PlistAllTraits, FlistAllTraits, Q,
                                        c0, maxClutchSize=10,
                                        maxLRO=15,
                                        traitDist)

test_that("probTraitCondLROPostBreedingNoEnv gives the same answer as probTraitCondLROPostBreeding", {
  expect_equal (probTraitCondLROPostBreedingNoEnv (Plist, Flist, c0,
                                                   maxClutchSize=10,
                                                   maxLRO=15,
                                                   traitDist),
                correct)
})


