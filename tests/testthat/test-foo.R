## BUG: This test passes if I just paste this section into the console but
## does not pass if I run test_file("test-probTraitCondLRO.R") (It
## says "object 'Plist' not found").  But if I do the test_file run
## after running by pasting the following into the console, all is
## well.  Why is Plist not recognized by test_file?

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
