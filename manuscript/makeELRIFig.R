rm(list=ls(all=TRUE)); graphics.off();

library ("luckieR")

species = "ELRI"

colorVecAll6 = c("mediumpurple1", "black", "orange", "forest green",
                 "turquoise", "red")
colorVec = colorVecAll6[2:5]

numTraits=2
numEnv = 13
maxAge = 15
maxLRO = 16

fname = paste0 (species, "MaxLRO", maxLRO, ".rda")
load (fname)

R25 = min(which(cumsum(probR) >= 0.25))
R50 = min(which(cumsum(probR) >= 0.5))
R75 = min(which(cumsum(probR) >= 0.75))
R90 = min(which(cumsum(probR) >= 0.9))
R95 = min(which(cumsum(probR) >= 0.95))
R99 = min(which(cumsum(probR) >= 0.99))

dev.new ()
centralPanel = matrix(1, 7, 6)
horizPanel = matrix (2, 1, 6)
mat = rbind(centralPanel, horizPanel)
layout (mat=mat)
par(mgp=c(2,1,0), cex.axis=1.6, cex.lab=1.8, cex=1.1,
    cex.main=1.8, bty="n");
par (mar=c(5,6,3,2), mgp=c(4,1,0))

matplot (0:(R99-1), probXCondR[1:R99,], xlab="LRO",
         ylab="Probability", type="l", col=c("black", "red"), lwd=2,
         lty=1, ylim=c(0,1))
legend (x="topleft", col=c("black", "red"),
        legend=c("With endophytes", "Without endophytes"),
        lwd=2, bty="n", cex=1.8)

par (mar=c(2,5,0,1))  
plot (0:(R99-1), probR[1:R99], axes=FALSE, type="l", lwd=2,
      xlab="", ylab="")
axis (side=1, at=c(R75, R90, R95, R99),
      label=c("75th%", "90th%", "95th%", "99th%"))

dev.copy2pdf (file="probXCondR.pdf")

load ("grassVarSkewnessResults.rds")

lifespans = apply (lifespanCondX, 1, max)

dev.new ()
par (bty="n", cex.axis=1.6, cex.lab=1.8, cex.main=1.8,
     mar=c(5.1,4.5,4.1,2.1))

load ("modalTimeToFirstRepro.rds")

index=2
ymax = max(varTrajecResults[index,1,,])
matplot (1:lifespans[index],
         t(varTrajecResults[index,1,,1:lifespans[index]]),
         xlab="Age", ylab="Contrib. to Var(R)", type="l",
         lwd=2, col=colorVec, lty=1, ylim=c(0, ymax))
lines (rep(modalTimeToFirstReproCondX[index, 1], 2),
       c(0, ymax), col="purple", lty=2, lwd=3)

legend (x="topright", lwd=2, col=c(colorVec, "purple"),
        legend=c("Survival", "Growth", "Environment",
                 "Fecundity", "Modal age at 1st repro."),
        lty=c(rep(1, 4), 2), cex=1.4, bty="n")

dev.copy2pdf (file="luckTrajecELRI.pdf")


maxClutchSize = 24
maxLRO = 24
groupDist = rep (0.5, 2)

## read in matrices
allFecMatrices = readRDS ("GrassEndo_list_of_mean_matrices.rds")
allPMatrices = readRDS ("GrassEndo_mean_list_of_Pmatrices.rds")

## Work with endophyte-positive individuals for this example.
rightEndoFec = allFecMatrices$Eplus
rightEndoP = allPMatrices$Eplus
rightEndoFec = rightEndoFec$Elymus_villosus
rightEndoP = rightEndoP$Elymus_villosus

## find number of year types
numEnv = length(PlistAllEndo[[1]])

## Assume that we get random sequences of year types, drawing with
## replacement
## List of matrices for each trait value (endophyte status)
Flist = Plist = list ()

for (q in 1:numEnv) {
  ## Strip the $iter1 silliness off of the entries of Plist and Flist
  Flist[[q]] = rightEndoFec[[q]]$iter1
  Plist[[q]] = rightEndoP[[q]]$iter1
}

## Assume that all env. are equally likely and the sequence is
## temporally uncorrelated.  You can make other choices if you have
## more information.
Q = matrix (1/numEnv, numEnv, numEnv)

## How many states do we have in each env.?
mz = dim(Plist[[1]])[1]

bigmz = mz * numEnv

## Initial state distrib.
c0 = rep (0, mz); c0[1] = 1
## Initial env. distrib.  Should probably be the stationary distribution.
u0 = rep (1/numEnv, numEnv)
## Initial cross-classified state (size x env.)
m0 = matrix (outer (c0, as.vector(u0)), bigmz, 1)

## Maximum age
maxAge = 15

aveSizeCondBreedAndSurv = aveSizeCondNeverBreedAndSurv =
  fracSurvCondBreed = fracSurvCondNeverBreed = rep (0, maxAge+1)

## Calculate fecundity and size/growth transition matrices for
## cross-classified states (size x env.)
F = M = matrix (0, bigmz, bigmz)
for (i in 1:numEnv) {
  for (j in 1:numEnv) {
    M[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Plist[[j]]*Q[i,j]
    F[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Flist[[j]]*Q[i,j]
  }
}

## Find M conditional on having at least one offspring
out = makePCondBreedDef3 (M, F, m0)

## Kernel conditional on breeding
MCondBreed = out$PCondBreed
## Kernel conditional on never breeding
MCondNeverBreed = out$PCondNeverBreed
## Initial state conditional on breeding
esM0CondBreed = out$bigc0CondBreed
## Initial state conditional on neverbreeding
esM0CondNeverBreed = out$bigc0CondNeverBreed

## State distrib. at age a. (Here a = 0.)
MaM0CondBreed = esM0CondBreed
MaM0CondNeverBreed = esM0CondNeverBreed

sizeVec = rep(0:(mz-1), numEnv)
esSizeVec = rep(sizeVec, 3)
## Loop over ages
for (a in 0:maxAge) {
  ## What fraction of the population is alive?
  fracSurvCondBreed[a+1] = sum(MaM0CondBreed)
  fracSurvCondNeverBreed[a+1] = sum(MaM0CondNeverBreed)
  ## What is the average size?
  aveSizeCondBreedAndSurv[a+1] = sum(MaM0CondBreed * esSizeVec) /
    fracSurvCondBreed[a+1]
  aveSizeCondNeverBreedAndSurv[a+1] = sum(MaM0CondNeverBreed * sizeVec) /
    fracSurvCondNeverBreed[a+1]
  ## Update state distribution
  MaM0CondBreed = MCondBreed %*% MaM0CondBreed
  MaM0CondNeverBreed = MCondNeverBreed %*% MaM0CondNeverBreed
}

## which age is closest to having 10% of the population left?
fewLeftCondBreed = which.min (abs (fracSurvCondBreed - 0.1))
fewLeftCondNeverBreed = which.min (abs (fracSurvCondNeverBreed -
                                             0.1))

xmax = max(fewLeftCondNeverBreed, fewLeftCondBreed)
dev.new (height=6, width=12)
par (cex.lab=1.9, cex.axis=1.8, lwd=2, mfrow=c(1,2),
     oma=c(0,0,1,0), mar=c(5.1,4.5,4.1,2.1), bty="n")
plot (0:(xmax-1), fracSurvCondNeverBreed[1:xmax], type="l", xlab="Age",
      ylab="Survival to age a", xaxt="n")
lines (0:(xmax-1), fracSurvCondBreed[1:xmax], col="red")
axis (side=1, at=0:(xmax-1), labels=0:(xmax-1))
legend (x="topright", col=c("black", "red"),
        legend=c("Never breed", "Breed"), lwd=2, bty="n",
        cex=1.8)

ymax = max (aveSizeCondNeverBreedAndSurv, aveSizeCondBreedAndSurv)
plot (0:(xmax - 1), aveSizeCondNeverBreedAndSurv[1:xmax],
      type="l", xlab="Age", xaxt="n",
      ylab="Ave. size if survive",
      ylim=c(0, ymax))
axis (side=1, at=0:(xmax-1), labels=0:(xmax-1))
lines (0:(xmax-1), aveSizeCondBreedAndSurv[1:xmax],
       type="l", col="red")

dev.copy2pdf (file="condKernelPlots.pdf")
