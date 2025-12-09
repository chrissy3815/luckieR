##################################################################
## R translation of formulas in 
##     Cochran, M. E. and S. Ellner.  1992.  Simple methods for 
##     calculating age-based life history parameters for stage-
##     structured populations.  Ecological Monographs 62:345-364.
## It is intended as a replacement for the StageCoach f77 program 
## that was originally provided as an online supplement to 
## Cochran & Ellner (1992). 
##
## Most original coding by Erin Souder, some edits by S.P. Ellner 
##################################################################

## Matrices are named as follows 
##    A is the population projection matrix (sexual and clonal births + state transitions) 
##    B is the sexual births matrix 
##    F is the Fission matrix 
##    P is the state transition ("survival/growth") matrix, not including splitting 
##    C = F + P, transition matrix where splitting is viewed as a transition 
##
##   Relation to Com(p)adre notation: 
##      StageCoach A = Com(p)adre matA
##    	StageCoach B = Com(p)adre matF
##      StageCoach P = Com(p)adre matU 
##      StageCoach F = Com(p)adre matC 
##      StageCoach C = Com(p)adre matU + matC 
##       
##	  We have A = P + B + F = matU + matF + matC 	  

## setwd("<the path to the folder containing luckieR")
if (Sys.info()["user"] == "Ellner") setwd("c:/repos") 

setwd("luckieR/stagecoach")

library(readxl)
library(Matrix)
library(matrixcalc)
library(expm)

Caswell_A <- read_excel("Caswell_A.xlsx",sheet = 1)
Caswell_P <- read_excel("Caswell_P.xlsx")
Caswell_B <- read_excel("Caswell_B.xlsx")
Fission <- read_excel("Caswell_F.xlsx")

A <- as.matrix(Caswell_P + Caswell_B + Fission)
B <- as.matrix(Caswell_B)
C <- as.matrix(Caswell_P + Fission)
P <- as.matrix(Caswell_P)

#################################################################
##                   Stage Based Information                   ##
#################################################################

#################################################################
##  Function to calculate the population growth rate           ##
##  Argument/s:                                                ##
##  A: Population projection matrix A                          ##
##  Return Value/s:                                            ##
##  Population growth rate for matrix A                        ##
##  Author/s:                                                  ##
##  Erin Souder                                                ##
##  Date: 03/01/2023                                           ##
#################################################################
pop_growth <- function(A) {
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  
   results <- Re(eigen(A)$values[1]) # need first real eigenvalue
   return(results) 
}

##################################################################
##  Function to calculate the stable stage distribution          ##
##  Argument/s:                                                 ##
##      A: Transition matrix A                                  ##
##  Return Value/s:                                             ##
##      Stable stage distribution                               ##
##  Author/s:                                                   ##
##  Erin Souder                                                 ##
##  Date: 03/01/2023                                            ##
  ##################################################################
stable_stage_dist <- function(A){
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  
  allvectors <- eigen(A) 
  num <- Re(allvectors$vectors[,1]) # need the first column of eigenvectors
  results <- num/sum(num)           # scaled eigenvector
  return(results) 
}

#################################################################
##  Function to compute the reproductive value vector v        ##
##  scaled so that v[1] = 1                                    ##
##  Argument/s:                                                ##
##      A: Transition matrix A                                 ##
##  Return Value/s:                                            ##
##      Relative reproductive value                                     ##
##  Author/s:                                                  ##
##  Erin Souder                                                ##
##  Date: 03/01/2023                                           ##
#################################################################
reproductive_value <- function(A){
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  
  num <- Re(eigen(t(A))$vectors[,1]) # transpose of matrix A
   # Column one of eigenvectors
  results = num/num[1]; 
  return(results) 
}

#################################################################
##  Function to determine sensitivity matrix                   ##
##          note how the function name is spelled              ##
##  Argument/s:                                                ##
##      A: Projection matrix A                                 ##
##  Return Value/s:                                            ##
##     Sensitivity matrix                                      ##
##  Author/s:                                                  ##
##  Erin Souder                                                ##
##  Date: 03/01/2023                                           ##
#################################################################
sensitivy_mat <- function(A){

  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  
  num <- reproductive_value(A) %*% t(stable_stage_dist(A)) 
  #reproductive value multiplied by transpose of stable stage distribution
  den <- as.numeric(reproductive_value(A) %*% stable_stage_dist(A))
  results <- Re(num/den)
  results
}

#################################################################
##  Function to calculate elasticity matrix                    ##
##  Argument/s:                                                ##
##      A: Projection matrix A                                 ##
##  Return Value/s:                                            ##
##      Elasticity matrix                                      ##
##  Author/s:                                                  ##
##  Erin Souder                                                ##
##  Date: 03/01/2023                                           ##
#################################################################
elasticity_mat <- function(A){
  

  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  
  x <- 1/pop_growth(A) 
  y <- sensitivy_mat(A)
  results <- (x * y) * A # sensitivity matrix multiplied by 1/growth rate
  Re(results)
}

##################################################################
##                END OF STAGE BASED INFORMATION                ##
##################################################################



#################################################################
##                BEGINNING OF AGE BASED INFORMATION           ##
#################################################################

######################################################################################
##  Function to calculate stage frequency of newborns at stable stage distribution  ##
##  Equation 19                                                                     ##
##  Argument/s:                                                                     ##
##      A: Transition matrix A                                                      ##
##      B: birth matrix B                                                           ##
##  Return Value/s:                                                                 ##
##      frequency of newborns at stable stage distribution for each newbornType     ##
##  Author/s:                                                                       ##
##  Erin Souder                                                                     ##
##  Date: 03/01/2023                                                                ##
######################################################################################
n_bj <- function(A,B){
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  
  num <- B %*% stable_stage_dist(A) # need to use matrix multiplication
  den <- sum(num)
  results <-  num/den
  results
}

#################################################################
##  Function to calculate the average age in stage i           ##
##  Equation 23 in Cochran & Ellner (1992)                     ##
##  Argument/s:                                                ##
##     A: Transition matrix A                                  ##
##     B: Birth matrix B                                       ##
##     C: Survival matrix P plus fission matrix F              ##
##  Return Value/s:                                            ##
##     average age of individuals in each stage                ##
##  Author/s:                                                  ##
##  Erin Souder                                                ##
##  Date: 03/01/2023                                           ##
#################################################################
age_in_stage <- function(A, B, C){

  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
  Imat <- diag(dim(A)[1]) # Identity matrix
  lambda <- pop_growth(A)
  bj <- n_bj(A,B)
  num <- rowSums(solve((Imat - (C/lambda)) %^% 2) %*% bj) 
  den <- rowSums(solve(Imat - (C/lambda)) %*% bj)
 
  results <- num/den
  results
}

##################################################################################
##  Function to calculate the standard deviation of age in stage i              ##
##  Equation 24 in Cochran and Ellner 1992                                      ##
##  Argument/s:                                                                 ##
##      A: Transition matrix A                                                  ##
##      B: Birth matrix B                                                       ##
##      C: Survival matrix P + fission matrix C                                 ##
##  Return Value/s:                                                             ##
##       Standard deviation of age in each stage                                ##
##  Author/s:                                                                   ##
##  Erin Souder                                                                 ##
##  Date: 03/01/2023                                                            ##
##################################################################################
age_in_stage_SD <- function(A, B, C){
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
  Imat <- diag(dim(A)[1]) # Identity matrix
  lam = pop_growth(A); 
  mat1 = solve(Imat + C/lam); 
  mat2 = solve(Imat - C/lam); 
  wts = n_bj(A,B); 
  
  num <- rowSums((mat1%*%mat2%*%mat2%*%mat2) %*% wts)
  den <- rowSums(mat2 %*% wts)
  
  results <- sqrt((num/den) - (age_in_stage(A,B,C))^2)
  return(results)
}

  #################################################################
  ##  Function to calculate fraction of age-a individuals in     ##
  ##      stage i, p_{i,a}                                       ##    
  ##  Equation 22 in Cochran and Ellner 1992                     ##
  ##  Argument/s:                                                ##
  ##  A: Transition matrix A                                     ##
  ##  B: Birth matrix B                                          ##
  ##  C: Survival matrix P + Fission matrix F                    ##
  ##  maxAge: number of years to use default 50                  ##
  ##  Return Value/s:                                            ##
  ##  Fraction of age-a individuals in stage class i             ##
  ##  Author/s:                                                  ##
  ##  Dr. Stephen Ellner                                         ##
  ##  Dr. Simone Blomberg                                        ##
  ##  Erin Souder                                                ##
  ##  Date: 03/01/2023                                           ##
  #################################################################
pit <- function(A,B,C, maxAge = 50){

  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}

  lam <- pop_growth(A)
  #lambda
  b <- n_bj(A,B)
  Imat <- diag(x = 1, dim(C))
  # Identity matrix
  den <- rowSums(solve(Imat - (C/lam)) %*% b)
  # 
  results <- matrix(NA, nrow = maxAge + 1, ncol = nrow(A))
  Ca = Imat;
  for (a in 0 : (maxAge)) {
   num <- (lam^(-a)) %*% rowSums((Ca) %*% b)
   x <- num/den
   results[a + 1,] <- x
   Ca <- C %*% Ca
  }
  return(results) 
}

########################################################################################
##  Function to calculate the annual number of newborns generated by an individual    ##
##      in state i, weighted by relative reproductive value                           ##
##  Equation 12 in Cochran and Ellner (1992)                                          ##
##  Argument/s:                                                                       ##
##      A: Transition matrix A                                                        ##
##      B: Birth matrix B                                                             ##
##  Return Value/s:                                                                   ##
##      Number of newborns generated by individuals in stage i, weighted              ##
##  Author/s:                                                                         ##
##  Erin Souder                                                                       ##
##  Date: 03/01/2023                                                                  ##
########################################################################################
gam_i <- function(A,B){

  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  
   scaled_rep_value <- function(A){
    # Scaling the reproductive value so the first one is 1
    x <- reproductive_value(A)[1]
    results <- reproductive_value(A)/x
    results
  }
  
  v <- scaled_rep_value(A)
  results <- t(B) %*% v
  results
}


bet_i <- function(B){
#######################################################################################
##  Function to calculate the annual number of newborns generated by an individual   ##
##      in stage i, NOT weighted by relative reproductive value                      ##
##  Argument/s:                                                                      ##
##     B: Birth matrix B                                                             ##
##  Return Value/s:                                                                  ##
##     Average total number of offspring for individuals in stage i, unweighted      ##
##  Author/s:                                                                        ##
##  Erin Souder                                                                      ##
##  Date: 03/01/2023                                                                 ##
####################################################################################### 
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  
  results <- colSums(B);  
  return(results)
}

##########################################################################
##  Generation time measure \bar{A}: mean age of parents for all        ##
##    newborns in one year, for a stable population                     ##  
##  Equation 26 in Cochran and Ellner 1992                              ##
##  Argument/s:                                                         ##
##    A: Transition matrix A                                            ##
##    B: Birth matrix B                                                 ##
##    C: Survival matrix P + Fission matrix F                           ##
##    weighted: if TRUE, offspring are weighted by reproductive value   ##
##  Return Value/s:                                                     ##
##  Average age of parents produced in current time frame               ##
##  Author/s:                                                           ##
##  Erin Souder                                                         ## 
##  Steve Ellner ('weighted' or not option)                             ##
##  Date: 03/01/2023                                                    ##
##########################################################################
pop_gen_time <- function(A, B, C,weighted=FALSE){

  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")} 
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
  if(weighted) kids_vector = gam_i(A,B) 
  if(!weighted) kids_vector = bet_i(A,B) 
  num <- sum(age_in_stage(A,B,C) * stable_stage_dist(A) * kids_vector) 
  den <- sum(stable_stage_dist(A) * kids_vector)
  results <- num/den
  return(results)
}

 
stable_age_distribution <- function(A,B,C, maxAge = 10){
  
  #################################################################
  ##  Function to determine stable age distribution              ##
  ##  Equation 31                                                ##
  ##  Argument/s:                                                ##
  ##  A: Transition matrix A                                     ##
  ##  B: Birth matrix B                                          ##
  ##  C: Survival matrix P + Fission matrix F                    ##
  ##  maxAge: number of years to use default 10                  ##
  ##  Return Value/s:                                            ##
  ##  Stable age distribution                                    ##
  ##  Author/s:                                                  ##
  ##  Dr. Stephen Ellner                                         ##
  ##  Dr. Simone Blomberg                                        ##
  ##  Erin Souder                                                ##
  ##  Date: 03/01/2023                                           ##
  #################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
  pp <- pit(A,B,C,maxAge)
  w <- stable_stage_dist(A)
 
  num <- pp %*% w
  den <- sum(w)
  results <- num/den
  results
}


life_expectancy <- function(P){
    
  ##################################################################
  ##  Function to determine life expectancy for each stage        ##
  ##  Equation 3                                                  ##
  ##    Argument/s:                                               ##
  ##    P: Survival matrix P                                      ##
  ##       OR                                                     ##
  ##     C: Survival matrix P + Fission matrix F                  ##
  ##  Return Value/s:                                             ##
  ##     Remaining Life expectancy for each stage                 ##
  ##  Author/s:                                                   ##
  ##  Erin Souder                                                 ##
  ##  Date: 03/01/2023                                            ##
  ##################################################################
  
  if(!is.matrix(P)) {stop("This is not a matrix")}
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  
  Imat <- diag(dim(P)[1]) # Identity matrix
  
  results <- colSums(solve(Imat - P)) # take the inverse of I-c in order to determine life expectancy
  results
}


life_expectancy_SD <- function(P){
    
  #######################################################################################
  ##  Function to determine the standard deviation for life expectancy for each stage  ##
  ##  Equation 5 in Cochran and Ellner 1992                                            ##
  ##  Argument/s:                                                                      ##
  ##  P: Survival matrix P                                                             ##
  ##  OR                                                                               ##
  ##  C: Survial matrix P + Fission matrix F                                           ##
  ##  Return Value/s:                                                                  ##
  ##  Standard deviation for the life expectancy for each stage                        ##
  ##  Author/s:                                                                        ##
  ##  Erin Souder                                                                      ##
  ##  Date: 03/01/2023                                                                 ##
  #######################################################################################
  
  if(!is.matrix(P)) {stop("This is not a matrix")}
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  
  Imat <- diag(dim(P)[1]) # Identity matrix
  #take the inverse and multiply it through twice and multiply by the sum of (I + C)
  y <- colSums((Imat + P) %*% (solve(Imat - P) %^% 2))
  #take the square root of (this answer minus life expectancy squared) for the SD
  results <- sqrt(y - life_expectancy(P)^2) 
  results
}



Di_mat <- function(P){

  ########################################################################
  ##  Function to calculate modified transition matrix D_i for stage i, ##
  ##      where individuals entering stage i stay there for one time    ##
  ##      step and then die (100% mortality)                            ##
  ##  Equation 8 in Cochran and Ellner 1992                             ##
  ##  Argument/s:                                                       ##
  ##  P: Survival matrix P                                              ##
  ##  Return Value/s:                                                   ##
  ##     Array of D of transition matrices, such that D[,,i] is D_i     ##
  ##  Author/s:                                                         ##
  ##  Dr. Stephen Ellner                                                ##
  ##  Erin Souder                                                       ##
  ##  Date: 03/01/2023                                                  ##
  ########################################################################
  
  if(!is.matrix(P)) {stop("This is not a matrix")}
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  
  results <- array(0, dim = c(dim(P),ncol(P)))
  for (i  in 1:ncol(P)) {
    results[,,i] = P; 
    results[,i,i] = 0; 
  }
  return(results) 
}
 
meantime <- function(P){
  ############################################################################
  ##  Function to compute the average time to reach stage i from stage j,   ##
  ##         conditional on that happening                                  ## 
  ##  Equation 9 in Cochran and Ellner 1992                                 ##
  ##  Argument/s:                                                           ##
  ##  P: Survival matrix P                                                  ##
  ##  Return Value/s:                                                       ##
  ##  The average time to reach stage i from stage j                        ##
  ##  Author/s:                                                             ##
  ##  Dr. Stephen Ellner                                                    ##
  ##  Dr. Simone Blomberg                                                   ##
  ##  Erin Souder                                                           ##
  ##  Date: 03/01/2023                                                      ##
  ############################################################################
  
  if(!is.matrix(P)) {stop("This is not a matrix")}
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  
  D <- Di_mat(P)
  results <- matrix(0, nrow= nrow(P), ncol= ncol(P))
  # Make empty matrix 
  
  for (i in 1:dim(results)[1]) {
    Imat <- diag(x=1, dim(D[,,i])[1])
    den <- solve(Imat - D[,,i])  
    num <- den %*% den
    results[i,] <- (num[i,]/den[i,]) - 1
  }
  # for loop for each Di_mat
  results[is.nan(results)] <- NA
  # NA used because can't get to this state from previous state
  
  return(results)
  
}

meantime_SD <- function(P){
    
  ######################################################################################
  ##  Function to determine the SD of the time to reach stage i from stage j          ##
  ##     conditional on that happening.                                               ##  
  ##  Equation 10 in Cochran and Ellner 1992                                          ##
  ##  Argument/s:                                                                     ##
  ##  P: Survival matrix P                                                            ##
  ##  Return Value/s:                                                                 ##
  ##  The standard deviation for the time to reach stage i from stage j,              ##
  ##           conditional on that happening 
  ##  Author/s:                                                                       ##
  ##  Dr. Stephen Ellner                                                              ##
  ##  Dr. Simone Blomberg                                                             ##
  ##  Erin Souder                                                                     ##
  ##  Date: 03/01/2023                                                                ##
  ######################################################################################
  
  if(!is.matrix(P)) {stop("This is not a matrix")} 
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  
  #OUTPUT DOESN'T MATCH STAGECOACH
  D <- Di_mat(P)
  m0 <- meantime(P)
  m1 <- m0 + 1
  # m1 used to match equations from paper
  results <- matrix(0, nrow = nrow(P), ncol= ncol(P))
  # Make new matrix 
  
  for (i in 1:dim(results)[1]) {
    Imat <- diag(x=1, dim(D[,,i])[1])
    num <- (Imat + D[,,i]) %*% solve((Imat - D[,,i]) %*% (Imat - D[,,i]) %*% (Imat - D[,,i])) 
    den <- solve(Imat - D[,,i]) 
    results[i,] <- (num[i,]/den[i,])
    results[i,] <- sqrt(results[i,] - (m1[i,]) ^ 2)
  }
  results[is.nan(results)] <- NA
  results[m0 < 0] <- NA
  # Transition not defined for this so it's NA
 return(results)
  
}

total_lifeSpan <- function(P){
  
  
  #########################################################################
  ##  Function to determine total conditional lifespan for each stage i  ##
  ##  Equation 6                                                         ##
  ##  Argument/s:                                                        ##
  ##  P: Survival matrix P                                               ##
  ##  Return Value/s:                                                    ##
  ##  The total lifespan for each stage i                                ##
  ##  Author/s:                                                          ##
  ##  Dr. Stephen Ellner                                                 ##
  ##  Dr. Simone Blomberg                                                ##
  ##  Erin Souder                                                        ##
  ##  Date: 03/01/2023                                                   ##
  #########################################################################
  
  if(!is.matrix(P)) {stop("This is not a matrix")}
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  
  LE <- life_expectancy(P)
 # must add each stage number to the correct vector
 MT <- meantime(P)
 results <- matrix(0,nrow = nrow(P),ncol = ncol(P))
 for (i in 1:dim(P)[1]) {
   for (j in 1:dim(P)[2]) {
     results[i,j] <- MT[i,j] + LE[i] + 1
     if(i < j){
       results[i,j] <- -1
     }
     else{
       results[i,j] <- MT[i,j] + LE[i] + 1
   }
     }
 }
 results[results < 0] < NA
 # NA used because can't get to state from previous state
  return(results)
}


total_lifeSpan_SD <- function(P){
  
  
  ###################################################################################
  ##  Function to determine the SD of total conditional lifespan for each stage i  ##
  ##  Equation 7                                                                   ##
  ##  Argument/s:                                                                  ##
  ##  P: Survival matrix P                                                         ##
  ##  Return Value/s:                                                              ##
  ##  The standard deviation of total lifespan for each stage i                    ##
  ##  Author/s:                                                                    ##
  ##  Dr. Stephen Ellner                                                           ##
  ##  Dr. Simone Blomberg                                                          ##
  ##  Erin Souder                                                                  ##
  ##  Date: 03/01/2023                                                             ##
  ###################################################################################
  
  if(!is.matrix(P)) {stop("This is not a matrix")}
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  
  m <- meantime_SD(P)
  l <- life_expectancy_SD(P)
  results <- matrix(0,nrow = nrow(P), ncol = ncol(P))
  for (i in 1:nrow(P)) {
    
      results[i,] <- l[i] ^ 2 + m[i,] ^ 2  
    }
    
  return(sqrt(results))
  
}


lx <- function (P, newbornTypes = NULL, max = 20) {
  
  ##########################################################################################
  ##  Function to determine the probability of survival to age x for type j-newborn (lx)  ##
  ##  Equation 2                                                                          ##
  ##  Argument/s:                                                                         ##
  ##  P: Survival matrix P                                                                ##
  ##  newbornTypes: type j newborn default NULL                                           ##
  ##  max: max number of repitions through equation default 20                            ##
  ##  Return Value/s:                                                                     ##
  ##  The probability of survival to age x for j type newborn                             ##
  ##  Author/s:                                                                           ##
  ##  Dr. Stephen Ellner                                                                  ##
  ##  Dr. Simone Blomberg                                                                 ##
  ##  Erin Souder                                                                         ##
  ##  Date: 03/01/2023                                                                    ##
  ##########################################################################################
  
  if(!is.matrix(P)) {stop("This is not a matrix")}
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  if(is.null(newbornTypes)) newbornTypes = c(1:ncol(P));
  
  results <- matrix(1,max,ncol(P))
  
  for (x in 2:(max)) {
   results[x,] = colSums(P %^% (x - 1)) 
    
  }
  # for loop 
  return(results[,newbornTypes])
  
}


lx_pop <- function (A,B, P, max = 20) {

  ##################################################################################
  ##  Function to determine the population survival probability (lx pop)          ##
  ##  Table 2                                                                     ##
  ##  Argument/s:                                                                 ##
  ##  A: Transition matrix                                                        ##
  ##  P: Survival matrix P                                                        ##
  ##  B: Birth matrix                                                             ##
  ##  max: max number of times to run through equation default 20                 ##
  ##  Return Value/s:                                                             ##
  ##  The probability of survival to age x for j type newborn for the population  ##
  ##  Author/s:                                                                   ##
  ##  Dr. Stephen Ellner                                                          ##
  ##  Dr. Simone Blomberg                                                         ##
  ##  Erin Souder                                                                 ##
  ##  Date: 03/01/2023                                                            ##
  ##################################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(P)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  
  l <- lx(P,max = max)
  b <- n_bj(A,B)
  results = l %*% b;
 
  return(results)
}


fx_weighted <- function(A,B,C, newbornTypes = NULL, max= 20 ){
  
  
  #################################################################################
  ##  Function to determine weighted reproduction at age x (fx weighted)          ##
  ##  Equation 13                                                                ##
  ##  Argument/s:                                                                ##
  ##  A: Transition matrix                                                       ##
  ##  C: Survival matrix P + Fission matrix F                                    ##
  ##  B: Birth matrix                                                            ##
  ##  newbornType: type j newborn default NULL                                   ##
  ##  max: max number of times to run through equation default 20                ##
  ##  Return Value/s:                                                            ##
  ##  The weighted average number of offspring produced by individuals at age x  ##
  ##  Author/s:                                                                  ##
  ##  Dr. Stephen Ellner                                                         ##
  ##  Dr. Simone Blomberg                                                        ##
  ##  Erin Souder                                                                ##
  ##  Date: 03/01/2023                                                           ##
  #################################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  if(is.null(newbornTypes)) newbornTypes = c(1:ncol(P));
  
  results <- matrix(NA, max, ncol(P));
  gamma_i = gam_i(A,B);
  
  for (x in 1:max) {
    Cx1 = C %^% (x - 1)
    num = t(Cx1) %*% gamma_i
    # Transpose used here
    den = colSums(Cx1);
    results[x,] = num/den;
    results[x,den == 0] = 0
  }
  #for loop in order to calculate the total maternity (f(x))
  return(results[,newbornTypes])
}

fx_unweighted <- function(A,B,C,newbornTypes = NULL, max = 20){
  
  
  ###################################################################################
  ##  Function to determine unweighted reproduction at age x (fx unweighted)       ##
  ##  Equation 13                                                                  ##
  ##  Argument/s:                                                                  ##
  ##  A: Transition matrix                                                         ##
  ##  C: Survival matrix P + Fission matrix F                                      ##
  ##  B: Birth matrix                                                              ##
  ##  newbornType: type j newborn default NULL                                     ##
  ##  max: max number of times to run through equation default 20                  ##
  ##  Return Value/s:                                                              ##
  ##  The unweighted average number of offspring produced by individuals at age x  ##
  ##  Author/s:                                                                    ##
  ##  Dr. Stephen Ellner                                                           ##
  ##  Dr. Simone Blomberg                                                          ##
  ##  Erin Souder                                                                  ##
  ##  Date: 03/01/2023                                                             ##
  ###################################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  if(is.null(newbornTypes)) newbornTypes = c(1:ncol(P));
  # newbornTypes if only want stages that can reproduce
  
  results <- matrix(NA, max, ncol(P));
  beta_i <- bet_i(B);
  
  for(x in 1:max){
    Cx1 <- C %^% (x - 1)
    num <- t(Cx1) %*% beta_i
    den <- colSums(Cx1);
    results[x,] <- num/den;
    results[x, den == 0] = 0
  }
  return(results[,newbornTypes])
}

fx_pop_weighted <- function(A,B,C, max = 20){
  
  
  #################################################################################
  ##  Function to determine population reproduction at time t (fx weighted pop)  ##
  ##  Table 2                                                                    ##
  ##  Argument/s:                                                                ##
  ##  A: Transition matrix                                                       ##
  ##  C: Survival matrix P + Fission matrix F                                    ##
  ##  B: Birth matrix                                                            ##
  ##  max: max number of times to run through equation default 20                ##
  ##  Return Value/s:                                                            ##
  ##  The weighted population reproduction at time t                             ##
  ##  Author/s:                                                                  ##
  ##  Dr. Stephen Ellner                                                         ##
  ##  Dr. Simone Blomberg                                                        ##
  ##  Erin Souder                                                                ##
  ##  Date: 03/01/2023                                                           ##
  #################################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
   f <- fx_weighted(A,B,C)
  b <- n_bj(A,B)
 results <- f %*% b
 return(results)
}

fx_pop_unweighted <- function(A,B,C, max = 20){
  
  
  ##############################################################################################
  ##  Function to determine unweighted population reproduction at time t (fx unweighted pop)  ##
  ##  Table 2                                                                                 ##
  ##  Argument/s:                                                                             ##
  ##  A: Transition matrix                                                                    ##
  ##  C: Survival matrix P + Fission matrix F                                                 ##
  ##  B: Birth matrix                                                                         ##
  ##  max: max number of times to run through equation default 20                             ##
  ##  Return Value/s:                                                                         ##
  ##  The unweighted population reproduction at time t                                        ##
  ##  Author/s:                                                                               ##
  ##  Dr. Stephen Ellner                                                                      ##
  ##  Dr. Simone Blomberg                                                                     ##
  ##  Erin Souder                                                                             ##
  ##  Date: 03/01/2023                                                                        ##
  ##############################################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
  f <-fx_unweighted(A,B,C)
  b <- n_bj(A,B)
  results <- f %*% b
  return(results)
}


Vx_V1 <- function(A,B,C,newbornTypes = NULL, max = 20){
  
  
  #########################################################################
  ##  Function to determine the age specific reproductive value (Vx/V1)  ##
  ##  Equation 32                                                        ##
  ##  Argument/s:                                                        ##
  ##  A: Transition matrix                                               ##
  ##  C: Survival matrix P + Fission matrix F                            ##
  ##  B: Birth matrix                                                    ##
  ##  newbornType: type j newborn default NULL                           ##
  ##  max: max number of times to run through equation default 20        ##
  ##  Return Value/s:                                                    ##
  ##  The  age specific reproductive value for each stage                ##
  ##  Author/s:                                                          ##
  ##  Dr. Simone Blomberg                                                ##
  ##  Erin Souder                                                        ##
  ##  Date: 02/03/2023                                                   ##
  #########################################################################
  
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  if(is.null(newbornTypes)) newbornTypes = c(1:ncol(C));
  
  age <- function(A,B){
    # scaling the reproductive value so sum v * bj = 1
    #this becomes first value in the loop for Vx_V1
    num <-  repro_value(A)
    den <- sum(repro_value(A) * n_bj(A,B))
    results <- (num/den) 
    return(results)
  }
  
  
  v <- age(A,B)
  results <- matrix(v,max,ncol(C))
  
  for (x in 2:max) {
    Cx1 <- C %^% (x - 1)
    num <- t(Cx1) %*% v
    # Use transpose here so multiplication takes place
    den <- colSums(Cx1)
    results[x,] <- num/den
    results[x, den == 0] = 0
    
  }
  return(results)
}


Vx_V1_pop <- function(A,B,C, max = 20){
  
  
  ############################################################################################
  ##  Function to determine the age specific reproductive value for population (Vx/V1 pop)  ##
  ##  Equation 33                                                                           ##
  ##  Argument/s:                                                                           ##
  ##  A: Transition matrix                                                                  ##
  ##  C: Survival matrix P + Fission matrix F                                               ##
  ##  B: Birth matrix                                                                       ##
  ##  newbornType: type j newborn default NULL                                              ##
  ##  max: max number of times to run through equation default 20                           ##
  ##  Return Value/s:                                                                       ##
  ##  The average age specific reproductive value for the population                        ##
  ##  Author/s:                                                                             ##
  ##  Dr. Simone Blomberg                                                                   ##
  ##  Erin Souder                                                                           ##
  ##  Date: 14/03/2023                                                                      ##
  ############################################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
 
  vx <- Vx_V1(A,B,C,max)
  b <- n_bj(A,B)
  
  results  = vx %*% b
 return(results)
}



net_rep <- function(B,C,newbornTypes = NULL){
  
  
  #######################################################################
  ##  Function to calculate the net reproductive rate (R0)              ##
  ##  Equation 17                                                      ##
  ##  Argument/s:                                                      ##
  ##  C: Survival matrix P + Fission matrix F                          ##
  ##  B: Birth matrix                                                  ##
  ##  newbornType: type j newborn                                      ##
  ##  Return Value/s:                                                  ##
  ##  The average number of offspring during an individual's lifetime  ##
  ##  Author/s:                                                        ##
  ##  Erin Souder                                                      ##
  ##  Date: 14/03/2023                                                 ##
  #######################################################################
  
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  if(is.null(newbornTypes)) newbornTypes = c(1:ncol(C));
  
  
  Imat <- diag(dim(C)[1]) # Identity matrix
  beta_i <- bet_i(B)
  results <- t(solve(Imat - C)) %*% beta_i 
  # need to use transpose here
  results[newbornTypes,]
}


net_rep_pop <- function(A,B,C){
  
  
  ##########################################################################################
  ##  Function to calculate the net reproductive rate for the population (R0 pop)          ##
  ##  Table 2                                                                             ##
  ##  Argument/s:                                                                         ##
  ##  A: Transition matrix A                                                              ##
  ##  C: Survival matrix P + Fission matrix F                                             ##
  ##  B: Birth matrix                                                                     ##
  ##  newbornType: type j newborn                                                         ##
  ##  Return Value/s:                                                                     ##
  ##  The average number of offspring during an individual's lifetime for the population  ##
  ##  Author/s:                                                                           ##
  ##  Erin Souder                                                                         ##
  ##  Date: 14/03/2023                                                                    ##
  ##########################################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
  R <- net_rep(B,C)
  b <- n_bj(A,B)
  results <- R %*% b
 results
}


average_age_production_unweighted <- function(A,B,C,newbornTypes = NULL){
  
  
  ###################################################################
  ##  Function to determine the average age of first reproduction  ##
  ##  Equation 27                                                  ##
  ##  Argument/s:                                                  ##
  ##  A: Transition matrix A                                       ##
  ##  C: Survival matrix P + Fission matrix F                      ##
  ##  B: Birth matrix                                              ##
  ##  newbornType: type j newborn                                  ##
  ##  Return Value/s:                                              ##
  ##  The average age of first reproduction for each stage j       ##
  ##  Author/s:                                                    ##
  ##  Dr. Simone Blomberg                                          ##
  ##  Erin Souder                                                  ##
  ##  Date: 28/03/2023                                             ##
  ###################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  if(is.null(newbornTypes)) newbornTypes = c(1:ncol(C));
  
  Imat <- diag(dim(C)[1]) # Identity matrix
  betai <- bet_i(B)
  num1 <- solve((Imat - C) %*% (Imat - C))
  num <- t(num1) %*% betai
  # Use beta for unweighted numbers
  den <- net_rep(B,C)
  results <- num/den
  results
}

average_age_production_weighted <- function(A,B,C,newbornTypes = NULL){
  
  
  ###################################################################
  ##  Function to determine the average age of first reproduction  ##
  ##  Equation 27                                                  ##
  ##  Argument/s:                                                  ##
  ##  A: Transition matrix A                                       ##
  ##  C: Survival matrix P + Fission matrix F                      ##
  ##  B: Birth matrix                                              ##
  ##  newbornType: type j newborn                                  ##
  ##  Return Value/s:                                              ##
  ##  The average age of first reproduction for each stage j       ##
  ##  Author/s:                                                    ##
  ##  Dr. Simone Blomberg                                          ##
  ##  Erin Souder                                                  ##
  ##  Date: 28/03/2023                                             ##
  ###################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  if(is.null(newbornTypes)) newbornTypes = c(1:ncol(C));
  
  Imat <- diag(dim(C)[1]) # Identity matrix
  gammai <- gam_i(A,B)
  num1 <- solve((Imat - C) %*% (Imat - C))
  num <- t(num1) %*% gammai
  # Use gam for weighted numbers
  den <- net_rep(B,C)
  results <- num/den
  results
}

average_age_production_SD_unweighted <- function(A,B,C,newbornTypes = NULL){
  
  
  #############################################################################################
  ##  Function to determine the standard deviation of the average age at first reproduction  ##
  ##  Equation 28                                                                            ##
  ##  Argument/s:                                                                            ##
  ##  A: Transition matrix A                                                                 ##
  ##  C: Survival matrix P + Fission matrix F                                                ##
  ##  B: Birth matrix                                                                        ##
  ##  newbornType: type j newborn                                                            ##
  ##  Return Value/s:                                                                        ##
  ##  The standard deviation of the average age of first reproduction for each stage j       ##
  ##  Author/s:                                                                              ##
  ##  Dr. Simone Blomberg                                                                    ##
  ##  Erin Souder                                                                            ##
  ##  Date: 28/03/2023                                                                       ##
  #############################################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  if(is.null(newbornTypes)) newbornTypes = c(1:ncol(C));
  
  Imat <- diag(dim(C)[1]) # Identity matrix
  betai <- bet_i(B)
  # Use beta for unweighted numbers
  a <- average_age_production(A,B,C)
  n1 <- (I + C) %*% solve((I - C) %*% (I - C) %*% (I - C))
  num <- t(n1) %*% betai
  den <- net_rep(B,C)
  results <- sqrt((num/den) - (a ^ 2))
  results
}

  #############################################################################################
  ##  Function to determine the standard deviation of age at first reproduction              ##
  ##  Equation 28                                                                            ##
  ##  Argument/s:                                                                            ##
  ##  A: Transition matrix A                                                                 ##
  ##  C: Survival matrix P + Fission matrix F                                                ##
  ##  B: Birth matrix                                                                        ##
  ##  newbornType: type j newborn                                                            ##
  ##  Return Value/s:                                                                        ##
  ##  The standard deviation of the average age of first reproduction for each stage j       ##
  ##  Author/s:                                                                              ##
  ##  Dr. Simone Blomberg                                                                    ##
  ##  Erin Souder                                                                            ##
  ##  Date: 28/03/2023                                                                       ##
  #############################################################################################
average_age_production_SD_weighted <- function(A,B,C,newbornTypes = NULL){
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  if(is.null(newbornTypes)) newbornTypes = c(1:ncol(C));
  
  Imat <- diag(dim(C)[1]) # Identity matrix
  gammai <- gam_i(A,B)
  # Use gamma for weighted numbers
  a <- average_age_production(A,B,C)
  n1 <- (I + C) %*% solve((I - C) %*% (I - C) %*% (I - C))
  num <- t(n1) %*% gammai
  den <- net_rep(B,C)
  results <- sqrt((num/den) - (a ^ 2))
  results
}

average_age_production_pop_unweighted <- function(A,B,C){
  
  
  ###################################################################
  ##  Function to determine the average age of first reproduction  ##
  ##  Table 2                                                      ##
  ##  Argument/s:                                                  ##
  ##  A: Transition matrix A                                       ##
  ##  C: Survival matrix P + Fission matrix F                      ##
  ##  B: Birth matrix                                              ##
  ##  newbornType: type j newborn                                  ##
  ##  Return Value/s:                                              ##
  ##  The average age of first reproduction for each stage j       ##
  ##  Author/s:                                                    ##
  ##  Dr. Simone Blomberg                                          ##
  ##  Erin Souder                                                  ##
  ##  Date: 28/03/2023                                             ##
  ###################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
  
  ##################
  # I am not sure if this is the correct number so if someone else could check the math
  ##################

  Imat <- diag(dim(C)[1]) # Identity matrix
  betai<- bet_i(B)
  bj <- n_bj(A,B)
  n1 <- solve((Imat - C) %*% (Imat - C))
  num <- t(n1) %*% betai 
  # Use beta for unweighted numbers
  den <- net_rep(B,C)
  results <- t(num/den) %*% bj
  results
}

average_age_production_pop_weighted <- function(A,B,C){
  
  
  ###################################################################
  ##  Function to determine the average age of first reproduction  ##
  ##  Table 2                                                      ##
  ##  Argument/s:                                                  ##
  ##  A: Transition matrix A                                       ##
  ##  C: Survival matrix P + Fission matrix F                      ##
  ##  B: Birth matrix                                              ##
  ##  newbornType: type j newborn                                  ##
  ##  Return Value/s:                                              ##
  ##  The average age of first reproduction for each stage j       ##
  ##  Author/s:                                                    ##
  ##  Dr. Simone Blomberg                                          ##
  ##  Erin Souder                                                  ##
  ##  Date: 28/03/2023                                             ##
  ###################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
  
  ##################
  # I am not sure if this is the correct number so if someone else could check the math
  ##################
  
  Imat <- diag(dim(C)[1]) # Identity matrix
  gammai <- gam_i(A,B)
  bj <- n_bj(A,B)
  n1 <- solve((Imat - C) %*% (Imat - C))
  num <- t(n1) %*% gammai 
  # Use gamma for weighted numbers
  den <- net_rep(B,C)
  results <- t(num/den) %*% bj
  results
}

average_age_production_SD_pop_unweighted <- function(A,B,C){
  
  
  #############################################################################################
  ##  Function to determine the standard deviation of the average age at first reproduction  ##
  ##  Equation 28                                                                            ##
  ##  Argument/s:                                                                            ##
  ##  A: Transition matrix A                                                                 ##
  ##  C: Survival matrix P + Fission matrix F                                                ##
  ##  B: Birth matrix                                                                        ##
  ##  newbornType: type j newborn                                                            ##
  ##  Return Value/s:                                                                        ##
  ##  The standard deviation of the average age of first reproduction for each stage j       ##
  ##  Author/s:                                                                              ##
  ##  Dr. Simone Blomberg                                                                    ##
  ##  Erin Souder                                                                            ##
  ##  Date: 28/03/2023                                                                       ##
  #############################################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
  ########### 
  # Can someone check math? Not sure if this is right #
  ###########
  
  Imat <- diag(dim(C)[1]) # Identity matrix
  betai <- bet_i(B)
  # Use beta for weighted numbers
  bj <- n_bj(A,B)
  a <- average_age_production_pop_unweighted(A,B,C)
  n1 <- (I + C) %*% solve((I - C) %*% (I - C) %*% (I - C))
  num <- (t(n1) %*% betai) 
  den <- net_rep(B,C)
  r1 <- t(num/den) %*% bj 
  results <- sqrt((r1) - (a ^ 2))
  results 
}

average_age_production_SD_pop_weighted <- function(A,B,C){
  
  
  #############################################################################################
  ##  Function to determine the standard deviation of the average age at first reproduction  ##
  ##  Equation 28                                                                            ##
  ##  Argument/s:                                                                            ##
  ##  A: Transition matrix A                                                                 ##
  ##  C: Survival matrix P + Fission matrix F                                                ##
  ##  B: Birth matrix                                                                        ##
  ##  newbornType: type j newborn                                                            ##
  ##  Return Value/s:                                                                        ##
  ##  The standard deviation of the average age of first reproduction for each stage j       ##
  ##  Author/s:                                                                              ##
  ##  Dr. Simone Blomberg                                                                    ##
  ##  Erin Souder                                                                            ##
  ##  Date: 28/03/2023                                                                       ##
  #############################################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
  ##########
  # Can someone check the math to check if this is right? #
  ##########
  
  Imat <- diag(dim(C)[1]) # Identity matrix
  gammai <- gam_i(A,B)
  # Use gamma for weighted numbers
  bj <- n_bj(A,B)
  a <- average_age_production_pop_weighted(A,B,C)
  n1 <- (I + C) %*% solve((I - C) %*% (I - C) %*% (I - C))
  num <- (t(n1) %*% gammai)
  den <- net_rep(B,C)
  r1 <- t(num/den) %*% bj
  results <- sqrt(abs(r1 - (a ^2)))
  results 
}

mean_age_residence <- function(C){
  
  
  #################################################################
  ##  Function to determine the average age of stage j (Si)      ##
  ##  Equation 29                                                ##
  ##  Argument/s:                                                ##
  ##  C: Survival matrix P + Fission matrix F                    ##
  ##  newbornType: type j newborn                                ##
  ##  Return Value/s:                                            ##
  ##  The average age of stage j                                 ##
  ##  Author/s:                                                  ##
  ##  Dr. Stephen Ellner                                         ##
  ##  Dr. Simone Blomberg                                        ##
  ##  Erin Souder                                                ##
  ##  Date: 03/01/2023                                           ##
  #################################################################

  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
#Need to do by newbornType and make NaN = 0
Imat <- diag(dim(C)[1]) #Identity matrix
# Identity matrix
num <- solve((Imat - C) %^% 2)
den <- solve(Imat - C)
results <- (num/den)
results[is.nan(results)] <- 0
results
}

mean_age_residence_SD <- function(C){
  
  
  #######################################################################################
  ##  Function to determine the standard deviation of the average age of stage j (Si)  ##
  ##  Equation 30                                                                      ##
  ##  Argument/s:                                                                      ##
  ##  C: Survival matrix P + Fission matrix F                                          ##
  ##  newbornType: type j newborn                                                      ##
  ##  Return Value/s:                                                                  ##
  ##  The standard deviation of the average age of stage j                             ##
  ##  Author/s:                                                                        ##
  ##  Dr. Simone Blomberg                                                              ##
  ##  Erin Souder                                                                      ##
  ##  Date: 03/01/2023                                                                 ##
  #######################################################################################
  
 
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(C)) {stop("This is not numeric")}
   #Need to fix this so it just gives the newbornTypes 1,3,4,5
  
  S <- mean_age_residence(C)
  Imat <- diag(dim(C)[1])  # identity matrix
 
  num <- (Imat + C) %*% solve((Imat - C) %^% 3)
  den <- solve(Imat - C)
  results <- abs(sqrt((num / den) - S ^ 2))
  results[is.nan(results)] <- 0
 results
  
}

mean_age_residence_pop <- function(A,B,C){
  
  
  ############################################################################
  ##  Function to determine the population average age of stage j (Si pop)  ##
  ##  Table 2                                                               ##
  ##  Argument/s:                                                           ##
  ##  A: Transition matrix A                                                ##
  ##  B: Birth matrix B                                                     ##
  ##  C: Survival matrix P + Fission matrix F                               ##
  ##  Return Value/s:                                                       ##
  ##  The population average age of stage j                                 ##
  ##  Author/s:                                                             ##
  ##  Dr. Simone Blomberg                                                   ##
  ##  Erin Souder                                                           ##
  ##  Date: 03/01/2023                                                      ##
  ############################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
  bj <- n_bj(A,B)
  Imat <- diag(dim(C)[1]) # Identity matrix
  num <- solve((Imat - C) %^% 2) %*% bj
  den <- solve(Imat - C) %*% bj
  results <- num/den
  results
  
}


mean_age_residence_pop_SD <- function(A,B,C){
  

  ######################################################################################
  ##  Function to determine the population SD of the average age of stage j (Si pop)  ##
  ##  Table 2                                                                         ##
  ##  Argument/s:                                                                     ##
  ##  A: Transition matrix A                                                          ##
  ##  B: Birth matrix B                                                               ##
  ##  C: Survival matrix P + Fission matrix F                                         ##
  ##  Return Value/s:                                                                 ##
  ##  The population standard deviation of the average age of stage j                 ##
  ##  Author/s:                                                                       ##
  ##  Dr. Simone Blomberg                                                             ##
  ##  Erin Souder                                                                     ##
  ##  Date: 03/01/2023                                                                ##
  ######################################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  
  #DOESN'T MATCH OUTPUT
  s <- mean_age_residence_pop(A,B,C)
  bj <- n_bj(A,B)
  Imat <- diag(x = 1, dim(C)) # Identity matrix
  # identity matrix
  num <- ((Imat + C) %*% solve(Imat - C) %^% (Imat - C) %*% (Imat - C)) %*% bj
  den <- (solve(Imat - C)) %*% bj
  results <- sqrt(abs((num/den) - s ^ 2))
  results
}


Q_mat <- function(P,B){
  
  
  ###################################################################
  ##  Function to calculate new transition matrix for reproduction  ##
  ##  Equation 14                                                  ##
  ##  Argument/s:                                                  ##
  ##  P: Survival matrix P                                         ##
  ##  stage: stage at which reproduction occurs (j)                ##
  ##  Return Value/s:                                              ##
  ##  Transition martrix                                           ##
  ##  Author/s:                                                    ##
  ##  Dr. Simone Blomberg                                          ##
  ##  Erin Souder                                                  ##
  ##  Date: 03/01/2023                                             ##
  ###################################################################
  
  if(!is.matrix(P)) {stop("This is not a matrix")}
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  
  results <- array(0, nrow = nrow(P), ncol = ncol(P))
  for (i in 1:dim(P)[1]) {
    for (j in  1: dim(P)[2]){
      if(j == B[2]){
        results[i,j] <- 0
      }
      else{
        results[i,j] <- P[i,j]
      }
    }
    }
  
 results
}  


maturity_age <- function(P, Q_mat, B, stage, newbornTypes = NULL ){
  
  
  #################################################################
  ##  Function to determine average age of maturity (E)          ##
  ##  Equation 15                                                ##
  ##  Argument/s:                                                ##
  ##  P: Survival matrix P                                       ##
  ##  Q_mat: Equation 14                                         ##
  ##  stage: stage at which reproduction occurs [i,]             ##
  ##  newbornTypes: newborn type j                               ##
  ##  Return Value/s:                                            ##
  ##  Average age of maturity at each stage j                    ##
  ##  Author/s:                                                  ##
  ##  Dr. Simone Blomberg                                        ##
  ##  Erin Souder                                                ##
  ##  Date: 03/01/2023                                           ##
  #################################################################
  
  if(!is.matrix(P)) {stop("This is not a matrix")}
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  if(is.null(newbornTypes)) newbornTypes = c(1:ncol(P));
  
  ### HOW DO I PULL JUST [6,] ####
  
  Q <- Q_mat(P,B)
  Imat <- diag(x = 1, dim(P)) # Identity matrix
  # identity matrix
  num <- (solve((Imat - Q) %*% (Imat - Q)))
  den <- (solve(Imat - Q))
  
  results <- num/den  
  
  results[stage,]
}


maturity_age_SD <- function(P, Q_mat, B, stage,newbornTypes = NULL){
  
  
  ###################################################################################
  ##  Function to determine standard deviation of the average age of maturity (E)  ##
  ##  Equation 16                                                                  ##
  ##  Argument/s:                                                                  ##
  ##  P: Survival matrix P                                                         ##
  ##  Q_mat: Equation 14                                                           ##
  ##  stage: stage at which reproduction occurs                                    ##
  ##  newbornType: newborn type j                                                  ##
  ##  Return Value/s:                                                              ##
  ##  Standard deviation of the average age of maturity at each stage j            ##
  ##  Author/s:                                                                    ##
  ##  Dr. Simone Blomberg                                                          ##
  ##  Erin Souder                                                                  ##
  ##  Date: 03/01/2023                                                             ##
  ###################################################################################

  if(!is.matrix(P)) {stop("This is not a matrix")}
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  if(is.null(newbornTypes)) newbornTypes = c(1:ncol(P));
  
  ###### HOW DO I PULL JUST [6,] ####
  
 Q <- Q_mat(P,B)
 m <- maturity_age(P,Q_mat,B,stage)
 Imat <- diag(dim(P)[1]) # Identity matrix
 num <- (Imat + Q) %*% (solve((Imat - Q) %*% (Imat - Q) %*% (Imat - Q)))
 den <- solve(Imat - Q)
 results <- sqrt(abs(num/den) - m^2)
 results[stage,]
}

maturity_age_pop <- function(A, P, Q_mat, B){
  
  
  ########################################################################
  ##  Function to determine average age of maturity (E) for population  ##
  ##  Table 2                                                           ##
  ##  Argument/s:                                                       ##
  ##  P: Survival matrix P                                              ##
  ##  A: Transition matrix A                                            ##
  ##  B: Birth matrix B                                                 ##
  ##  Q_mat: Equation 14                                                ##
  ##  stage: stage at which reproduction occurs                         ##
  ##  newbornType: newborn type j                                       ##
  ##  Return Value/s:                                                   ##
  ##  Average age of maturity at each stage j                           ##
  ##  Author/s:                                                         ##
  ##  Dr. Simone Blomberg                                               ##
  ##  Erin Souder                                                       ##
  ##  Date: 17/04/2023                                                  ##
  ########################################################################
  
  
  if(!is.matrix(P)) {stop("This is not a matrix")}
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  
  #### RESULTS DO NOT MATCH PAPER ####
  
  Q <- Q_mat(P,B)
  Imat <- diag(x = 1, dim(P)) # Identity matrix
  # identity matrix
  bj <- n_bj(A,B)
  num <- (solve((Imat - Q) %*% (Imat - Q))) %*% bj
  den <- (solve(Imat - Q)) %*% bj
  
  results <- (num/den)
  
  results
}


generation_time <- function(A,B,C,newbornTypes = NULL){
  
  
  #################################################################
  ##  Function to determine generation time                      ##
  ##  Argument/s:                                                ##
  ##  A: Transition matrix A                                     ##
  ##  B: Birth matrix B                                          ##
  ##  C: Survival matrix P + Fission matrix F                    ##
  ##  newbornType: newborn type j                                ##
  ##  Return Value/s:                                            ##
  ##  Generation time for each stage j                           ##
  ##  Author/s:                                                  ##
  ##  Erin Souder                                                ##
  ##  Date: 03/01/2023                                           ##
  #################################################################
  
  if(!is.matrix(A)) {stop("This is not a matrix")}
  if(!is.matrix(B)) {stop("This is not a matrix")}
  if(!is.matrix(C)) {stop("This is not a matrix")}
  if(!is.square.matrix(A)) {stop("This is not a square matrix")}
  if(!is.square.matrix(B)) {stop("This is not a square matrix")}
  if(!is.square.matrix(C)) {stop("This is not a square matrix")}
  if(!is.numeric(A)) {stop("This is not numeric")}
  if(!is.numeric(B)) {stop("This is not numeric")}
  if(!is.numeric(C)) {stop("This is not numeric")}
  if(is.null(newbornTypes)) newbornTypes = c(1:ncol(C));
  
  lam <- pop_growth(A)
  R <- sum(net_rep_pop(A,B,C))
  results <- log(R) / log(lam)
  results
}

#################################################################
##                       RUNNING PROGRAM                       ##
#################################################################

pop_growth(A)

stable_stage_dist(A)

reproductive_value(A)

sensitivy_mat(A)

elasticity_mat(A)

n_bj(A,B)

age_in_stage(A,B,C)

age_in_stage_SD(A,B,C)

gam_i(A,B)

bet_i(B)

pop_gen_time(A,B,C)

pit(A,B,C, maxAge = 10)

stable_age_distribution(A,B,C, maxAge = 10)

life_expectancy(P)

life_expectancy_SD(P)

Di_mat(P)

meantime(P)

meantime_SD(P)

total_lifeSpan(P)

total_lifeSpan_SD(P)

lx(P,newbornTypes = NULL, max = 20)

lx_pop(A,B,P, max = 20)

fx_weighted(A,B,C,newbornTypes = NULL, max = 20)

fx_unweighted(A,B,C,newbornTypes = NULL, max = 20)

fx_pop_weighted(A,B,C, max = 20)

fx_pop_unweighted(A,B,C, max = 20)

Vx_V1(A,B,C,newbornType = NULL, max = 20)

Vx_V1_pop(A,B,C, max = 20)

net_rep(B,C, newbornType = c(1,3,4,5))

net_rep_pop(A,B,C, newbornType = c(1,3,4,5))

average_age_production(A,B,C,c(1,3,4,5))

average_age_production_SD(A,B,C,c(1,3,4,5))

mean_age_residence(C)

mean_age_residence_SD(C)

mean_age_residence_pop(A,B,C)

mean_age_residence_pop_SD(A,B,C)

Q_mat(P, c(6))

maturity_age(P,Q_mat,stage = c(6),newbornType = c(1,3,4,5))

maturity_age_SD(P,Q_mat,stage = c(6),newbornType = c(1,3,4,5))

generation_time(A,B,C,newbornType = c(1,3,4,5))