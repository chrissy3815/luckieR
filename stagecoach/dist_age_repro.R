###############################################################################
## Function for mean age at first reproduction, and the distribution of stage
## at first reproduction, both conditional on reproducion at least once. 
##  
## Formulas are from IPMbook, chapter 3, esp. page 71 for distribution. 
## Those formulas assume the model includes an explicit breeding probability
## p_b(z). Here we just input a COM(P)ADRE-style Umat and Fmat, and infer
## a 'breeding probability' that is actually the probability of having at least 
## one recruit that survives to be censused. Options for the inferred breeding
## probability are 
##   	Poisson, p_b(z) = 1 - dpois(0,lambda=meanFecundity(z)), 
##		Bernoulli, p_b(z) = meanFecundity(z)
## 		Fixed, p_b(z) = 1 if meanFecundity(z)>0, zero otherwise.  
##
##  NOTE, metaLuck_fns.R uses a different definition of Fixed. I think
##     this one is right.  
##
## We do not include costs of reproduction, so in IPMbook notation, 
## pi_0 = pi_b = U in eqn. (3.2.9)
##
## NOTE CAREFULLY: 
## If the model is a Leslie matrix in the strict sense (finite maximum age)
## then the distribution of stages at first reproduction implies the distribution
## of ages at first reproduction. 
##
## For other models, the function it gives the mean AGE at first reproduction
## correctly, but the stage distribution is the distribution of STAGES 
## at first reproduction, which may or may not be biologically meaningful.  
##
## NOTE ALSO that mixdist does not make sense for an age-structured model
## because everyone is the same age at birth.  
################################################################################

dist_age_repro<- function(Umat, Fmat, repro_var = 'poisson', mixdist=NULL){
  ## take the column sum of reproduction:
  betabar<- colSums(Fmat) #betabar is the average offspring production
  
  ## Breeding probability values, p_b:
  if (repro_var %in% c("poisson", "Poisson")){
    p_b<- 1-dpois(0, lambda=betabar) # poisson: prob(x>0)
  } else if (repro_var %in% c("Bernoulli", "bernoulli")){
    p_b<- betabar # bernoulli: prob(occurrence)
  } else if (repro_var %in% c("fixed", "Fixed")){
    p_b<- ifelse(betabar>0,1,0) # all individuals reproduce, if any do 
  }
  
  ## Modified kernel in which breeding is a second kind of death,
  ## IPMbook p.70 after eqn. 3.2.10. 
  nz = nrow(Umat); 
  P_0 = matrix(NA,nz,nz); 
  for(z in 1:nz) P_0[,z] = (1-p_b[z])*Umat[,z] 
  N_0<- solve(diag(nz) - P_0) # matrix_inverse(I-P_0)
  
  ## Probability of reproducing at least once (i.e, dies by breeding) 
  ## IPMbook p.70 
  B<- p_b%*%N_0
  
  ## Survival kernel P_b = P_0 conditioned on breeding at least once
  ## IPMbook eqn. 3.2.11 
  P_b = matrix(NA,nz,nz);  
  for(z in 1:nz) P_b[,z] = P_0[,z]*B/B[z]
  
  ## Calculate the expected age at first reproduction, conditional on
  ## state at birth. Note, we use the convention that age at birth=0
  ## so we have to subtract 1 from the usual colSums calculation.   
  e<- rep(1,nz)
  N_b<- solve(diag(nz) - P_b) # matrix_inverse(I-P_b)
  expArCond_z<- colSums(N_b)-1; 
  
  mean_age_repro = expArCond_z
  if(is.null(mixdist)){
    mean_age_repro_cohort = NA
  } else {
    ## Page 71 of IPM book, note that mean age of repro for a cohort 
    ## has to be weighted by the different probabilities of reproducing.
    mixdist<- t(mixdist)*B/sum(t(mixdist)*B)
    expAr<- sum(expArCond_z*mixdist)
    mean_age_repro_cohort = expAr
  }
  
  ## State-at-death distribution for P_b, thus the state at first 
  ## reproduction distribution, conditional on reproducting at least
  ## once, for the original model. IPMbook eqn. 3.2.8  
  s_b = colSums(P_b); 
  Omega = matrix(NA,nz,nz) 
  for(z in 1:nz) Omega[,z] = (1-s_b)*N_b[,z]; 
  
  if(is.null(mixdist)){
    dist_cohort=NA
    } else {
    dist_cohort=Omega%*%matrix(mixdist,ncol=1); 
   }	
  
  return(list(B=B, mean_age = mean_age_repro, mean_age_cohort=mean_age_repro_cohort,
				dist_stage = Omega, dist_cohort=dist_cohort))

}  

########### TESTING STUFF 
  
Umat = matrix( c(0,   0,   0,
                 0.5, 0,   0,  
				 0,   0.5, 0), 3,3, byrow=TRUE) 
				 
Fmat = matrix(0,3,3); 
Fmat[1,] = c(0, 1, 2); 
dist_age_repro(Umat,Fmat,repro_var='fixed'); 
## No breeding at age 0, breeding is certain at ages 1,2 
## All results check out OK. 

Fmat[1,] = c(0, 1, 5000); 
dist_age_repro(Umat,Fmat,repro_var='poisson',mixdist=c(1,0,0)); 
## Here breeding at age 2 is certain, breeding at age 1
## has probability 1-dpois(0,lambda=1) = 0.632, and there 
## is no chance of breeding at age 0. 
## 
## First, what is B? 
## B[3]=1, correct 
## B[2] = (breed in stage 2) + (don't breed in 2)(live to 3) = 
##         0.632 + (1-0.632)*0.5 = 0.816, correct. 
## B[1] = 0.5*B[2], correct 
##
## Conditional on breeding, what is the mean age? 
## B[3] = 0, you breed in your birth year, correct. 
## B[1] = 1 + B[2], as it should be. 
## If you are born into stage 2, the outcomes in which
## you breed are: 
##  - breed in stage 2 then die or don't: prob = 0.632, age = 0 
##  - don't breed in state 2, live to stage 3: prob = (1-0.632)*0.5 = 0.184
## The total probability is 0.816, so conditional probs are unconditional/0.816
## Mean age is then (0.632/0.816)*0 + (0.184/0.816)*1 = 0.225, correct. 
## 
## The returned conditional stage distribution is 
##  [1,] 0.0000000 0.0000000    0
##  [2,] 0.7746003 0.7746003    0
##  [3,] 0.2253997 0.2253997    1
## which matches the calculations for mean age. 
 
  