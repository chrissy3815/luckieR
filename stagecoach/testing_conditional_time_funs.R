rm(list=ls(all=TRUE)); 
u = Sys.info(); 
if(u["user"]=="Ellner") setwd("c:/repos/luckieR/stagecoach"); 

source("conditional_time_funs.R"); 

#######################################################
#### For testing conditional times 
#######################################################
M1 = matrix(c(0, 0,    0,  0,
              1, 0,    0,  0,
			  0, 0.5,  0,  0,
			  0, 0.5,  1,  0), 4,4,byrow=TRUE) 
N1 = solve(diag(4)-M1); colSums(N1); 			  
## 3.5 2.5 2.0 1.0  

M2 = matrix(c(0, 0,    0,  0,
              1, 0,    0,  0,
			  0, 0.25, 0,  0,
			  0, 0.25, 1,  0), 4,4,byrow=TRUE) 
N2 = solve(diag(4)-M2); colSums(N2); 			  
## 2.75 1.75 2.00 1.00


M7 = matrix(0,7,7); 
M7[2,1] = M7[6,3] = M7[7,5] = M7[7,6] = 1 ## states 1,3,5,6 
M7[3,2] = M7[4,2] = M7[5,2] = 1/3; # state 2 
M7[3,4] = M7[7,4] = 0.5; # state 4 

mean_lifespan(M7); 

var_lifespan(M7);
# calculations for stage 2 
probs = c(1/3,1/3,1/6,1/6); vals = c(4,3,3,5);  
EX = sum(probs*vals); 
EX2  = sum(probs*vals^2); 
VarX = EX2 - EX^2; 



#######################################################
#### Testing the population moment wrappers
#######################################################
par(mfrow=c(5,4),mar=c(3,3,1,1)); 
for(rep in 1:20) {
	N = 5000*(1 + rpois(3,9)); 
	X1 = rnorm(N[1],4*runif(1)-2,1+5*runif(1)); 
	X2 = rpois(N[2],1 + rpois(10,1)) 
	X3 = c(runif(N[3]/2,-10,-5), runif(N[3]/2,5,10)) 
	X = c(X1,X2,X3); hist(X,main=""); 
	mixdist = N/sum(N);

    popVar = function(x) mean( (x-mean(x))^2 ); 
	popmu3 = function(x) mean( (x - mean(x))^3);   
	mean_by_type = c(mean(X1),mean(X2),mean(X3)); 
    var_by_type = c(popVar(X1),popVar(X2),popVar(X3)); 
	mu3_by_type = c(popmu3(X1), popmu3(X2), popmu3(X3)); 
	
	out = pop_mean_var(mean_by_type,var_by_type,mixdist);
	out2 = pop_mu3(mean_by_type,var_by_type,mu3_by_type,mixdist); 
	
	cat(mean(X), popVar(X), popmu3(X),"\n"); 
	cat(out$pop_mean,out$pop_var,out2,"\n"); 
}
	

	