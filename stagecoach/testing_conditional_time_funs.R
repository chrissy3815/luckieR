
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
	

	
	