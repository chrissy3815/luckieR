#########################################################################################
## Function to convert life table columns lx and mx into a Leslie matrix.  
##
## We construct the matrix for a post-breeding census of the population, so that
## individuals currently age a must survive to reach age a+1, at which time 
## they produce m_a offspring which are immediately censused.   
##
## NOTE: if m0>0, the life table has been constructed with nonstandard definitions;
## m0>0 with the usual definitions would imply that newborns immediately have offspring,  
## which in turn immediately have offspring, which in turn...  
##
## By convention, the lx column extends to the maximum age in the model. This can be 
## either the maximum possible age in the population, or the last age at which mx>0, 
## with post-reproductive ages removed from the model. 
#########################################################################################  

lifetable_to_Leslie = function(lx, mx) {

  # Basic checks
  if (!is.numeric(lx) || !is.numeric(mx)) stop("lx and mx must be numeric vectors.")
  n = length(mx); n_lx = length(lx)
  if(n!=n_lx) stop("lx and mx must be the same length"); 

  # Normalize so lx[1] = 1; this allows lx to be counts rather than relative frequencies. 
  if (lx[1] <= 0) stop("lx[1] must be > 0 to normalize.")
  lx = lx/lx[1]

  # Monotonicity check (non-increasing)
  if (any(diff(lx) > 1e-12)) {
    warning("lx is not non-increasing; check your life table.")
  }
  
  # Check that m0 = 0 
  if(mx[1]!=0) stop("You've got tribbles, m_0 is not zero")

  # Now make the matrix ######################################################### 
  
  # Compute survival rates p_x = l_{x+1} / l_x
  px = lx[2:n]/lx[1:(n-1)]

  # Build U subdiagonal
  U = matrix(0, n-1, n-1)
  U[cbind(2:(n-1), 1:(n - 2))] = px[1:(n-2)]; 

  # Build F. For prebreeding census, if the survival of newborns to their first census, 
  # denoted s0, is not supplied, then it is assumed to equal the earliest observed survival. 
  F = matrix(0, n-1, n-1)
  Fx = mx[-1]*px; # live to age (a+1), have m_{a+1} offspring   
  F[1,] = Fx

  return(list(U = U, F = F, L = U + F, px = px))
}
 
lx = c(1.00, 0.8, 0.6, 0.3, 0.1);      # ages 0,1,2,3,4
mx = c(0,    1,    2,   2,   1)
res = lifetable_to_Leslie(lx, mx)
res; 

require(Rage); 
out = mpm_to_table(res$U,res$F,start=1); 
cbind(out$lx,out$mx); ## does not agree. 

