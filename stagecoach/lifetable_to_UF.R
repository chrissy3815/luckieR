#########################################################################################
## Function to convert life table columns lx and mx into a Leslie matrix 
## For prebreeding census, s0 is the survival from birth to first census of newborns. 
## If not supplied, s0 is assumed to be lx[2]/lx[1], the earliest observed survival. 
#########################################################################################  

lifetable_to_UF = function(lx, mx, census = c("postbreeding", "prebreeding"),s0=NULL) {

  # Basic checks
  if (!is.numeric(lx) || !is.numeric(mx)) stop("lx and mx must be numeric vectors.")
  n = length(mx); n_lx = length(lx)
  # if(n!=n_lx) stop("lx and mx must be the same length"); 

  # Normalize so lx[1] = 1; this allows lx to be counts rather than relative frequencies. 
  if (lx[1] <= 0) stop("lx[1] must be > 0 to normalize.")
  lx = lx/lx[1]

  # Monotonicity check (non-increasing)
  if (any(diff(lx) > 1e-12)) {
    warning("lx is not non-increasing; check your life table.")
  }

  # Compute all adjacent survival ratios p_all = l_{x+1} / l_x over what's available
  p_all = lx[-1] / lx[-length(lx)]
  p_all[!is.finite(p_all)] = 0

  # Build U subdiagonal
  U = matrix(0, n, n)
  if (n >= 2) {
     if (length(p_all) < n - 1) stop("lx is too short to compute U")
     U[cbind(2:n, 1:(n - 1))] = p_all
  }

  # Build F. For prebreeding census, if the survival of newborns to their first census, 
  # denoted s0, is not supplied, then it is assumed to equal the earliest observed survival. 
  F = matrix(0, n, n)
  if (census == "postbreeding") {
    F[1, ] = mx
  } else {
	if(is.null(s0)) s0 = p_all[1] 
    F[1, ] = mx*p_all[1]
  }

  return(list(U = U, F = F, L = U + F, px = p_all))
}
 
lx = c(1.00, 0.75, 0.55, 0.30, 0.12)     # ages 0,1,2,3,4
mx = c(0.00, 0.10, 0.60, 1.20, 0.50)
res = lifetable_to_UF(lx, mx, census = "prebreeding")
res; 

