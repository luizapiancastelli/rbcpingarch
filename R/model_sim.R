#' Simulate a Bivariate Conditional Poisson INGARCH(1,1) process
#'
#' @export
#' @param A 2x2 matrix of parameters
#' @param B 2x2 matrix of parameters
#' @param omega
#' @param phi
#' @param n
rBCPINGARCH = function(A,B, omega, phi, n){

  burn_in = 300
  Y = matrix(NA, ncol =2, nrow = n+burn_in)
  cor_matrix = matrix(NA, ncol =1, nrow = n+burn_in)
  lambda = matrix(NA, ncol =2, nrow = n+burn_in)

  l_t1 = solve(diag(2) - A - B)%*%omega
  Y[1,] = sim_pair(l_t1[1],l_t1[2], phi)

  lambda[1,] = l_t1

  for(t in 2:(n+burn_in)){
    lambda[t, ] = omega + A%*%lambda[t-1,] + B%*%Y[t-1, ]
    Y[t,] = sim_pair(lambda[t,1],lambda[t,2], phi)
  }

  return(Y[(burn_in+1):(n+burn_in),])
}

sim_pair= function(l1, l2, phi){
  y1_i = rpois(1, l1)
  mu2 = l2*exp(-l1*(exp(phi)-1))
  y2_i = rpois(1, mu2*exp(phi*y1_i))
  return(c(y1_i, y2_i))
}

