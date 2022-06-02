#' Fit a Bivariate Conditional Poisson INGARCH(1,1) model
#'
#' @export
#' @param Y integer matrix of bivariate count data.
#' @param A.diag boolean specifying if the A matrix should be diagonal. Defaults to TRUE.
#' @param B.diag boolean specifying if the B matrix should be diagonal. Defaults to FALSE.
#' @return Named list containing the parameter estimates $par, standard errors estimated from the hessian matrix $se,
#' initial values used by the algorithm $initial_values, and value of optimized log-likelihood function.
fit_BCP_INGARCH = function(Y, A.diag = TRUE, B.diag = FALSE){
  inits = list();
  inits$omega = matrix(colMeans(Y), ncol = 1); inits$phi = 0;

  if(A.diag == TRUE & B.diag == FALSE){
    inits$A =  rep(0.1,2)
    inits$B = matrix(rep(0.1,4),ncol = 2)
    op = rstan::optimizing(stanmodels$modelAdiag, data = list(Y = Y, n = nrow(Y), med_y1 = median(Y[,1]), med_y2 = median(Y[,2])),verbose = FALSE, init = inits)
    par = op$par;
    a_names = c('A11', 'A22')
    b_names = c('B11', 'B21', 'B12', 'B22')

  } else if(A.diag == TRUE & B.diag == TRUE){

    inits$A = rep(0.1,2)
    inits$B = rep(0.1,2)
    op = rstan::optimizing(stanmodels$modelABdiag, data = list(Y = Y, n = nrow(Y),med_y1 = median(Y[,1]), med_y2 = median(Y[,2])),verbose = FALSE,init = inits)
    par = op$par;
    a_names = c('A11', 'A22')
    b_names = c('B11', 'B22')

  } else if(A.diag == FALSE & B.diag == TRUE){

    inits$A = matrix(rep(0.1,4),ncol = 2)
    inits$B = rep(0.1,2)
    op = rstan::optimizing(stanmodels$modelBdiag , data = list(Y = Y, n = nrow(Y),med_y1 = median(Y[,1]), med_y2 = median(Y[,2])),verbose = FALSE,init = inits)
    par = op$par;
    a_names = c('A11', 'A21', 'A12', 'A22')
    b_names = c('B11', 'B22')

  } else if(A.diag == FALSE & B.diag == FALSE){

    inits$A = inits$B = matrix(rep(0.1,4),ncol = 2)
    op = rstan::optimizing(stanmodels$modelfull, data = list(Y = Y, n = nrow(Y),med_y1 = median(Y[,1]), med_y2 = median(Y[,2])),verbose = FALSE,init = inits)
    par = op$par;
    a_names = c('A11', 'A21', 'A12', 'A22')
    b_names = c('B11', 'B21', 'B12', 'B22')

  }
  par_names = names(op$par)
  A = op$par[stringr::str_detect(par_names, "A\\[" )]
  B = op$par[stringr::str_detect(par_names, "B\\[" )]
  phi =op$par["phi"]
  omega = op$par[stringr::str_detect(par_names, "omega" )]
  lambda =  op$par[stringr::str_detect(par_names, "lambda" )]

  hes = numDeriv::hessian(loglik_BCP_INGARCH, c(A, B, omega, phi), Y=Y, B.diag = B.diag)
  ses = sqrt(-diag(solve(hes)))
  names(ses) = c(a_names, b_names, 'omega1', 'omega2', 'phi')

  output = list();
  output$par = c(A, B, omega, phi)
  output$se = ses
  output$lambda= matrix(lambda, ncol = 2)
  output$initial_values = inits
  output$A.diag = A.diag;
  output$B.diag = B.diag;
  output$loglik = op$value
  return(output)
}


#' Log-likelihod value from model fit
loglik_BCP_INGARCH = function(params, Y, A.diag, B.diag){

  n = nrow(Y);
  A = params[stringr::str_detect(names(params), "A" )]
  B = params[stringr::str_detect(names(params), "B" )]

  if(length(A)==2){
    A_mat = diag(A);
  } else {
    A_mat = matrix(A, ncol =2)
  }
  if(length(B)==2){
    B_mat = diag(B);
  } else {
    B_mat = matrix(B, ncol =2)
  }
  omega = params[stringr::str_detect(names(params), "omega" )]
  phi =  params['phi']
  lambda = matrix(NA, ncol =2, nrow = n)
  lambda[1,1] = omega[1] + A_mat[1,1]*mean(Y[,1]) + B_mat[1,1]*median(Y[,1]) + B_mat[1,2]*median(Y[,2])
  lambda[1,2] = omega[2] + A_mat[2,2]*mean(Y[,2]) + B_mat[2,2]*median(Y[,2]) + B_mat[2,1]*median(Y[,1])

  for(t in 2:n){
    lambda[t,] = omega + A_mat%*%lambda[t-1,] + B_mat%*%as.numeric(Y[t-1, ])
  }

  y1 = Y[2:n,1]; y2 = Y[2:n,2];
  l1 = lambda[2:n, 1]; l2 = lambda[2:n, 2];
  loglik = sum(y1*log(l1)) +  sum(y2*log(l2))
  loglik = loglik  - sum(l1*(1 + y2*(exp(phi) -1)) )
  loglik = loglik - sum(l2*exp( -l1*(exp(phi)-1) + phi*y1))
  loglik = loglik +  phi*sum(y1*y2)
  names(loglik) = NULL
  return(loglik)
}


#' Model information criteria for BCP-INGARCH(1,1) fits.
#'
#' @param fit fitted model obtained with the function fit_BCP_INGARCH. 
#' @return AIC anc BIC of fitted models.
#' @export
model_information = function(fit){
  
  params_A = ifelse(fit$A.diag,2,4)
  params_B = ifelse(fit$B.diag,2,4)
  k = params_A + params_B + 3
  n = nrow(fit$lambda)
  aic = 2*k - 2*fit$loglik
  bic = k*log(n)- 2*fit$loglik
  
  criteria = c(aic, bic)
  names(criteria) = c("AIC", "BIC")
  return(criteria)
}

