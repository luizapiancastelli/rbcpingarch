#' Bootstrap standard errors of BCP-INGARCH(1,1) CMLEs
#'
#' @export
#' @param model_fit Fitted BCP-INGARCH(1,1) model
#' @param Y observed bivariate count data
#' @param replicas Number of bootstrap replicas
get_boot_se = function(model_fit, Y, replicas){

  A.diag = model_fit$A.diag
  B.diag = model_fit$B.diag

  MLEs = list();
  MLEs$omega= model_fit$omega
  MLEs$phi = model_fit$phi
  if(A.diag == TRUE & B.diag == TRUE){
    MLEs$A = diag(model_fit$A)
    MLEs$B =  diag(model_fit$B)
    model  = stanmodels$modelABdiag
  } else if(A.diag == TRUE & B.diag == FALSE){
    MLEs$A = diag(model_fit$A)
    MLEs$B = matrix(model_fit$B, ncol = 2)
    model = stanmodels$modelAdiag
  } else if(A.diag == FALSE & B.diag == TRUE){
    model = stanmodels$modelBdiag
    
  } else {
    model = stanmodels$modelfull
    
  }

  boot_estimates = furrr::future_map(1:replicas, ~ boot_replica(MLEs, nrow(Y), model , model_fit$initial_values))
  boot_estimates = do.call(rbind, boot_estimates)
  boot_sd = apply(boot_estimates,2, sd)
  return(boot_sd)

}

boot_replica = function(MLEs, n, inits){

  data_boot = rBCPINGARCH(MLEs$A,MLEs$B,MLEs$omega, MLEs$phi, n)
  Y_boot = data_boot$Y
  inits$omega = matrix(colMeans(Y_boot), ncol = 1)
  op_boot = optimizing(model, data = list(Y = Y_boot, n = nrow(Y_boot), med_y1 = median(Y_boot[,1]), med_y2 = median(Y_boot[,2])),
                       verbose = FALSE,init = inits)
  if(op_boot$return_code != 0){

    while(op_boot$return_code != 0){
      data_boot = rBCPINGARCH(MLEs$A,MLEs$B,MLEs$omega, MLEs$phi, n)
      Y_boot = data_boot$Y
      inits$omega = matrix(colMeans(Y_boot), ncol = 1)
      op_boot = rstan::optimizing(model, data = list(Y = Y_boot, n = nrow(Y_boot), med_y1 = median(Y_boot[,1]), med_y2 = median(Y_boot[,2])),
                           verbose = FALSE,init = inits)
    }
  }
  par_names = names(op_boot$par)
  A_boot = op_boot$par[str_detect(par_names, "A\\[" )]
  B_boot = op_boot$par[str_detect(par_names, "B\\[" )]
  phi_boot =op_boot$par["phi"]
  omega_boot = op_boot$par[str_detect(par_names, "omega" )]

  return(c(A_boot, B_boot, omega_boot, phi_boot))

}

