


fit_monotone_PS <-  function(Y,
                             U,D,
                             P_l,lambda_l,
                             P_m,lambda_m,
                             lambda_ridge,
                             P_l_shape,lambda_l_shape,
                             max_iter){
  t <- 0
  W <- list()
  W[[t+1]] <- matrix(data=0,nrow=nrow(P_l_shape),ncol=nrow(P_l_shape))
  ## allocate storage for sequence of coefficients
  beta_mat <- matrix(data=NA,
                     nrow=ncol(B.),
                     ncol=max.iter)
  fit_list <- list()
  ## estimate beta[t] on W=W[t-1]
  not_converged <- TRUE
  while (not_converged & t < max.iter) {
    t <- t+1
    fit_list[[t]] <- fit_cholesky_PS(Y,
                                     y_aug=NULL,
                                     U,
                                     D,
                                     Pl,
                                     lambda_l,
                                     Pm,
                                     lambda_m,
                                     P_l_shape=W[[t]]%*%P_l_shape,
                                     lambda_l_shape,
                                     lambda_ridge=0)
    beta_mat[,t] <- fit_list[[t]]$coef
    W[[t+1]] <- diag(as.numeric((Pl_shape%*%beta_mat[,t]) > 0))
    
    ## if the weight vector for the monotonicity penalty doesnt change, exit
    if ( min(diag(W[[t+1]])==diag(W[[t]])) == 1 ) {
      not_converged <- FALSE
    }
  }
  
  #beta_mat <- beta_mat[,colSums(is.na(beta_mat))==0]
  #beta_hat <- beta_mat[,ncol(beta_mat)]
  fit <- fit_list[[t]]
  fit
}


