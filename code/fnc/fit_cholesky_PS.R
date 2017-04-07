

fit_cholesky_PS <- function(yVec,U,
                            P_l,lambda_l,
                            P_m,lambda_m,
                            lambda_ridge){
      Pen <- rbind(lambda_l*P_l,
                   lambda_m*P_m,
                   lambda_ridge*diag(ncol(P_l)))
      n.col <- ncol(U)
      nix <- rep(0,nrow(Pen))
      
      nix.ridge <- rep(0, ncol(U))
      coef.est <- rep(1, ncol(U))
      
      mu <- rep(mean(yVec), length(yVec))
      it <- 0
      repeat {
            if(it == 0) {
                  eta <- mu
            }
            
            it <- it + 1
            if(it > 25)
                  break
            
            mu <- eta
            h.prime <- 1
            w <- rep(1, length(y_vec))
            u <- (y_vec - mu)/h.prime + eta
            
            startTS <- Sys.time()
            f <- lsfit(rbind(U,Pen), c(u, nix), wt = c(w, nix + 1) *c(w, nix + 1), intercept = F)
            endTS <- Sys.time()
            endTS-startTS
            
            coef.old <- coef.est
            coef.est <- as.vector(f$coef)
            
            d.coef <- max(abs((coef.est[coef.old>0] - coef.old[coef.old>0])/coef.old[coef.old>0]))
            if(d.coef < 1e-008)
                  break
            print(c(it, d.coef))
            eta <- U %*% coef.est
      }
      
      if(it > 24) {
            warning(paste("parameter estimates did NOT converge in 25 iterations"
            ))
      }      
      
      # H <- hat(f$qr, intercept = F)[1:(m-1)]
      # trace <- eff.dim <- sum(H)
      # 
      # cv <- press.mu <- press.e <- var.c <- NULL
      # dev <- sum(f$residuals[1:(m-1)]^2)
      # dispersion.parm <- dev/((m-1) - trace)
      # press.e <- f$residuals[1:(m-1)]/(1 - H)
      # cv <- sqrt(sum((press.e)^2)/(m-1))
      # press.mu <- y_vec - press.e
      # 
      # aic <- dev + 2 * trace
      # aic
      #f$coefficients      
      f
}
