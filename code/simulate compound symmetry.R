

library(MASS)

m <- 30
grid <- expand.grid(s=(1:m),t=(1:m)) %>%
      subset(.,s>t) %>%
      transform(l=(s-t),m=s+t) %>%
      transform(.,phi=(l/max(grid$l))^2)


Sigma <- matrix(0.7,nrow=m,ncol=m) + diag(rep(0.3),m)
Omega <- solve(Sigma)


C <- t(chol(Sigma))
D <- diag(diag(C))
L <- C%*%solve(D)
T_mod <- solve(L)



matplot(Y_chol <- L %*% t(mvrnorm(n=100,mu=rep(0,m),Sigma=I(D^2))),col="pink",type="l")
matlines(t(mvrnorm(n=100,mu=rep(0,m),Sigma=Sigma)),col="blue")

grid <- expand.grid(t=1:m,s=1:m) %>%
      transform(l=t-s,m=s+t)

persp(grid$s, grid$t, z=diag(rep(1,nrow(T_mod)))-T_mod, col = color[facetcol], phi = 30, theta = -30)



      data.frame(grid,
                 phi=as.vector(diag(rep(1,nrow(T_mod)))-T_mod)) %>%
            subset(.,t>s) %>%
      orderBy(~t+s,.) %>%
      wireframe(phi~t+s,
                data=.,
                screen=list(z=-10,x=-65),
                light.source = c(5,20,10),
                pretty=TRUE,
                scales = list(arrows = FALSE),
                drape=FALSE,
                cex=0.15,
                col="grey",
                par.settings = list(axis.line = list(col = "transparent")))

      
      ggrid <-  data.frame(grid,
                           phi=as.vector(diag(rep(1,nrow(T_mod)))-T_mod))
      ggrid$phi[ggrid$l<=0] <- 0
      plot(subset(ggrid, m==31)$l,subset(ggrid, m==31)$phi,type="l",
           xlab="l",ylab=expression(phi))
      for (grid.m in 16:46){
            lines(subset(ggrid, m==grid.m)$l,subset(ggrid, m==grid.m)$phi,
                  col=grid.m)      
      }
      plot(subset(ggrid, l==1)$m,subset(ggrid, l==1)$phi,type="l",
           xlab="l",ylab=expression(phi))
      for (grid.l in 2:31){
            lines(subset(ggrid, l==grid.l)$m,subset(ggrid, l==grid.l)$phi,
                  col=grid.l)      
      }
      
      transform(expand.grid(t=1:30,s=1:30),
                l=t-s,m=s+t)
      
      wireframe(phi~s+t,
                data=ggrid,
                light.source = c(5,20,10),
                pretty=TRUE,
                scales = list(arrows = FALSE),
                drape=FALSE,
                screen=list(z = 80, x = -70, y = 3),
                row.values=transform(expand.grid(t=1:30,s=1:30),
                                           l=t-s,m=s+t)$s,
                col.values=transform(expand.grid(t=1:30,s=1:30),
                                     l=t-s,m=s+t)$t,
                cex=0.15,
                col="grey",
                par.settings = list(axis.line = list(col = "transparent")))
      
      
      
      
      
      
      
      
      
      
      library(magrittr)
      library(rlist)
      library(plyr)
      library(dplyr)
      M <- m <- 30
      N <- 50
      grid <- expand.grid(t=1:M,s=1:M) %>% subset(.,t>s) %>%
            transform(.,l=t-s,
                      m=(t+s)/2)
      
      bPars <- rbind(c(0,1,50,3,100,2),
                     c(0,1,60,3,100,3))
      
      y <- mvrnorm(n=N,mu=rep(0,m),Sigma=Sigma)
      y_vec <- as.vector(t(y[,-1]))
      
      X <- matrix(data=0,nrow=N*(M-1),ncol=choose(M,2))
      no.skip <- 0
      for (t in 2:M){
            X[((0:(N-1))*(M-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- y[,1:(t-1)]
            no.skip <- no.skip + t - 1
      }
      
      
      
      
      
      l. <- as.vector(outer(rep(1,length(unique(grid$m))),unique(grid$l)))
      Bl <- bsplbase(l./max(grid$l), bPars[1,  ])$base
      knots.l <- bsplbase(l./max(grid$l), bPars[1,  ])$knots
      knots.l <- knots.l[1:(length(knots.l)-(bPars[1,4]+1))]
      
      
      m. <- as.vector(outer(unique(grid$m),rep(1,length(unique(grid$l)))))
      Bm <- bsplbase(m./max(grid$m), bPars[2,  ])$base
      knots.m <- bsplbase(m./max(grid$m), bPars[2,  ])$knots
      knots.m <- knots.m[1:(length(knots.m)-(bPars[2,4]+1))]
      
      knot_grid <- expand.grid(m=knots.m,l=knots.l)[,2:1]
      keep <- (knot_grid$m >= 0.5) & (knot_grid$m < max(knot_grid$l)-0.5*knot_grid$l) |
            (knot_grid$m < 0.5) & (knot_grid$m > min(knot_grid$l)+0.5*knot_grid$l)
      knot_grid <- transform(knot_grid,keep=factor(keep))
      knot_grid <- data.frame(knot_grid,
                              expand.grid(m_index=1:length(knots.m),
                                          l_index=(1:length(knots.l)))[,2:1])
      
      
      
      
      n1 <- ncol(Bl)
      n2 <- ncol(Bm)  # Compute tensor products
      
      Bl. <- kronecker(Bl, t(rep(1, n2)))
      Bm. <- kronecker(t(rep(1, n1)), Bm)
      
      
      grid <- orderBy(~ l+m,grid)
      discard_rows <- left_join(data.frame(l=l.,m=m.,frame="big"),
                                data.frame(grid[,c("l","m")],frame="little")
                                ,by=c("l","m"))$frame.y %>% is.na %>% which
      
      B. <- Bl.*Bm.
      basisKeepIndex <- knot_grid$keep[!duplicated(knot_grid[,c("l_index","m_index")])] %>% equals(TRUE) %>% which
      B. <- B.[-discard_rows,basisKeepIndex]
      
      
      
      
      
      
      
      
      
      knot_grid <- subset(knot_grid,keep==TRUE)
      dl <- bPars[1, 6]
      Pl <- matrix(data=0,nrow=sum(knot_grid$m==unique(knot_grid$m)[1])-dl,
                   ncol=nrow(knot_grid))
      Pl[,which(knot_grid$m==unique(knot_grid$m)[1])] <- diff(diag(sum(knot_grid$m==unique(knot_grid$m)[1])),
                                                              differences = dl)
      for(i in 2:length(unique(knot_grid$m))){
            pl <- matrix(data=0,nrow=sum(knot_grid$m==unique(knot_grid$m)[i])-dl,
                         ncol=nrow(knot_grid))
            pl[,which(knot_grid$m==unique(knot_grid$m)[i])] <- diff(diag(sum(knot_grid$m==unique(knot_grid$m)[i])),
                                                                    differences = dl)      
            Pl <- rbind(Pl,pl)
      }
      lambdal <- bPars[1, 5]
      
      
      
      
      
      
      
      
      
      dm <- bPars[2, 6]
      Pm <- matrix(data=0,nrow=sum(knot_grid$l==unique(knot_grid$l)[1])-dm,
                   ncol=nrow(knot_grid))
      Pm[,which(knot_grid$l==unique(knot_grid$l)[1])] <- diff(diag(sum(knot_grid$l==unique(knot_grid$l)[1])),
                                                              differences = dm)
      for(i in 2:length(unique(knot_grid$l))){
            pm <- matrix(data=0,nrow=sum(knot_grid$l==unique(knot_grid$l)[i])-dm,
                         ncol=nrow(knot_grid))
            pm[,which(knot_grid$l==unique(knot_grid$l)[i])] <- diff(diag(sum(knot_grid$l==unique(knot_grid$l)[i])),
                                                                    differences = dm)      
            Pm <- rbind(Pm,pm)
      }
      lambdam <- bPars[2, 5]
      
      
      U. <- X%*%B.
      
      
      
      
      
      
      
      knot_grid <- subset(knot_grid,keep==TRUE)
      dl <- bPars[1, 6]
      if(dl>0){
            if(sum(knot_grid$m==unique(knot_grid$m)[1])>dl){
                  Pl <- matrix(data=0,nrow=sum(knot_grid$m==unique(knot_grid$m)[1])-dl,
                               ncol=nrow(knot_grid))
                  Pl[,which(knot_grid$m==unique(knot_grid$m)[1])] <- diff(diag(sum(knot_grid$m==unique(knot_grid$m)[1])),
                                                                          differences = dl)            
                  for(i in 2:length(unique(knot_grid$m))){
                        pl <- matrix(data=0,nrow=sum(knot_grid$m==unique(knot_grid$m)[i])-dl,
                                     ncol=nrow(knot_grid))
                        pl[,which(knot_grid$m==unique(knot_grid$m)[i])] <- diff(diag(sum(knot_grid$m==unique(knot_grid$m)[i])),
                                                                                differences = dl)      
                        Pl <- rbind(Pl,pl)
                  }
            }
            
            if((sum(knot_grid$m==unique(knot_grid$m)[1])<=dl) &(sum(knot_grid$m==unique(knot_grid$m)[2])>dl)){
                  Pl <- matrix(data=0,nrow=sum(knot_grid$m==unique(knot_grid$m)[2])-dl,
                               ncol=nrow(knot_grid))
                  Pl[,which(knot_grid$m==unique(knot_grid$m)[2])] <- diff(diag(sum(knot_grid$m==unique(knot_grid$m)[2])),
                                                                          differences = dl)            
                  for(i in 3:length(unique(knot_grid$m))){
                        if(sum(knot_grid$m==unique(knot_grid$m)[i])>dl){
                              pl <- matrix(data=0,nrow=sum(knot_grid$m==unique(knot_grid$m)[i])-dl,
                                           ncol=nrow(knot_grid))
                              pl[,which(knot_grid$m==unique(knot_grid$m)[i])] <- diff(diag(sum(knot_grid$m==unique(knot_grid$m)[i])),
                                                                                      differences = dl)      
                              Pl <- rbind(Pl,pl)      
                        }
                  }
            }
            
      }
      if(dl==0){
            Pl <- diag(nrow(knot_grid))
      }
      lambdal <- bPars[1, 5]
      
      
      
      
      
      
      
      
      
      dm <- bPars[2, 6]
      dm <- 2
      if(dm>0){
            if((sum(knot_grid$m==unique(knot_grid$m)[1])>dl)){
                  Pm <- matrix(data=0,nrow=sum(knot_grid$l==unique(knot_grid$l)[1])-dm,
                               ncol=nrow(knot_grid))
                  Pm[,which(knot_grid$l==unique(knot_grid$l)[1])] <- diff(diag(sum(knot_grid$l==unique(knot_grid$l)[1])),
                                                                          differences = dm)
                  for(i in 2:length(unique(knot_grid$l))){
                        if(sum(knot_grid$m==unique(knot_grid$m)[i])>dl){
                              pm <- matrix(data=0,nrow=sum(knot_grid$l==unique(knot_grid$l)[i])-dm,
                                           ncol=nrow(knot_grid))
                              pm[,which(knot_grid$l==unique(knot_grid$l)[i])] <- diff(diag(sum(knot_grid$l==unique(knot_grid$l)[i])),
                                                                                      differences = dm)      
                              Pm <- rbind(Pm,pm)
                        }
                  }  
            }
            
            if((sum(knot_grid$m==unique(knot_grid$m)[1])<=dl) &(sum(knot_grid$m==unique(knot_grid$m)[2])>dl)){
                  Pm <- matrix(data=0,nrow=sum(knot_grid$l==unique(knot_grid$l)[2])-dm,
                               ncol=nrow(knot_grid))
                  Pm[,which(knot_grid$l==unique(knot_grid$l)[2])] <- diff(diag(sum(knot_grid$l==unique(knot_grid$l)[2])),
                                                                          differences = dm)
                  for(i in 3:length(unique(knot_grid$l))){
                        if(sum(knot_grid$m==unique(knot_grid$m)[i])>dl){
                              pm <- matrix(data=0,nrow=sum(knot_grid$l==unique(knot_grid$l)[i])-dm,
                                           ncol=nrow(knot_grid))
                              pm[,which(knot_grid$l==unique(knot_grid$l)[i])] <- diff(diag(sum(knot_grid$l==unique(knot_grid$l)[i])),
                                                                                      differences = dm)      
                              Pm <- rbind(Pm,pm)
                        }
                  }  
            }
      }
      if(dm==0){
            Pm <- diag(nrow(knot_grid))
      }
      
      lambdam <- bPars[2, 5]
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
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
            f$coefficients      
      }
      
      lambdas <- expand.grid(lam_l=exp(seq(-2,7,length.out=15)),
                               lam_m=exp(seq(-2,7,length.out=15)))
      
      coef_list <- list.zip(lam_l=lambdas$lam_m,
                              lam_m=lambdas$lam_l) %>%
            lapply(.,function(l){
                  fit_cholesky_PS(y_vec,U.,Pl,l$lam_l,
                                  Pm,l$lam_m,
                                  0.12)    
            })
      
      setwd("/Users/taylerblake/Documents/Dissertation/code/simulations")
      save(coef_list,file="cpdSymm_coef_list_2.Rdata")
      
      
      
      
      
      
      
      gg <- expand.grid(s=(1:m),t=(1:m)) %>% subset(.,t>s)
      gg <- transform(gg,l=t-s,m=s+t)
      
      l. <- gg$l
      Bl <- bsplbase(l./max(gg$l), bPars[1,  ])$base
      knots.l <- bsplbase(l./max(gg$l), bPars[1,  ])$knots
      knots.l <- knots.l[1:(length(knots.l)-(bPars[1,4]+1))]
      
      m. <- gg$m
      #m. <- as.vector(outer(unique(grid$m),rep(1,length(unique(grid$l)))))
      Bm <- bsplbase(m./max(gg$m), bPars[2,  ])$base
      knots.m <- bsplbase(m./max(gg$m), bPars[2,  ])$knots
      knots.m <- knots.m[1:(length(knots.m)-(bPars[2,4]+1))]
      
      knot_grid <- expand.grid(m=knots.m,l=knots.l)[,2:1]
      keep <- (knot_grid$m >= 0.5) & (knot_grid$m < max(knot_grid$l)-0.5*knot_grid$l) |
            (knot_grid$m < 0.5) & (knot_grid$m > min(knot_grid$l)+0.5*knot_grid$l)
      knot_grid <- transform(knot_grid,keep=factor(keep))
      knot_grid <- data.frame(knot_grid,
                              expand.grid(m_index=1:length(knots.m),
                                          l_index=(1:length(knots.l)))[,2:1])
      
      B_lm <- t(sapply(1:nrow(Bl),function(i){as.vector(outer(Bm[i,],Bl[i,]))}))
      basisKeepIndex <- knot_grid$keep[!duplicated(knot_grid[,c("l_index","m_index")])] %>% equals(TRUE) %>% which
      knot_grid <- subset(knot_grid,keep==TRUE)
      B_lm <- B_lm[,basisKeepIndex]
      
      
      
      
      
      
      phi_list <- lapply(coef_list,function(beta){
            B_lm%*%beta
      })
      i <- 1
      
      
      #Phi_hat <- B_lm %*% f$coefficients
     # par(ask=TRUE)
      #for (i in 76:nrow(lambdas)){
            laml <- round(lambdas$lam_l[i],3)
            lamm <- round(lambdas$lam_m[i],3)
            my_title <- c(as.expression(bquote(lambda[l] == .(laml))), as.expression(bquote(lambda[m] == .(lamm))))
            Phi_hat <- phi_list[[i]]
            data.frame(phi=Phi_hat,gg) %>%
                  subset(.,l<(max(gg$l)-2)) %>%
                  wireframe(phi~s*t,
                            data=.,
                            screen=list(z=55,x=-85,y=10),
                            light.source = c(5,20,10),
                            pretty=TRUE,
                            scales = list(arrows = FALSE),
                            drape=FALSE,
                            cex=0.15,
                            col="grey",
                            par.settings = list(axis.line = list(col = "transparent")),
                            main=my_title
                  )            
     # }
      i <- i+1
      
