
library(plyr)
library(dplyr)
library(doBy)
library(tidyr)
library(splines)
library(magrittr)
library(ggplot2)









###########################################################################################################
###########################################################################################################


M <- m <- 20

grid <- expand.grid(t=1:M,s=1:M) %>% subset(.,t>s) %>%
      transform(.,l=t-s,
                m=(t+s)/2)

bPars <- rbind(c(0,1,50,3,100,2),
               c(0,1,60,3,100,3))



phi <- function(t1.index,t2.index,m){
      if(t1.index-t2.index==1){
            garp <- (2*((t1.index)/m)^2)-0.5
      }
      if(t1.index-t2.index!=1){
            garp <- 0
      }
      garp
}

true_phi <- sapply(1:nrow(grid),function(row.i){
      phi(grid$t[row.i],grid$s[row.i],m=m)
})

grid <- transform(grid,true_phi=true_phi)
grid$true_phi


indices_of_nonzeros <- as.matrix(expand.grid(t=(1:m),s=(1:m)) %>% subset(.,(t-s)==1))
nonzero_phis <- (2*((2:m)/m)^2)-0.5

T_mat <- diag(rep(1,m))
phis <- as.vector(rep(0,sum(lower.tri(T_mat))))
T_mat[indices_of_nonzeros] <- -nonzero_phis
phis <- -T_mat[lower.tri(T_mat)]


N <- 3


y <- t(solve(T_mat)%*%matrix(data=rnorm(N*m,mean=0,sd=1),
                             nrow=m,
                             ncol=N))
#iv_errors <- rep(0.1,m)
#D <- diag(iv_errors)
D <- diag(m) 
     
true_Sigma <- solve(T_mat)%*%D%*%t(solve(T_mat))
true_Omega <- t(T_mat)%*%solve(D)%*%T_mat


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








reference_grid <- data.frame(l=l.,m=m.)[-discard_rows,]
match_row <- function(row,df){
      match.ind <- match(df[,1],row[1])
      for(j in 2:ncol(df)){
            match.ind <- match.ind*match(df[,j],row[j])
      }
      which(!is.na(match.ind))
}
x_col_ind <-   sapply(1:nrow(grid),function(row.i){
      match_row(grid[row.i,c("l","m")],reference_grid)
})
grid$column_index <- x_col_ind

X <- matrix(data=0,nrow=N*(m-1),ncol=nrow(reference_grid))
for(row.i in 1:nrow(grid)){
      X[((0:(N-1))*(M-1)) +grid$t[row.i]-1,grid$column_index[row.i]] <- as.vector(y)[grid$s[row.i]]
}

U. <- X%*%B.















Pen <- rbind(lambdal*Pl, lambdam*Pm)
n.col <- ncol(U.)
nix <- rep(0,nrow(Pen))

nix.ridge <- rep(0, ncol(U.))
coef.est <- rep(1, ncol(U.))

mu <- rep(mean(y_vec), length(y_vec))
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
      #w <- rep(1, length(y[-1]))
      w <- rep(10, length(y_vec))
      #u <- (y[-1] - mu)/h.prime + eta
      u <- (y_vec - mu)/h.prime + eta

      startTS <- Sys.time()
      f <- lsfit(rbind(U.,Pen), c(u, nix), wt = c(w, nix + 1) *c(w, nix + 1), intercept = F)
      endTS <- Sys.time()
      endTS-startTS

      coef.old <- coef.est
      coef.est <- as.vector(f$coef)

      d.coef <- max(abs((coef.est[coef.old>0] - coef.old[coef.old>0])/coef.old[coef.old>0]))
      if(d.coef < 1e-008)
            break
      print(c(it, d.coef))
      eta <- U. %*% coef.est
      }

      if(it > 24) {
            warning(paste("parameter estimates did NOT converge in 25 iterations"
            ))
      }

H <- hat(f$qr, intercept = F)[1:(m-1)]
trace <- eff.dim <- sum(H)

cv <- press.mu <- press.e <- var.c <- NULL
dev <- sum(f$residuals[1:(m-1)]^2)
dispersion.parm <- dev/((m-1) - trace)
press.e <- f$residuals[1:(m-1)]/(1 - H)
cv <- sqrt(sum((press.e)^2)/(m-1))
press.mu <- y_vec - press.e

aic <- dev + 2 * trace
aic








##-----------------------------------------------------------------------------


gg <- expand.grid(s=(1:50),t=(1:50)) %>% subset(.,t>s)
gg <- transform(gg,l=t-s,m=s+t)

l. <- as.vector(outer(rep(1,length(unique(gg$m))),unique(gg$l)))
Bl <- bsplbase(l./max(gg$l), bPars[1,  ])$base
knots.l <- bsplbase(l./max(gg$l), bPars[1,  ])$knots
knots.l <- knots.l[1:(length(knots.l)-(bPars[1,4]+1))]


m. <- as.vector(outer(unique(gg$m),rep(1,length(unique(gg$l)))))
Bm <- bsplbase(m./max(gg$m), bPars[2,  ])$base
knots.m <- bsplbase(m./max(gg$m), bPars[2,  ])$knots
knots.m <- knots.m[1:(length(knots.m)-(bPars[2,4]+1))]

knot_gg <- expand.grid(m=knots.m,l=knots.l)[,2:1]
keep <- (knot_gg$m >= 0.5) & (knot_gg$m < max(knot_gg$l)-0.5*knot_gg$l) |
      (knot_gg$m < 0.5) & (knot_gg$m > min(knot_gg$l)+0.5*knot_gg$l)
knot_gg <- transform(knot_gg,keep=factor(keep))
knot_gg <- data.frame(knot_gg,
                      expand.grid(m_index=1:length(knots.m),
                                  l_index=(1:length(knots.l)))[,2:1])

Bl. <- kronecker(Bl, t(rep(1, n2)))
Bm. <- kronecker(t(rep(1, n1)), Bm)

gg <- orderBy(~ l+m,gg)
discard_rows <- left_join(data.frame(l=l.,m=m.,frame="big"),
                          data.frame(gg[,c("l","m")],frame="little")
                          ,by=c("l","m"))$frame.y %>% is.na %>% which

B. <- Bl.*Bm.
basisKeepIndex <- knot_grid$keep[!duplicated(knot_gg[,c("l_index","m_index")])] %>% equals(TRUE) %>% which
B. <- B.[-discard_rows,basisKeepIndex]







Phi_hat <- B. %*% f$coefficients
data.frame(phi=Phi_hat,gg) %>%
      wireframe(phi~s*t,
                data=.,
                screen=list(z=80,x=-75),
                light.source = c(5,20,10),
                pretty=TRUE,
                scales = list(arrows = FALSE),
                drape=FALSE,
                cex=0.15,
                col="grey",
                par.settings = list(axis.line = list(col = "transparent")))







true_phi <- 2*(gg$t/max(gg$t))^2-0.5
true_phi[gg$l >1] <- 0
data.frame(true_phi=true_phi,gg) %>%
      wireframe(true_phi~s*t,
                data=.,
                screen=list(z=80,x=-75),
                light.source = c(5,20,10),
                pretty=TRUE,
                scales = list(arrows = FALSE),
                drape=FALSE,
                cex=0.15,
                col="grey",
                par.settings = list(axis.line = list(col = "transparent")))







###########################################################################################################
###########################################################################################################
###########################################################################################################






















l. <- grid$l
Bl <- bsplbase(l./max(grid$l), bPars[1,  ])$base
knots.l <- bsplbase(l./max(grid$l), bPars[1,  ])$knots
knots.l <- knots.l[1:(length(knots.l)-(bPars[1,4]+1))]

m. <- grid$m

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



B_lm <- t(sapply(1:nrow(Bl),function(i){as.vector(outer(Bm[i,],Bl[i,]))}))
basisKeepIndex <- knot_grid$keep[!duplicated(knot_grid[,c("l_index","m_index")])] %>% equals(TRUE) %>% which
B_lm <- B_lm[,basisKeepIndex]

U_lm <- X%*%B_lm







knot_grid <- subset(knot_grid,keep==TRUE)
dl <- bPars[1, 6]
if(dl>0){
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
if(dl==0){
      Pl <- diag(nrow(knot_grid))
}
lambdal <- bPars[1, 5]









dm <- bPars[2, 6]
if(dm>0){
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
}
if(dm==0){
      Pm <- diag(nrow(knot_grid))
}

lambdam <- bPars[2, 5]
















Pen <- rbind(lambdal*Pl, lambdam*Pm)#,lambdaridge*Pridge)
n.col <- ncol(U_lm)
nix <- rep(0,nrow(Pen))

nix.ridge <- rep(0, ncol(U_lm))
coef.est <- rep(1, ncol(U_lm))

mu <- rep(mean(y_vec), length(y_vec))
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
      #w <- rep(10, length(y_vec))
      #u <- (y[-1] - mu)/h.prime + eta
      u <- (y_vec - mu)/h.prime + eta
      
      startTS <- Sys.time()
      f <- lsfit(rbind(U_lm,Pen), c(u, nix), wt = c(w, nix + 1) *c(w, nix + 1), intercept = F)
      endTS <- Sys.time()
      endTS-startTS
      
      coef.old <- coef.est
      coef.est <- as.vector(f$coef)
      
      d.coef <- max(abs((coef.est[coef.old>0] - coef.old[coef.old>0])/coef.old[coef.old>0]))
      if(d.coef < 1e-008)
            break
      print(c(it, d.coef))
      eta <- U_lm %*% coef.est
}

if(it > 24) {
      warning(paste("parameter estimates did NOT converge in 25 iterations"
      ))
}

H <- hat(f$qr, intercept = F)[1:(m-1)]
trace <- eff.dim <- sum(H)

cv <- press.mu <- press.e <- var.c <- NULL
dev <- sum(f$residuals[1:(m-1)]^2)
dispersion.parm <- dev/((m-1) - trace)
press.e <- f$residuals[1:(m-1)]/(1 - H)
cv <- sqrt(sum((press.e)^2)/(m-1))
press.mu <- y_vec - press.e

aic <- dev + 2 * trace
aic








##-----------------------------------------------------------------------------


gg <- expand.grid(s=(1:50),t=(1:50)) %>% subset(.,t>s)
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









Phi_hat <- B_lm %*% f$coefficients
data.frame(phi=Phi_hat,gg) %>%
      wireframe(phi~s*t,
                data=.,
                screen=list(z=80,x=-75),
                light.source = c(5,20,10),
                pretty=TRUE,
                scales = list(arrows = FALSE),
                drape=FALSE,
                cex=0.15,
                col="grey",
                par.settings = list(axis.line = list(col = "transparent")))







true_phi <- 2*(gg$t/max(gg$t))^2-0.5
true_phi[gg$l >1] <- 0
data.frame(true_phi=true_phi,gg) %>%
      wireframe(true_phi~s*t,
                data=.,
                screen=list(z=80,x=-75),
                light.source = c(5,20,10),
                pretty=TRUE,
                scales = list(arrows = FALSE),
                drape=FALSE,
                cex=0.15,
                col="grey",
                par.settings = list(axis.line = list(col = "transparent")))


















###########################################################################################################
###########################################################################################################
###########################################################################################################









bPars <- rbind(c(0,1,70,3,100,2),
               c(0,1,70,3,100,3))



phi <- function(t1.index,t2.index,m){
      if(t1.index>t2.index){
            t_ij <- -exp(-2*(t1.index-t2.index)/((m-1)))             
      }
      if(t1.index==t2.index){
            t_ij <- 1             
      }
      if(t1.index<t2.index){
            t_ij <- 0             
      }
      t_ij
}




full_grid <- expand.grid(t=1:m,s=1:m) %>% orderBy(~t+s,.)
true_T <- sapply(1:nrow(full_grid),function(row.i){
      phi(full_grid$t[row.i],full_grid$s[row.i],m=m)
}) %>%unlist


T_mat <- matrix(data=true_T,nrow=m,ncol=m,byrow=TRUE)
T_mat


N <- 10
y <- t(solve(T_mat)%*%matrix(data=rnorm(N*m,mean=0,sd=1),
                             nrow=m,
                             ncol=N))
#iv_errors <- rep(0.1,m)
#D <- diag(iv_errors)
D <- diag(m) 

true_Sigma <- solve(T_mat)%*%D%*%t(solve(T_mat))
true_Omega <- t(T_mat)%*%solve(D)%*%T_mat


y_vec <- as.vector(t(y[,-1]))

X <- matrix(data=0,nrow=N*(M-1),ncol=choose(M,2))
no.skip <- 0
for (t in 2:M){
      X[((0:(N-1))*(M-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- y[,1:(t-1)]
      no.skip <- no.skip + t - 1
}













l. <- grid$l
Bl <- bsplbase(l./max(grid$l), bPars[1,  ])$base
knots.l <- bsplbase(l./max(grid$l), bPars[1,  ])$knots
knots.l <- knots.l[1:(length(knots.l)-(bPars[1,4]+1))]

m. <- grid$m

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



B_lm <- t(sapply(1:nrow(Bl),function(i){as.vector(outer(Bm[i,],Bl[i,]))}))
basisKeepIndex <- knot_grid$keep[!duplicated(knot_grid[,c("l_index","m_index")])] %>% equals(TRUE) %>% which
B_lm <- B_lm[,basisKeepIndex]

U_lm <- X%*%B_lm







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
dm <- 1
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

lambdas_2 <- expand.grid(lam_l=exp(seq(-3,3.33,length.out=10)),
                       lam_m=exp(seq(-3,3.33,length.out=10)))

lambdas_3 <- expand.grid(lam_l=exp(seq(-3,3.33,length.out=10)),
                         lam_m=exp(seq(-2,5,length.out=10)))

coef_list_3 <- list.zip(lam_l=lambdas_3$lam_m,
            lam_m=lambdas_3$lam_l) %>%
      lapply(.,function(l){
            fit_cholesky_PS(y_vec,U_lm,Pl,l$lam_l,
                            Pm,l$lam_m,
                            0.12)    
      })



setwd("/Users/taylerblake/Documents/Dissertation/code/simulations")
save(coef_list,file="coef_list.Rdata")
save(coef_list_2,file="coef_list_2.Rdata")
save(coef_list_2,file="coef_list_3.Rdata")
##-----------------------------------------------------------------------------


gg <- expand.grid(s=(1:30),t=(1:30)) %>% subset(.,t>s)
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
phi_list_3 <- lapply(coef_list_3,function(beta){
      B_lm%*%beta
})
i <- 1


#Phi_hat <- B_lm %*% f$coefficients

laml <- round(lambdas_3$lam_l[i],3)
lamm <- round(lambdas_3$lam_m[i],3)
my_title <- c(as.expression(bquote(lambda[l] == .(laml))), as.expression(bquote(lambda[m] == .(lamm))))
Phi_hat <- phi_list_2[[i]]
data.frame(phi=Phi_hat,gg) %>%
      subset(.,l<(max(gg$l)-2)) %>%
      wireframe(phi~s*t,
                data=.,
                screen=list(z=65,x=-85),
                light.source = c(5,20,10),
                pretty=TRUE,
                scales = list(arrows = FALSE),
                drape=FALSE,
                cex=0.15,
                col="grey",
                par.settings = list(axis.line = list(col = "transparent")),
                main=my_title
                )

i <- i+1



true_phi <- apply(gg,MARGIN = 1,FUN=function(row){
      -phi(row[2],row[1],m=m)
})

data.frame(true_phi=true_phi,gg) %>%
      wireframe(true_phi~s*t,
                data=.,
                screen=list(z=65,x=-85),
                light.source = c(5,20,10),
                pretty=TRUE,
                scales = list(arrows = FALSE),
                drape=FALSE,
                cex=0.15,
                col="grey",
                par.settings = list(axis.line = list(col = "transparent")))





lambdas <- expand.grid(lam_l=exp(seq(3.33,4.5,length.out=6)),
                       lam_m=exp(seq(3.33,4.5,length.out=6)))


coef_list_2 <- list.zip(lam_l=lambdas$lam_m,
                      lam_m=lambdas$lam_l) %>%
      lapply(.,function(l){
            fit_cholesky_PS(y_vec,U_lm,Pl,l$lam_l,
                            Pm,l$lam_m,
                            0.12)    
      })



phi_list_2 <- lapply(coef_list_2,function(beta){
      B_lm%*%beta
})


list.append(phi_list,phi_list_2) %>% length





































library(inline)
library(Rcpp)
library(RcppArmadillo)
rcpp_inc <- 'using namespace Rcpp;
using namespace arma;
'

## Rcpp code to calculate the matrix pseudo inverse 
src <- ' mat m1 = as<mat>(m1in);
//return(wrap(m1.i()));

mat m2 = inv(m1);
return(wrap(m2));
'
fn <- cxxfunction(signature(m1in="numeric"), src, plugin='RcppArmadillo', rcpp_inc)

m3 <- matrix(runif(16), nr=4)
res <- fn(m3)
test <- solve(m3)
all.equal(test, res)



###########################################################################