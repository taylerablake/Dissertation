
library(splines)
library(magrittr)
library(ggplot2)
library(plyr)
library(dplyr)

































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


N <- 1
y <- solve(T_mat)%*%matrix(data=rnorm(N*m,mean=0,sd=1),nrow=m,ncol=N)
true_Sigma <- solve(T_mat)%*%t(solve(T_mat))
true_Omega <- t(T_mat)%*%T_mat






l. <- as.vector(outer(rep(1,length(unique(grid$m))),unique(grid$l)))
Bl <- bsplbase(l./max(grid$l), bPars[1,  ])$base
knots.l <- bsplbase(l./max(grid$l), bPars[1,  ])$knots
knots.l <- knots.l[1:(length(knots.l)-(bPars[1,4]+1))]

m. <- as.vector(outer(unique(grid$m),rep(1,length(unique(grid$l)))))
Bm <- bsplbase(m./max(grid$m), bPars[2,  ])$base
knots.m <- bsplbase(m./max(grid$m), bPars[2,  ])$knots
knots.m <- knots.m[1:(length(knots.m)-(bPars[2,4]+1))]

knotGrid <- expand.grid(m=knots.m,l=knots.l)[,2:1]
keep <- (knotGrid$m >= 0.5) & (knotGrid$m < max(knotGrid$l)-0.5*knotGrid$l) |
      (knotGrid$m < 0.5) & (knotGrid$m > min(knotGrid$l)+0.5*knotGrid$l)
knotGrid <- transform(knotGrid,keep=factor(keep))
knotGrid <- data.frame(knotGrid,
                       expand.grid(m_index=1:length(knots.m),l_index=(1:length(knots.l)))[,2:1])






ggplot(knotGrid,aes(l,m)) + geom_point(aes(colour=keep))

knotGrid <- knotGrid %>% orderBy(~l_index + m_index,.)

n1 <- ncol(Bl)
n2 <- ncol(Bm)  # Compute tensor products

Bl. <- kronecker(Bl, t(rep(1, n2)))
Bm. <- kronecker(t(rep(1, n1)), Bm)




# grid <- orderBy(~ l+m,grid)
# discard_rows <- left_join(data.frame(l=l.,m=m.,frame="big"),
#                           data.frame(grid[,c("l","m")],frame="little")
#                           ,by=c("l","m"))$frame.y %>% is.na %>% which
# 
# B. <- Bl.*Bm.
# #basisKeepIndex <- knotGrid$keep[!duplicated(knotGrid[,c("l_index","m_index")])] %>% equals(TRUE) %>% which
# #B. <- B.[-discard_rows,basisKeepIndex]
# B. <- B.[-discard_rows,]
# 
# 
# 
# 
# ##
## TO DO: for the knots that lie outside the support of Phi, allocate a weight of 0 to these
## observations to use the following lsfit approach

# 
# 
# 
# 
# 
# #dl <- bPars[1, 6]
# dl <- 1
# Dl <- diag(n1)
# if(dl != 0) {
#       for(j in 1:dl) {
#             Dl <- diff(Dl)
#       }
# }
# lambdal <- bPars[1, 5]
# Pl <-   sqrt(lambdal)*kronecker(Dl, diag(n2))
# 
# 
# dm <- bPars[2, 6]
# #d2 <- 1
# Dm <- diag(n2)
# if(dm != 0) {
#       for(j in 1:dm) {
#             Dm <- diff(Dm)
#       }
# }
# lambdam <- bPars[2, 5]
# Pm <- sqrt(lambdam) * kronecker(diag(n1), Dm)
# Pen <- rbind(Pl, Pm)
# 
# 
# 
# 
# 
# 
# reference_grid <- data.frame(l=l.,m=m.)[-discard_rows,]
# match_row <- function(row,df){
#       match.ind <- match(df[,1],row[1])
#       for(j in 2:ncol(df)){
#             match.ind <- match.ind*match(df[,j],row[j])
#       }
#       which(!is.na(match.ind))
# }
# x_col_ind <-   sapply(1:nrow(grid),function(row.i){
#       match_row(grid[row.i,c("l","m")],reference_grid)
# })
# 
# grid$column_index <- x_col_ind
# 
# X <- matrix(data=0,nrow=m-1,ncol=nrow(reference_grid))
# for(row.i in 1:nrow(grid)){
#       X[grid$t[row.i]-1,grid$column_index[row.i]] <- as.vector(y)[grid$s[row.i]] 
# }
# U. <- X%*%B.
# 
# 
# 
# 
# z1 <- rep(0, n2 * (n1 - dl))
# z2 <- rep(0, n1 * (n2 - dm))
# 
# n.col <- ncol(U.)
# nix <- c(z1, z2)
# 
# nix.ridge <- rep(0, n1 * n2)
# coef.est <- rep(1, ncol(U.))
# 
# mu <- rep(mean(y[-1]), length(y[-1]))
# it <- 0
# repeat {
#       if(it == 0) {
#             eta <- mu
#       }
#       
#       it <- it + 1
#       if(it > 25)
#             break
#       
#       mu <- eta
#       h.prime <- 1
#       w <- rep(1, length(y[-1]))
#       u <- (y[-1] - mu)/h.prime + eta
#       
#       
#       startTS <- Sys.time() 
#       f <- lsfit(rbind(U.,Pen), c(u, nix), wt = c(w, nix + 1) *c(w, nix + 1), intercept = F)
#       endTS <- Sys.time()
#       endTS-startTS
#       
#       startTS <- Sys.time() 
#       fFast <- fastLm(X=rbind(U.,Pen),y=c(u, nix))
#       endTS <- Sys.time()
#       endTS-startTS
#       
#       
#       
#       coef.old <- coef.est
#       coef.est <- as.vector(f$coef)
#       hist(as.vector(f$coef))
#       d.coef <- max(abs((coef.est - coef.old)/coef.old))
#       if(d.coef < 1e-008)
#             break
#       print(c(it, d.coef))
#       eta <- U. %*% coef.est
# }
# if(it > 24) {
#       warning(paste("parameter estimates did NOT converge in 25 iterations"
#       ))
# }
# 
# 
# H <- hat(f$qr, intercept = F)[1:(m-1)]
# trace <- eff.dim <- sum(H)
# 
# cv <- press.mu <- press.e <- var.c <- NULL
# dev <- sum(f$residuals[1:(m-1)]^2)
# dispersion.parm <- dev/((m-1) - trace)
# press.e <- f$residuals[1:(m-1)]/(1 - H)
# cv <- sqrt(sum((press.e)^2)/(m-1))
# press.mu <- y[-1] - press.e
# 
# aic <- dev + 2 * trace
# w.aug <- c(w, (c(z1, z2) + 1))
# 
# Phi_hat <- B. %*% f$coefficients
# pl <- length(unique(grid$l))
# pm <- length(unique(grid$m))
# 
# 
# Phi_mat <- matrix(Phi_hat, pl, pm, byrow = T)
# Phi_mat[!lower.tri(Phi_mat)] <- 0
# 
# i.1 <- 1:length(unique(grid$l))/max(grid$l)
# in1 <- 1:length(unique(grid$l))
# i.2 <- 1:length(unique(grid$m))/max(grid$m)
# in2 <- 1:length(unique(grid$m))
# if(length(pl) > 100) {
#       i.1 <- round(seq(from = min(grid$l)/max(grid$l), to = 1, length = 100))
#       in1 <- round(seq(from = 1, to = pl, length = 50))
# }
# if(length(M2.index) > 100) {
#       i.2 <- round(seq(from = min(grid$m)/max(grid$m), to = 1, length = 100))
#       in2 <- round(seq(from = 1, to = pm, length = 50))
# }
# 
# persp(sort((unique(grid$m)/max(grid$m))[in2]), 
#       sort((unique(grid$l)/max(grid$l))[in1]),
#       t(Phi_mat[in1, in2]),
#       theta=35,
#       phi=25,
#       xlab = "y.lab",
#       ylab = "x.lab", zlab = "phi")
# 
# 
# 
# grid <- orderBy(~l + m,grid)
# cloud(true_phi~ l + m,data=grid)
# 
# persp(grid$m,grid$l,grid$true_phi,theta=35,phi=25)





###########################################################################################################








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