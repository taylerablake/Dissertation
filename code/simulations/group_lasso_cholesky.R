



library(plyr)
library(dplyr)
library(doBy)
library(tidyr)
library(splines)
library(magrittr)
library(ggplot2)
library(grpreg)








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


N <- 20
y <- t(solve(T_mat)%*%matrix(data=rnorm(N*m,mean=0,sd=1),
                             nrow=m,
                             ncol=N))
true_Sigma <- solve(T_mat)%*%t(solve(T_mat))
true_Omega <- t(T_mat)%*%T_mat


y_vec <- as.vector(t(y[,-1]))

X <- matrix(data=0,nrow=N*(M-1),ncol=choose(M,2))
no.skip <- 0
for (t in 2:M){
      X[((0:(N-1))*(M-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- y[,1:(t-1)]
      no.skip <- no.skip + t - 1
}

true_Sigma <- solve(T_mat)%*%t(solve(T_mat))
true_Omega <- t(T_mat)%*%T_mat




####################################################################################################
####################################################################################################

## build B-spline marginal bases for l,m, tensor product basis for l x m interaction



l. <- grid$l
Bl <- bsplbase(l./max(l.), bPars[1,  ])$base
knots.l <- bsplbase(l./max(l.), bPars[1,  ])$knots
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

B <- cbind(Bl[,colSums(abs(Bl)>0)],Bm[,colSums(abs(Bm)>0)],B_lm[,colSums(abs(B_lm)>0)])
active_knot_indices <- c(which(colSums(abs(Bl)>0)>0),
                         ncol(Bl)+which(colSums(abs(Bm)>0)>0),
                         ncol(Bl)+ncol(Bm)+which(colSums(abs(B_lm)>0)>0))
U <- X%*%B


 coef_group <- rep(1,ncol(B))
# coef_group[unique(which(abs(Bl[grid$l %in% sort(unique(grid$l))[1:3],colSums(abs(Bl)>0)]) > 0,arr.ind=TRUE)[,2])] <- 2
 #coef_group[sum(colSums(abs(Bl)>0)) + sum(colSums(abs(Bm)>0)) + unique(which(abs(B_lm[grid$l%in%sort(unique(grid$l))[1:3],]) > 0,arr.ind=TRUE)[,2])] <- 3
 coef_group[sum(colSums(abs(Bl)>0)) + ] <- 2
 coef_group[(sum(colSums(abs(Bl)>0))+1):(sum(colSums(abs(Bl)>0))+sum(colSums(abs(Bm)>0)))] <- 4
l_penalty_w <- (Bl[,colSums(abs(Bl)>0)>0]>0) %>% apply(.,MARGIN=2,FUN=which,arr.ind=TRUE) %>% lapply(.,function(x){
      max(grid$l[x])^2
}) %>% unlist


group_multiplier <- log(c(l_penalty_w,
                      rep(max(l_penalty_w),sum(colSums(abs(Bm)>0)>0)),
                      l_penalty_w[match(knot_grid$l[which(colSums(abs(B_lm)>0)>0)],knots.l)])) +1



fit_bridge <- gBridge(U,y_vec,alpha=0.8,
                      group=c(rep(1,sum(colSums(Bl>0)>0)),
                              rep(2,sum(colSums(Bm>0)>0)),
                              rep(3,sum(colSums(B_lm>0)>0))),
                      lambda.min = 0.0001)
fit_bridge <- gBridge(U,y_vec,alpha=0.7,lambda.min = 0.0001)
plot(fit_bridge)


best_bridge <- gBridge(U,y_vec,
                       group=c(rep(1,sum(colSums(Bl>0)>0)),
                               rep(2,sum(colSums(Bm>0)>0)),
                               rep(3,sum(colSums(B_lm>0)>0))),
                       alpha=0.8,
                       lambda=select(fit_bridge,criterion = "GCV")$lambda)
best_bridge <- gBridge(U,y_vec,
                       alpha=0.7,
                       lambda=select(fit_bridge,criterion = "GCV")$lambda)

y_hat <- cbind(rep(1,nrow(U)),U)%*%coefficients(best_bridge)
sum(abs(coefficients(best_bridge))>0)
plot((y_vec-y_hat)*(1/y_vec))




plot(best.fit)
      
coefficients(best.fit)[abs(coefficients(best.fit))>0]
(y_vec-cv.fit$Y[,match(cv.fit$lambda.min,cv.fit$lambda)])^2


Phi_hat <- cbind(rep(1,nrow(B)),B)%*%coefficients(best.fit)















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
B_lm <- B_lm[,basisKeepIndex]

B_gg <- cbind(Bl,Bm,B_lm)
B_gg <- B_gg[,active_knot_indices]

Phi_hat_gg <- cbind(rep(1,nrow(B_gg)),B_gg)%*%coefficients(best.fit)
data.frame(phi=Phi_hat_gg,gg) %>%
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




