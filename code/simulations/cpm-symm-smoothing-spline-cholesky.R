
library(doParallel)
library(foreach)
library(doBy)
library(lattice)
library(MASS)
library(magrittr)
library(rlist)
library(plyr)
library(stringr)
library(dplyr)
library(rlist)

M <- m <- 50
N <- 30

grid <- expand.grid(s=(1:m),t=(1:m)) %>%
  subset(.,s>t) %>%
  transform(l=(s-t),m=(s+t)/2)

Sigma <- matrix(0.7,nrow=m,ncol=m) + diag(rep(0.3),m)
Omega <- solve(Sigma)

C <- t(chol(Sigma))
D <- diag(diag(C))
L <- C%*%solve(D)
T_mod <- solve(L)


y <- t(solve(T_mat)%*%matrix(data=rnorm(N*m,mean=0,sd=0.3),
                             nrow=m,
                             ncol=N))
true_Sigma <- solve(T_mat)%*%t(solve(T_mat))
true_Omega <- t(T_mat)%*%T_mat

y <- mvrnorm(n=N,mu=rep(0,m),Sigma=Sigma)
y_vec <- as.vector(t(y[,-1]))

X <- matrix(data=0,nrow=N*(M-1),ncol=choose(M,2))
no.skip <- 0
for (t in 2:M){
  X[((0:(N-1))*(M-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- y[,1:(t-1)]
  no.skip <- no.skip + t - 1
}

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
## define basis functions, representers
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir(file.path(getwd(),"lib"))


k0 <- function(x){
      return(rep(1,length(x)))
}		
k1 <- function(x){
      return(x-0.5)
}
k2 <- function(x){
      return( 0.5*((k1(x))^2 - (1/(12))) )
}
k4 <- function(x){
      return(  (1/24)*( (k1(x))^4 -((k1(x))^2/2) + (7/240)) )
}




R1 <- function(l1,l2,m){
      
      if(m==1){
            representer <- k1(l1)*k1(l2) +  k2(l1)*k2(l2) - k4( abs(l1-l2)  )
      }
      if(m==2){
            representer <- k2(l1)*k2(l2) - k4( abs(l1-l2) ) 
      }
      representer
}





#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
## construct B, K
#-----------------------------------------------------------------------------------------



############################################################
## TODO: MAKE SURE YOURE CONSTRUCTING K CORRECTLY - 
##        SPECIFICALLY THE INTERACTION TERM R1_LM

R1_l <- function(l1,l2){R1(l1,l2,m=1)}
R1_L <- outer(grid$l/max(grid$l), grid$l/max(grid$l), "R1_l")

R1_m <- function(m1,m2){R1(m1,m2,m=1)}
R1_M <- outer(grid$m/max(grid$m), grid$m/max(grid$m), "R1_m")

R1_LM <- R1_L*R1_M + outer(k1(grid$l/max(grid$l)),k1(grid$l/max(grid$l)))*R1_M
K <- R1_L + R1_M + R1_LM

# B <- matrix(data=c(rep(1,nrow(grid)),k1(grid$l/max(grid$l))*k1(grid$m/max(grid$m))),
#             nrow=nrow(grid),ncol=2,byrow=FALSE)
B <- matrix(data=c(rep(1,nrow(grid))),ncol=1,nrow=nrow(grid),byrow=FALSE)
QR_B <- qr(B,complete=TRUE)
Q_B <- qr.Q(QR_B,complete=TRUE)
Q2_B <- Q_B[,(ncol(B)+1):ncol(Q_B)]
Q1_B <- Q_B[,1:ncol(B)]
R_B.big <- qr.R(QR_B,complete=TRUE)
R_B <- R_B.big[1:ncol(B),]
R_Binv <- solve(R_B)

Dinv <- diag(rep(1,length(y_vec)))







QR_X <- qr(X,complete=TRUE)
Q_X <- qr.Q(QR_X,complete=TRUE)
Q2_X <- Q_X[,(ncol(B)+1):ncol(Q_X)]
Q1_X <- Q_X[,1:ncol(X)]
R_X.big <- qr.R(QR_X,complete=TRUE)
R_X <- R_X.big[1:ncol(X),]
R_Xinv <- solve(R_X)

#-----------------------------------------------------------------------------------------
## Build solutions
#-----------------------------------------------------------------------------------------

lambdas <- as.list(exp(seq(-8,8,length.out=100)))
P <- solve(t(X)%*%Dinv%*%X)

Ms <- lapply(lambdas,function(l){
      M <- solve( t(Q2_B)%*%(K + l*P)%*%Q2_B )
      M
})
      

c <- lapply(Ms,function(mat){
      Q2_B%*%mat%*%t(Q2_B)%*% P %*%t(X)%*%Dinv%*%y_vec      
})

d <- lapply(list.zip(lam=lambdas,c=c),function(l){
      d <- R_Binv%*%t(Q1_B)%*%( P%*%t(X)%*%Dinv%*%y_vec - ( K + l$lam*P )%*%l$c )      
})
cholesky <- lapply(list.zip(c=c,d=d),function(l){
      Phi <- B%*%l$d + K%*%l$c
      T_hat <- diag(rep(1,m))
      T_hat[lower.tri(T_hat)] <- -Phi
      list(phi=Phi,T_mat=T_hat,omega=t(T_hat)%*%T_hat)
})


entropy_loss <- function(trueSigma, omegaHat){
      I_hat <- trueSigma%*%omegaHat
      sum(diag(I_hat)) -log(det(I_hat)) - ncol(omegaHat)
}

lapply(cholesky,function(l){
      entropy_loss(true_Sigma,l$omega)
}) %>% unlist %>% plot(x=log(unlist(lambdas)),
                       y=.,
                       type="l",
                       ylab=expression(Delta[1]),
                       xlab=expression(log(lambda)))



quadratic_loss <- function(trueSigma, omegaHat){
      I_hat <- trueSigma%*%omegaHat
      sum( diag(I_hat-diag(1,ncol(omegaHat)))^2  )
}


lapply(cholesky,function(l){
      quadratic_loss(true_Sigma,l$omega)
}) %>% unlist %>% plot(x=log(unlist(lambdas)),
                       y=.,
                       type="l",
                       ylab=expression(Delta[2]),
                       xlab=expression(log(lambda)))


















Rl_gg <- sapply(grid$l/max(grid$l),function(grid_l) {sapply(seq(0,1,length.out=200),function(pred_l){R1_l(pred_l,grid_l)})})
Bl_gg <- matrix(data=c( rep(1,200), k1(seq(0,1,length.out=200))),nrow=200,ncol=2,byrow=FALSE)

l_smooth <- list.zip(c=c,d=d) %>% lapply(.,function(l){
      Rl_gg%*%l$c #+ Bl_gg%*%l$d
      })


Rm_gg <- sapply(grid$m/max(grid$m),function(grid_m) {sapply(seq(0,1,length.out=200),function(pred_m){R1_m(pred_m,grid_m)})})
Bm_gg <- matrix(data=c( rep(1,200), k1(seq(0,1,length.out=200))),nrow=200,ncol=2,byrow=FALSE)
m_smooth <- list.zip(c=c,d=d) %>% lapply(.,function(l){
      Rm_gg%*%l$c #+Bm_gg%*%l$d
})

l_smooth <- list.cbind(l_smooth) 
m_smooth <- list.cbind(m_smooth) 

matplot(seq(0,1,length.out=200),l_smooth,
        col=terrain.colors(100,alpha=0.7),type="l",
        xlab="l",
       ylab= expression(phi[l]))

matplot(seq(0,1,length.out=200),m_smooth,
        col=terrain.colors(100,alpha=0.7),type="l",
        xlab="m",
        ylab= expression(phi[m]))


gg <- expand.grid(l=seq(0,1,length.out=200),
                  m=seq(0,1,length.out=200))

Rl_gg <- sapply(grid$l/max(grid$l),
                function(grid_l){
                      sapply(gg$l,
                             function(pred_l){
                                   R1_l(pred_l,grid_l)
                                   })})

Rm_gg <- sapply(grid$m/max(grid$m),
                function(grid_m){
                      sapply(gg$m,
                             function(pred_m){
                                   R1_m(pred_m,grid_m)
                             })})

lm_smooth <- lapply(c,function(coef){
      as.vector((Rl_gg*Rm_gg)%*%coef)    
})









jet.colors <- colorRampPalette( c("deepskyblue2","green") )
nbcol <- 100
color <- jet.colors(nbcol)


lm_smooth[[lapply(cholesky,function(l){quadratic_loss(true_Sigma,l$omega)}) %>% unlist %>%which.min]] %>%
      data.frame(gg,phi_lm=.) %>% wireframe(phi_lm~l+m,
                                            data=.,
                                            screen=list(z=20,x=-75),
                                            light.source = c(5,20,10),
                                            col="grey",
                                            scales = list(arrows = FALSE),
                                            drape=FALSE,
                                            cex=0.15,
                                            colorkey=FALSE,
                                            par.settings = list(axis.line = list(col = "transparent")))
lm_smooth <- list.cbind(lm_smooth) %>% as.vector

library(ggplot2)
library(ggthemes)

best_phi_lm <- lm_smooth[[lapply(cholesky,function(l){quadratic_loss(true_Sigma,l$omega)}) %>% unlist %>%which.min]]
best_phi_m <- m_smooth[[lapply(cholesky,function(l){quadratic_loss(true_Sigma,l$omega)}) %>% unlist %>%which.min]]
true_phi_m <- data.frame(m=seq(0,1,length.out=200),
                  phi=2*(seq(0,1,length.out=200)^2 + seq(0,1,length.out=200) ) )


library(doBy)
data.frame(lambda=expand.grid(1:40000,lam=unlist(lambdas))$lam,
           phi=lm_smooth,
           l=rep(gg$l,length(lambdas)),
           m=rep(gg$m,length(lambdas))) %>%
      ggplot(.,aes(x=m,y=phi,group=lambda)) + geom_line(aes(colour=lambda)) +
      scale_color_continuous_tableau(palette = "Green") +
      theme_minimal() 







best_phi_lm <- lm_smooth[[lapply(cholesky,function(l){quadratic_loss(true_Sigma,l$omega)}) %>% unlist %>%which.min]]
best_phi_m <- m_smooth[[lapply(cholesky,function(l){quadratic_loss(true_Sigma,l$omega)}) %>% unlist %>%which.min]]
true_phi_m <- data.frame(m=seq(0,1,length.out=200),
                         phi=2*(seq(0,1,length.out=200)^2 + seq(0,1,length.out=200) ) )


data.frame(phi=best_phi_m+best_phi_lm,gg) %>%
      subset(.,l%in%seq(0,1,length.out=200)[c(1:5,seq(25,200,by=25))]) %>%
      ggplot(.,aes(x=m,y=phi,group=l)) + geom_line(aes(colour=l)) +
      scale_color_continuous_tableau(palette = "Green") +
      theme_minimal() +
      geom_line(data=true_phi_m,
                aes(x=m,y=phi),
                inherit.aes = FALSE,
                colour="red") 









