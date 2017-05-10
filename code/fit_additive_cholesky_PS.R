


## for varying penalties on the components corresponding
## to l and m (and the interaction), construct
##          - the corresponding null basis and 
##            representers of the penalized space
##          - cholesky factors T which belong to the
##            null space of J

## verify that as lambda -> infinity, the solution approaches
## the correct functional form (and also the true cholesky factor!)
for (helper.file in list.files(file.path(getwd(),"lib"))) {
  source(file.path(getwd(),"lib",helper.file))
}


cl <- makeCluster(detectCores()-3)
registerDoParallel(cl)
clusterCall(cl,function() {
  .libPaths("~/Rlibs/lib")
  library(doBy)
  library(lattice)
  library(MASS)
  library(magrittr)
  library(rlist)
  library(plyr)
  library(stringr)
  library(dplyr)
  library(doParallel)
})


M <- 20
N <- 100
grid <- build_grid(20)



Sigma <- matrix(0.7,nrow=m,ncol=m) + diag(rep(0.3),m)
Omega <- solve(Sigma)


C <- t(chol(Sigma))
D <- diag(diag(C))
L <- C%*%solve(D)
T_mod <- solve(L)
D <- diag(rep(diag(C)[-1],N))


bPars <- rbind(c(0,1,30,3,100,3),
               c(0,1,30,3,100,3))

Bl <- bsplbase(grid$l/max(grid$l),
               bPars[1,],outer.okay = TRUE)$base
Bm <- bsplbase(grid$m/max(grid$m),
               bPars[2,],outer.okay = TRUE)$base

n1 <- ncol(Bl)
n2 <- ncol(Bm)



y <- mvrnorm(n=N,mu=rep(0,M),Sigma=Sigma)
y_vec <- as.vector(t(y[,-1]))


X <- matrix(data=0,nrow=N*(M-1),ncol=choose(M,2))
no.skip <- 0
for (t in 2:M){
  X[((0:(N-1))*(M-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- y[,1:(t-1)]
  no.skip <- no.skip + t - 1
}

B <- cbind(as.vector(rep(1,nrow(Bl))),Bl,Bm)
U <- X %*% B









dl=3
Dl <- diff(diag(ncol(Bl)),
           differences=dl)

dm=3
Dm <- diff(diag(ncol(Bm)),
           differences=dm)

nl <- ncol(Dl)
nm <- ncol(Dm)

lambdas <- expand.grid(laml=exp(seq(-0.5,7,length.out=10)),
                       lamm=exp(seq(-0.5,7,length.out=10)))
ridge.adj <- 0.0000001


  




cl <- makeCluster(detectCores()-3)
registerDoParallel(cl)
clusterCall(cl,function() {
  .libPaths("~/Rlibs/lib")
  library(doBy)
  library(lattice)
  library(MASS)
  library(magrittr)
  library(rlist)
  library(plyr)
  library(stringr)
  library(dplyr)
  library(doParallel)
})



clusterExport(cl,c("grid",
                   "bsplbase",
                   "fit_cholesky_PS",
                   "Sigma",
                   "N",
                   "M",
                   "y",
                   "y_vec",
                   "y_aug",
                   "Bl",
                   "Bm",
                   "B.",
                   "Pl",
                   "Pm",
                   "lambdas",
                   "dl",
                   "dm",
                   "U.",
                   "Pl",
                   "Pm"))

fit_list <- foreach(lambda_pair=iter(lambdas,by="row")) %dopar% {
  
  Pen <- rbind(cbind(rep(0,nrow(Dl)),lambda_pair$laml*Dl,
                     matrix(data=0,nrow=nrow(Dl),ncol=ncol(Dm))),
               cbind(rep(0,nrow(Dl)),
                     matrix(data=0,nrow=nrow(Dm),ncol=ncol(Dl)),
                     lambda_pair$lamm*Dm))
  
  
    Pen <- rbind(Pen,
                 ridge.adj*cbind(rep(0,nl+nm),diag(nl+nm)))

  nix <- as.vector(rep(0,nrow(Pen)))
  w_aug <- as.vector(c(1/diag(D),rep(1,nrow(Pen))))
  y_aug <-  as.vector(c(as.vector(t(y[,-1])),rep(0,nrow(Pen)))) 
  
  fit <- lsfit(rbind(U,Pen), y_aug,
               wt = w_aug, intercept = F)
  fit
}



phi_list <- lapply(fit_list,function(l) {
  phi <- diag(M)
  phi[lower.tri(phi)] <- B %*% l$coefficients
  phi
})


i <- 1


laml <- round(lambdas$laml[i],3)
lamm <- round(lambdas$lamm[i],3)

my_title <- c(as.expression(bquote(lambda[l] == .(laml))),
              as.expression(bquote(lambda[m] == .(lamm))))
wireframe(phi_list[[i]],
          scales = list(arrows = FALSE),
          xlab="",
          ylab="",
          zlab=expression(phi),
          main=my_title)
i <- i+1























B <- cbind(as.vector(rep(1,nrow(Bl))),
           Bl,Bm, 
           grid$m/max(grid$m)*Bl,
           #((grid$m/max(grid$m))^2)*Bl,
           grid$l/max(grid$l)*Bm)
           #((grid$l/max(grid$l))^2)*Bm)
U <- X %*% B


clusterExport(cl,c("grid",
                   "bsplbase",
                   "fit_cholesky_PS",
                   "Sigma",
                   "N",
                   "M",
                   "y",
                   "y_vec",
                   "y_aug",
                   "Bl",
                   "Bm",
                   "B.",
                   "Pl",
                   "Pm",
                   "lambdas",
                   "dl",
                   "dm",
                   "U.",
                   "Pl",
                   "Pm"))

lambdas <- expand.grid(lam_l1=exp(seq(-2,2,length.out=10)),
                       lam_m1=exp(seq(-2,2,length.out=10)),
                       lam_l2=exp(seq(-2,2,length.out=10)),
                       lam_m2=exp(seq(-2,2,length.out=10)))
fit_list <- foreach(lambda_pair=iter(lambdas,by="row")) %dopar% {
  
  Pen <- rbind(cbind(rep(0,nrow(Dl)),
                     lambda_pair$lam_l1*Dl,
                     matrix(data=0,nrow=nrow(Dl),
                            ncol=((2*ncol(Dm))+ncol(Dl)))),
               ######################################################
               cbind(rep(0,nrow(Dl)),
                     matrix(data=0,nrow=nrow(Dm),
                            ncol=ncol(Dl)),
                     lambda_pair$lam_m1*Dm,
                     matrix(data=0,nrow=nrow(Dm),ncol=ncol(Dl)+ncol(Dm))),
               ######################################################
               cbind(rep(0,nrow(Dl)),
                     matrix(data=0,nrow=nrow(Dm),
                            ncol=ncol(Dl)+ncol(Dm)),
                     lambda_pair$lam_l2*Dl,
                     matrix(data=0,nrow=nrow(Dm),ncol=ncol(Dm))),
               ######################################################
               cbind(rep(0,nrow(Dl)),
                     matrix(data=0,nrow=nrow(Dm),
                            ncol=((2*ncol(Dl))+ncol(Dm))),
                     lambda_pair$lam_m2*Dm))
  
  
  Pen <- rbind(Pen,
               ridge.adj*cbind(rep(0,2*(nl+nm)),diag(2*(nl+nm))))
  
  nix <- as.vector(rep(0,nrow(Pen)))
  w_aug <- as.vector(c(1/diag(D),rep(1,nrow(Pen))))
  y_aug <-  as.vector(c(as.vector(t(y[,-1])),rep(0,nrow(Pen)))) 
  
  fit <- lsfit(rbind(U,Pen), y_aug,
               wt = w_aug, intercept = F)
  fit
}






phi_list <- lapply(fit_list,function(l) {
  phi <- diag(M)
  phi[lower.tri(phi)] <- B %*% l$coefficients
  phi
})
i <- 1


laml1 <- round(lambdas$lam_l1[i],2)
lamm1 <- round(lambdas$lam_m1[i],2)
laml2 <- round(lambdas$lam_l2[i],2)
lamm2 <- round(lambdas$lam_m2[i],2)
my_title <- c(as.expression(bquote(lambda[l1] == .(laml1))),
              as.expression(bquote(lambda[m1] == .(lamm1))),
              as.expression(bquote(lambda[l2] == .(laml2))),
              as.expression(bquote(lambda[m2] == .(lamm2))))
wireframe(phi_list[[i]],
          scales = list(arrows = FALSE),
          xlab="",
          ylab="",
          bty="n",
          zlab=expression(phi),
          main=my_title)
i <- i+1


