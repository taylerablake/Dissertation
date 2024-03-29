
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

setwd(file.path("/Users","taylerblake","Documents","Dissertation","code"))
source(normalizePath(file.path(getwd(),"fnc/bsplbase.R")))
source(normalizePath(file.path(getwd(),"fnc/fit_cholesky_PS.R")))

cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
clusterCall(cl,function() {
      .libPaths("~/Rlibs/lib")
      library(doBy)
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

N <- 30
M <- m <- 30
grid <- expand.grid(s=(1:m),t=(1:m)) %>%
      subset(.,s>t) %>%
      transform(l=(s-t),m=s+t)



theta <- 0.8
sig2 <- 0.1



Sigma <- diag( rep( sig2*(1+theta^2), M ) )
for (i in 1:nrow(Sigma)) {
      for (j in 1:ncol(Sigma)) {
            if (abs(i-j)==1)
                  Sigma[i,j] <- sig2*theta
      }
}

Omega <- solve(Sigma)

C <- t(chol(Sigma))
D <- diag(diag(C))
L <- C%*%solve(D)
T_mod <- solve(L)      

grid <- expand.grid(t=1:m,s=1:m) %>%
      subset(.,t>s) %>%
      transform(.,l=t-s,
                m=(t+s)/2) %>%
      orderBy(~ l+m,.)

bPars <- rbind(c(0,1,length(unique(grid$l))-5,3,100,3),
               c(0,1,length(unique(grid$m))-5,3,100,3))


Bl <- bsplbase(as.vector(grid$l)/max(grid$l), bPars[1,])$base
Bm <- bsplbase(as.vector(grid$m)/max(grid$m), bPars[2,])$base

n1 <- ncol(Bl)
n2 <- ncol(Bm)

knots.l <- bsplbase(as.vector(grid$l)/max(grid$l), bPars[1,])$knots
interior.knots.l <- knots.l[1:(length(knots.l)-(bPars[1,4]+1))]

knots.m <- bsplbase(as.vector(grid$m)/max(grid$m), bPars[2,])$knots
interior.knots.m <- knots.m[1:(length(knots.m)-(bPars[1,4]+1))]

B. <- kronecker(Bm,
                t(as.vector(rep(1,ncol(Bl))))) * kronecker(t(as.vector(rep(1,ncol(Bm)))),
                                                           Bl)


knot_grid <- expand.grid(m=interior.knots.m,l=interior.knots.l)[,2:1]
keep <- (knot_grid$m >= 0.5) & (knot_grid$m < max(knot_grid$l)-(0.5*knot_grid$l)) |
      (knot_grid$m <= 0.5) & (knot_grid$m > min(knot_grid$l)+(0.5*knot_grid$l))
knot_grid <- transform(knot_grid,keep=factor(keep))
knot_grid <- data.frame(knot_grid,
                        expand.grid(m_index=1:length(interior.knots.m),
                                    l_index=(1:length(interior.knots.l)))[,2:1])

basisKeepIndex <- knot_grid$keep[!duplicated(knot_grid[,c("l_index","m_index")])] %>% equals(TRUE) %>% which
B. <- B.[,basisKeepIndex]
knot_grid <- subset(knot_grid,keep==TRUE)


dl <- 3
if(dl>0){
      if(sum(knot_grid$m==unique(knot_grid$m)[1])>dl){
            Pl <- matrix(data=0,nrow=sum(knot_grid$m==unique(knot_grid$m)[1])-dl,
                         ncol=nrow(knot_grid))
            Pl[,which(knot_grid$m==unique(knot_grid$m)[1])] <- diff(diag(sum(knot_grid$m==unique(knot_grid$m)[1])),
                                                                    differences = dl)            
            for(i in 2:length(unique(knot_grid$m))){
                  if (sum(knot_grid$m==unique(knot_grid$m)[i]) > dl) {
                        pl <- matrix(data=0,nrow=sum(knot_grid$m==unique(knot_grid$m)[i])-dl,
                                     ncol=nrow(knot_grid))
                        pl[,which(knot_grid$m==unique(knot_grid$m)[i])] <- diff(diag(sum(knot_grid$m==unique(knot_grid$m)[i])),
                                                                                differences = dl)      
                        Pl <- rbind(Pl,pl)
                  } 
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

lambdas <- expand.grid(lam_l=exp(seq(-1,5.3,length.out=20)),
                       lam_m=exp(seq(-1,5.3,length.out=20)))

clusterExport(cl,c("grid",
                   "bsplbase",
                   "fit_cholesky_PS",
                   "Sigma",
                   "N",
                   "m",
                   "M",
                   "Bl",
                   "Bm",
                   "B.",
                   "Pl",
                   "Pm",
                   "lambdas",
                   "dl",
                   "dm"))

nsim <- 50
startTS <- Sys.time()
PS_fit_sim <- foreach(icount(nsim),.noexport = c("y",
                                                 "y_vec",
                                                 "X",
                                                 "U.")) %dopar% {
                                                       y <- mvrnorm(n=N,mu=rep(0,m),Sigma=Sigma)
                                                       y_vec <- as.vector(t(y[,-1]))
                                                       
                                                       X <- matrix(data=0,nrow=N*(M-1),ncol=choose(M,2))
                                                       no.skip <- 0
                                                       for (t in 2:M){
                                                             X[((0:(N-1))*(M-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- y[,1:(t-1)]
                                                             no.skip <- no.skip + t - 1
                                                       }
                                                       
                                                       U. <- X%*%B.
                                                       fit_list <- list.zip(lam_l=lambdas$lam_m,
                                                                            lam_m=lambdas$lam_l) %>%
                                                             lapply(.,function(l){
                                                                   fit_cholesky_PS(y,U.,
                                                                                   D=diag(diag(C)),
                                                                                   Pl,l$lam_l,
                                                                                   Pm,l$lam_m,
                                                                                   0.0000001)    
                                                             })
                                                       
                                                       
                                                       for(i in 1:length(fit_list)){
                                                             fit_list[[i]]$lam_l <- lambdas$lam_l[i]
                                                             fit_list[[i]]$lam_m <- lambdas$lam_m[i]
                                                       }
                                                       fit_list
                                                       
                                                 }
endTS <- Sys.time()
endTS-startTS





timeStamp <- Sys.time() %>% str_sub(.,start=1,end=19) %>% str_replace_all(.," ","_") %>% str_replace_all(.,":","-")
save(PS_fit_sim,
     file = file.path(getwd(),
                      "..",
                      "data",
                      "moving average",
                      paste0("MA1_fits_dl_",
                             dl,
                             "_dm_",
                             dm,"_N_",
                             N,"_M_",
                             M,
                             "_",
                             timeStamp,
                             ".RData")))



