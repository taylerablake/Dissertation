


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

N <- 100
M <- 20
grid <- expand.grid(s=(1:M),t=(1:M)) %>%
      subset(.,s>t) %>%
      transform(l=(s-t),m=(s+t)/2)


## start simple: phi(t,s) = 0.8*(t-s)
T_mod <- diag(M)
for (this_t in 2:nrow(T_mod)) {
      for (this_s in 1:(this_t-1)) {
            T_mod[this_t,this_s] <- -0.9*((M-(this_t-this_s)+1)/(4*M))    
      }
}

wireframe(diag(M)-T_mod)
L_mod <- solve(T_mod)

D <- diag(rep(0.1,N*(M-1)))
#D <- diag(rep(seq(0.1,0.01,length.out=M)[-1],N))
Sigma <- L_mod %*%(diag(rep(0.1,M))^2) %*% t(L_mod)


y <- mvrnorm(n=N,mu=rep(0,M),Sigma=Sigma)
y_vec <- as.vector(t(y[,-1]))


heatColors <- heat.colors(20)
plot(y[,19],y[,20],col=heatColors[1],
     xlab=expression(y[i]),
     ylab=expression(y[20]))
for (i in 2:15) {
      points(y[,20-i],y[,20],col=heatColors[i])      
}

plot(y[,14],
     y[,15],col=heatColors[1],
     xlab=expression(y[i]),
     ylab=expression(y[15]))
for (i in 2:10) {
      points(y[,15-i],y[,15],col=heatColors[i])      
}



X <- matrix(data=0,nrow=N*(M-1),ncol=choose(M,2))
no.skip <- 0
for (t in 2:M){
      X[((0:(N-1))*(M-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- y[,1:(t-1)]
      no.skip <- no.skip + t - 1
}


bPars <- rbind(c(0,1,length(unique(grid$l)),3,100,3),
               c(0,1,length(unique(grid$m)),3,100,3))


Bl <- bsplbase(as.vector(c(rep(0,M-1),grid$l)/max(grid$l)),
               bPars[1,],outer.okay = TRUE)$base
Bm <- bsplbase(as.vector(c(unique(grid$s)/max(grid$s),grid$m)/max(grid$m)),
               bPars[2,],outer.okay = TRUE)$base

n1 <- ncol(Bl)
n2 <- ncol(Bm)

knots.l <- bsplbase(as.vector(grid$l)/max(grid$l), bPars[1,])$knots
interior.knots.l <- knots.l[1:(length(knots.l)-(bPars[1,4]+1))]

knots.m <- bsplbase(as.vector(grid$m)/max(grid$m), bPars[2,])$knots
interior.knots.m <- knots.m[1:(length(knots.m)-(bPars[1,4]+1))]

B. <- kronecker(Bm,t(as.vector(rep(1,ncol(Bl))))) * kronecker(t(as.vector(rep(1,ncol(Bm)))),
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



dl <- 2
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


y_aug <- c(rep(1,M-1))
X_aug <- rbind(cbind(diag(M-1),
                     matrix(0,nrow=M-1,
                            ncol=(nrow(B.)-(M-1)))),
               cbind(matrix(0,nrow=nrow(X),
                            ncol=(nrow(B.)-ncol(X))),
                     X)) 
U. <- X_aug%*%B.




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
      fit_cholesky_PS(y,
                      y_aug=rep(1,M-1),
                      U.,
                      D=diag(rep(0.1,M)),
                      Pl,
                      lambda_pair$lam_l,
                      Pm,
                      lambda_pair$lam_m,
                      0.0000001)    
    }
  
  
  for(i in 1:length(fit_list)){
    fit_list[[i]]$lam_l <- lambdas$lam_l[i]
    fit_list[[i]]$lam_m <- lambdas$lam_m[i]
  }




gg <- expand.grid(s=(1:M),t=(1:M)) %>% subset(.,t>s)
gg <- transform(gg,l=t-s,m=s+t)

l. <- gg$l
Bl <- bsplbase(l./max(gg$l), bPars[1,  ])$base
knots.l <- bsplbase(l./max(gg$l), bPars[1,  ])$knots
knots.l <- knots.l[1:(length(knots.l)-(bPars[1,4]+1))]

m. <- gg$m
Bm <- bsplbase(m./max(gg$m), bPars[2,  ])$base
knots.m <- bsplbase(m./max(gg$m), bPars[2,  ])$knots
knots.m <- knots.m[1:(length(knots.m)-(bPars[2,4]+1))]

knot_grid <- expand.grid(m=knots.m,l=knots.l)[,2:1]
keep <- (knot_grid$m >= 0.5) & (knot_grid$m < max(knot_grid$l)-0.5*knot_grid$l) |
      (knot_grid$m <= 0.5) & (knot_grid$m > min(knot_grid$l)+0.5*knot_grid$l)
knot_grid <- transform(knot_grid,keep=factor(keep))
knot_grid <- data.frame(knot_grid,
                        expand.grid(m_index=1:length(knots.m),
                                    l_index=(1:length(knots.l)))[,2:1])

B_lm <- t(sapply(1:nrow(Bl),function(i){as.vector(outer(Bm[i,],Bl[i,]))}))
basisKeepIndex <- knot_grid$keep[!duplicated(knot_grid[,c("l_index","m_index")])] %>% equals(TRUE) %>% which
knot_grid <- subset(knot_grid,keep==TRUE)
B_lm <- B_lm[,basisKeepIndex]

phi_list <- lapply(fit_list,function(l){
      B_lm%*%l$coef
})
i <- 1

phi <- diag(M)
phi[lower.tri(phi)] <- phi_list[[i]]
laml <- round(lambdas$lam_l[i],3)
lamm <- round(lambdas$lam_m[i],3)
my_title <- c(as.expression(bquote(lambda[l] == .(laml))), as.expression(bquote(lambda[m] == .(lamm))))
wireframe(phi,
          scales = list(arrows = FALSE),
          par.settings = list(axis.line = list(col = "transparent")),
          main=my_title)

i <- i+1



