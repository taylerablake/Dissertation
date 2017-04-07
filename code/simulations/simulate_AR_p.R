


library(MASS)
library(magrittr)
library(rlist)
library(plyr)
library(dplyr)
M <- m <- 30
N <- 50
grid <- expand.grid(t=1:M,s=1:M) %>% subset(.,t>s) %>%
      transform(.,l=t-s,
                m=(t+s)/2)



list_Sigma <- list()
phi <- 0.8
sig2 <- 0.1

Sigma_AR_1 <- matrix(data=NA,nrow=M,
                  ncol=M)
for (i in 1:nrow(Sigma_AR_1)) {
      for (j in 1:ncol(Sigma_AR_1)) {
            Sigma_AR_1[i,j] <- sig2*( phi^(abs(i-j)) )/( 1 - phi^2 )
      }
}
list_Sigma[[2]] <- Sigma_AR_1

Omega <- solve(Sigma_AR_1)

C <- t(chol(Sigma_AR_1))
D <- diag(diag(C))
L <- C%*%solve(D)
T_mod <- solve(L)      

matplot(t(mvrnorm(n=20,mu=rep(0,m),Sigma=Sigma_AR_1)),col="pink",type="l")
matlines(sapply(1:20,function(i){
      arima.sim(n = M,
                list(ar = c(0.8),
                     ma = c(0)),
                sd = sqrt(0.1))
}),
col="blue")

bPars <- rbind(c(0,1,50,3,100,2),
               c(0,1,60,3,100,3))


y <- mvrnorm(n=N,mu=rep(0,m),Sigma=Sigma_AR_1)
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
dl <- 2
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










lambdas <- expand.grid(lam_l=exp(seq(-2,7,length.out=15)),
                       lam_m=exp(seq(-2,7,length.out=15)))

ar_1_coef_list <- list.zip(lam_l=lambdas$lam_m,
                      lam_m=lambdas$lam_l) %>%
      lapply(.,function(l){
            try(            fit_cholesky_PS(y_vec,U.,Pl,l$lam_l,
                                            Pm,l$lam_m,
                                            0.12)    
            )
      })

setwd("/Users/taylerblake/Documents/Dissertation/code/simulations")
save(ar_1_coef_list,file="stationary_ar1_coef_list.Rdata")



data.frame(expand.grid(t=1:M,s=1:M),
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






ar1_phi_list <- lapply(ar_1_coef_list,function(beta){
      B_lm%*%beta
})
i <- 1


#Phi_hat <- B_lm %*% f$coefficients
# par(ask=TRUE)
#for (i in 76:nrow(lambdas)){
laml <- round(lambdas$lam_l[i],3)
lamm <- round(lambdas$lam_m[i],3)
my_title <- c(as.expression(bquote(lambda[l] == .(laml))), as.expression(bquote(lambda[m] == .(lamm))))
Phi_hat <- ar1_phi_list[[i]]
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

