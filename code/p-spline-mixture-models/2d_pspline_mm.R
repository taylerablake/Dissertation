source(file.path(getwd(),"fnc","bsplbase.R"))


M <- 30

grid <- expand.grid(t=1:M,s=1:M) %>%
  subset(.,t>s) %>%
  transform(.,l=t-s,
            m=(t+s)/2)

bPars <- rbind(c(0,1,20,3,100,3),
               c(0,1,20,3,100,3))

Bl <- bsplbase(as.vector(grid$l)/max(grid$l), bPars[1,])$base
#Bl <- bs(as.vector(grid$l),df=length(unique(grid$l)))
#Bm <- bs(as.vector(grid$m)/max(grid$m),df=length(unique(grid$m)))
Bm <- bsplbase(as.vector(grid$m)/max(grid$m), bPars[2,])$base
knots_l <- bsplbase(as.vector(grid$l)/max(grid$l), bPars[1,])$knots
knots_m <- bsplbase(as.vector(grid$m)/max(grid$m), bPars[2,])$knots
B. <- kronecker(Bm,t(as.vector(rep(1,ncol(Bl))))) * kronecker(t(as.vector(rep(1,ncol(Bm)))),Bl)

dl <- 2
Dl <- diag(ncol(Bl))
if(dl != 0) {
      for(j in 1:dl) {
            Dl <- diff(Dl)
      }
}
lambdal <- 10
#Pl <- sqrt(lambdal) * kronecker(Dl, diag(ncol(Bm)))


dm <- 1
Dm <- diag(ncol(Bm))
if(dm != 0) {
      for(j in 1:dm) {
            Dm <- diff(Dm)
      }
}
lambdam <- 20
#Pm <- sqrt(lambdam) * kronecker(diag(ncol(Bl)), Dm)





Pl <- t(Dl)%*%Dl
Pl_SVD <- eigen(Pl,symmetric=TRUE)
U_ls <- Pl_SVD$vectors[,1:(ncol(Dl)-dl)]
U_ln <- Pl_SVD$vectors[,(ncol(Dl)-dl+1):ncol(Dl)]
  T_l <- cbind(U_ln,U_ls)

  X_l <- Bl %*% U_ln 
  Z_l <- Bl %*% U_ls

  
  
Pm_SVD <- eigen(t(Dm)%*%Dm)
U_ms <- Pm_SVD$vectors[,(dm+1):ncol(Dm)]
U_mn <- Pm_SVD$vectors[,1:dm]
      T_m <- cbind(U_mn,U_ms)
      X_m <- Bm %*% U_mn 
      Z_m <- Bm %*% U_ms


      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      l. <- as.vector(outer(rep(1,length(unique(grid$m))),unique(grid$l)))
      Bl <- bsplbase(l./max(grid$l), bPars[1,  ])$base
      knots.l <- bsplbase(l./max(grid$l), bPars[1,  ])$knots
      knots.l <- knots.l[1:(length(knots.l)-(bPars[1,4]+1))]
      
      
      
      
      m. <- as.vector(outer(unique(grid$m),rep(1,length(unique(grid$l)))))
      Bm <- bsplbase(m./max(grid$m), bPars[2,  ])$base
      knots.m <- bsplbase(m./max(grid$m), bPars[2,  ])$knots
      knots.m <- knots.m[1:(length(knots.m)-(bPars[2,4]+1))]
      
      
      
      
      
      
      
      