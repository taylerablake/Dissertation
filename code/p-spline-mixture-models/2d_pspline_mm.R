


Bl <- bsplbase(as.vector(grid$l)/max(grid$l), bPars[1,])$base
Bm <- bsplbase(as.vector(grid$m)/max(grid$m), bPars[2,])$base
B. <- kronecker(Bm,t(as.vector(rep(1,ncol(Bl))))) * kronecker(t(as.vector(rep(1,ncol(Bm)))),Bl)

dl <- 2
Dl <- diag(ncol(Bl))
if(dl != 0) {
      for(j in 1:dl) {
            Dl <- diff(Dl)
      }
}
lambdal <- 10
Pl <- sqrt(lambdal) * kronecker(Dl, diag(ncol(Bm)))


dm <- 1
Dm <- diag(ncol(Bm))
if(dm != 0) {
      for(j in 1:dm) {
            Dm <- diff(Dm)
      }
}
lambdam <- 20
Pm <- sqrt(lambdam) * kronecker(diag(ncol(Bl)), Dm)








Pl_SVD <- eigen(t(Dl)%*%Dl)
U_ls <- Pl_SVD$vectors[,which(Pl_SVD$values > sort(Pl_SVD$values)[2])]
U_ln <- Pl_SVD$vectors[,which(Pl_SVD$values <= sort(Pl_SVD$values)[2])]
      T_l <- cbind(U_ln,U_ls)

Pm_SVD <- eigen(t(Dm)%*%Dm)
U_ms <- Pm_SVD$vectors[,which(Pm_SVD$values > sort(Pm_SVD$values)[1])]
U_mn <- Pm_SVD$vectors[,which(Pm_SVD$values == sort(Pm_SVD$values)[1])]
      T_m <- cbind(U_mn,U_ms)



X_l <- Bl %*% U_ln 
Z_l <- Bl %*% U_ls

X_m <- Bm %*% U_mn 
Z_m <- Bm %*% U_ms

