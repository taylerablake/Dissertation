
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


N <- 100
M <- 20
grid <- expand.grid(s=(1:M),t=(1:M)) %>%
      subset(.,s>t) %>%
      transform(l=(s-t),m=(s+t)/2)


## start simple: phi(t,s) = 0.8*(t-s)
T_mod <- diag(M)
for (this_t in 1:nrow(T_mod)) {
      for (this_s in 1:(this_t-1)) {
       T_mod[this_t,this_s] <- -0.8*((M-(this_t-this_s)+1)/M)     
      }
}
T_mod[1,1] <- 1
wireframe(diag(M)-T_mod)
L_mod <- solve(T_mod)

D <- diag(rep(0.1,N*(M-1)))
Sigma <- L_mod %*%(diag(rep(0.1,M))^2) %*% t(L_mod)


y <- mvrnorm(n=N,mu=rep(0,M),Sigma=Sigma)
y_vec <- as.vector(t(y[,-1]))



heatColors <- heat.colors(15)
plot(y[,19],y[,20],col=heatColors[1],
     xlab=expression(y[i]),
     ylab=expression(y[20]))
for (i in 2:10) {
      points(y[,20-i],y[,20],col=heatColors[i])      
}

plot(y[,14],y[,15],col=heatColors[1],
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


R1_L <- matrix(0,nrow=nrow(grid),
               ncol=nrow(grid))
for ( this_row in 1:nrow(R1_L) ) {
      for ( this_col in 1:ncol(R1_L) ) {
            R1_L[this_row,this_col] <- R1(grid$l[this_row]/max(grid$l),
                                          grid$l[this_col]/max(grid$l),
                                          m=2)
      }
}


R1_M <- matrix(0,nrow=nrow(grid),
               ncol=nrow(grid))
for ( this_row in 1:nrow(R1_M) ) {
      for ( this_col in 1:ncol(R1_M) ) {
            R1_M[this_row,this_col] <- R1(grid$m[this_row]/max(grid$m),
                                          grid$m[this_col]/max(grid$m),
                                          m=1)
      }
}


R1_LM <- R1_L*R1_M + outer(k1(grid$l/max(grid$l)),k1(grid$l/max(grid$l)))*R1_M
K <- R1_L + R1_M + R1_LM

B <- matrix(data=c(rep(1,nrow(grid)),
                   k1(grid$l/max(grid$l))),
            ncol=2,
            nrow=nrow(grid),
            byrow=FALSE)
QR_B <- qr(B,complete=TRUE)
Q_B <- qr.Q(QR_B,complete=TRUE)
Q2_B <- Q_B[,(ncol(B)+1):ncol(Q_B)]
Q1_B <- Q_B[,1:ncol(B)]
R_B.big <- qr.R(QR_B,complete=TRUE)
R_B <- R_B.big[1:ncol(B),]
R_Binv <- solve(R_B)

Dinv <- diag(1/diag(D))



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

lambdas <- as.list(exp(seq(-1,8,length.out=100)))
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
      T_hat <- diag(rep(1,M))
      T_hat[lower.tri(T_hat)] <- -Phi
      list(phi=Phi,T_mat=T_hat,omega=t(T_hat)%*%T_hat)
})

pred_l <- seq(min(grid$l/max(grid$l) ),
              1,
              length.out=100)
pred_m <- seq(min(grid$m/max(grid$m) ),
              1,
              length.out=100)
pred_R1_L <- matrix(0,
                    nrow=length(pred_l),
                    ncol=nrow(grid))
for (this_col in 1:ncol(pred_R1_L)) {
    pred_R1_L[,this_col] <- R1_l(grid$l[this_col]/max(grid$l),pred_l) 
}

pred_R1_M <- matrix(0,
                    nrow=length(pred_m),
                    ncol=nrow(grid))
for (this_col in 1:ncol(pred_R1_M)) {
      pred_R1_M[,this_col] <- R1_m(grid$m[this_col]/max(grid$m),pred_m) 
}

pred_B_L <- matrix(data=c(rep(1,length(pred_l)),
                          k1(pred_l)),
                   nrow=length(pred_l),
                   ncol=2,
                   byrow=FALSE)

phi_l <- list.zip(d_hat=d,c_hat=c) %>%
      lapply(.,function(l) { 
            pred_R1_L %*% l$c_hat + pred_B_L%*%l$d_hat - l$d_hat[1,]
            })
list.cbind(phi_l) %>% matplot(pred_l,.,type="l",col=heat.colors(50))

phi_l0 <- list.zip(d_hat=d,c_hat=c) %>%
      lapply(.,function(l) { 
            pred_B_L%*%l$d_hat - l$d_hat[1,]
      })
list.cbind(phi_l0) %>% matplot(pred_l,.,type="l",col=heat.colors(100))



phi_m <- list.zip(d_hat=d,c_hat=c) %>%
      lapply(.,function(l) { 
            pred_R1_M %*% l$c_hat
      })
list.cbind(phi_m) %>% matplot(.,type="l",col=heat.colors(100))


pred_grid <- expand.grid(t=1:100,s=1:100) %>%
      subset(.,t>s) %>%
      transform(l=t-s,m=0.5*(t+s))
pred_R1_L <- matrix(0,
                    nrow=nrow(pred_grid),
                    ncol=nrow(grid))
for (this_col in 1:ncol(pred_R1_L)) {
      pred_R1_L[,this_col] <- R1_l(grid$l[this_col]/max(grid$l),pred_grid$l) 
}

pred_R1_M <- matrix(0,
                    nrow=nrow(pred_grid),
                    ncol=nrow(grid))
for (this_col in 1:ncol(pred_R1_M)) {
      pred_R1_M[,this_col] <- R1_m(grid$m[this_col]/max(grid$m),pred_grid$m) 
}

pred_R1_LM <- pred_R1_L*pred_R1_M + outer(k1(pred_grid$l/max(pred_grid$l)),k1(grid$l/max(grid$l)))*pred_R1_M
phi_lm <- list.zip(d_hat=d,c_hat=c) %>%
      lapply(.,function(l) { 
            pred_R1_LM %*% l$c_hat
      })
phi <- list.zip(d_hat=d,c_hat=c) %>%
      lapply(.,function(l) { 
            l$d_hat[1,] +
                  (l$d_hat[2,]* k1(pred_grid$l/max(pred_grid$l))) +
                  (pred_R1_L + pred_R1_LM +pred_R1_M) %*% l$c_hat
      })

par(mfrow=c(3,3))
Phi_LM <- diag(max(pred_grid$t))
Phi_LM[lower.tri(Phi_LM)] <- phi_lm[[10]]
wireframe(Phi_LM,
          row.values = 1:100,
          column.values = 1:100)

par(mfrow=c(3,3))
Phi <- diag(max(pred_grid$t))
Phi[lower.tri(Phi)] <- phi[[60]]
wireframe(Phi,
          row.values = 1:100,
          column.values = 1:100)








      this_lambda <- 1
      wireframe(diag(M)-cholesky[[this_lambda]]$T_mat,
                xlab="",
                ylab="",
                zlab="",
                scales=list(arrows=FALSE))
      
      this_lambda <- this_lambda+1

      
      


      
      
      
      


## l + 2m = t-s+t+s = 2t
T_mod <- diag(M)


