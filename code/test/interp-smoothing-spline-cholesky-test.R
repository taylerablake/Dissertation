


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
for (this_t in 2:nrow(T_mod)) {
  for (this_s in 1:(this_t-1)) {
    T_mod[this_t,this_s] <- -0.8*((M-(this_t-this_s)+1)/M)     
  }
}

wireframe(diag(M)-T_mod)
L_mod <- solve(T_mod)

D <- diag(rep(0.1,N*(M-1)))
Sigma <- L_mod %*%(diag(rep(0.1,M))^2) %*% t(L_mod)


y <- mvrnorm(n=N,mu=rep(0,M),Sigma=Sigma)
y_vec <- as.vector(t(y[,-1]))
y_aug <- as.vector(c(y_vec,
                     rep(1,
                         length(unique(grid$l)))))


heatColors <- heat.colors(15)
plot(y[,19],y[,20],col=heatColors[1],
     xlab=expression(y[i]),
     ylab=expression(y[20]))
for (i in 2:10) {
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
K11 <- R1_L + R1_M + R1_LM



### Construct K12 #####################################################################

R1_L <- matrix(0,nrow=nrow(grid),
               ncol=length(unique(grid$m)))
for ( this_row in 1:nrow(R1_L) ) {
  for ( this_col in 1:ncol(R1_L) ) {
    R1_L[this_row,this_col] <- R1(grid$l[this_row]/max(grid$l),
                                  0,
                                  m=2)
  }
}


R1_M <- matrix(0,nrow=nrow(grid),
               ncol=length(unique(grid$m)))
for ( this_row in 1:nrow(R1_M) ) {
  for ( this_col in 1:ncol(R1_M) ) {
    R1_M[this_row,this_col] <- R1(grid$m[this_row]/max(grid$m),
                                  unique(grid$m[this_col])/max(grid$m),
                                  m=1)
  }
}


R1_LM <- R1_L*R1_M + outer(k1(grid$l/max(grid$l)),k1(rep(0,length(unique(grid$m)))))*R1_M
K12 <- R1_L + R1_M + R1_LM
K21 <- t(K12)


### Construct K22 #####################################################################

R1_L <- matrix(0,nrow=length(unique(grid$m)),
               ncol=length(unique(grid$m)))
for ( this_row in 1:nrow(R1_L) ) {
  for ( this_col in 1:ncol(R1_L) ) {
    R1_L[this_row,this_col] <- R1(0,
                                  0,
                                  m=2)
  }
}


R1_M <- matrix(0,nrow=length(unique(grid$m)),
               ncol=length(unique(grid$m)))
for ( this_row in 1:nrow(R1_M) ) {
  for ( this_col in 1:ncol(R1_M) ) {
    R1_M[this_row,this_col] <- R1(unique(grid$m[this_row])/max(grid$m),
                                  unique(grid$m[this_col])/max(grid$m),
                                  m=1)
  }
}


R1_LM <- R1_L*R1_M + outer(k1(rep(0,length(unique(grid$m)))),
                           k1(rep(0,length(unique(grid$m)))))*R1_M
K22 <- R1_L + R1_M + R1_LM

K <- rbind( cbind(K11,K12),
            cbind(K21,K22))














B1 <- matrix(data=c(rep(1,nrow(grid)),
                   k1(grid$l/max(grid$l))),
            ncol=2,
            nrow=nrow(grid),
            byrow=FALSE)
B2 <- matrix(data=c(rep(1,length(unique(grid$l))),
                    k1(rep(0,length(unique(grid$l))))),
             ncol=2,
             nrow=nrow(grid),
             byrow=FALSE)
B <- rbind(B1,B2)


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










c <- my.Q2%*%solve( t(my.Q2)%*%(K + lambdas[lambda.i]*solve( t(Z.1)%*%my.D.1.inv%*%Z.1 ))%*%my.Q2 )%*%t(my.Q2)%*%solve( t(Z.1)%*%my.D.1.inv%*%Z.1 )%*%t(Z.1)%*%my.D.1.inv%*%Y1.aug
d <- my.R.inv%*%t(my.Q1)%*%( solve( t(Z.1)%*%my.D.1.inv%*%Z.1 )%*%t(Z.1)%*%my.D.1.inv%*%Y1.aug - ( K + lambdas[lambda.i]*solve( t(Z.1)%*%my.D.1.inv%*%Z.1 ) )%*%my.a1mat[,lambda.i] )







