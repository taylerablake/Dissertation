
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

X <- matrix(data=0,nrow=N*(M-1),ncol=choose(M,2))
no.skip <- 0
for (t in 2:M){
      X[((0:(N-1))*(M-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- y[,1:(t-1)]
      no.skip <- no.skip + t - 1
}

 