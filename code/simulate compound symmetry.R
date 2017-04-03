

theta <- 0.9
MA_samples <- sapply(1:100000,function(i){
      arima.sim(n = 20, list(ma = c(theta)),
                sd = .02)    
}) %>% t





expand.grid(s=(1:5),t=(1:5)) %>%
      subset(.,s>t) %>%
      transform(l=(s-t),m=s+t) %>%
      transform(.,phi=ifelse(l==1,exp(-0.5*l),0))
      
      
      
grid <- expand.grid(s=(1:5),t=(1:5)) %>%
      subset(.,s>t) %>%
      transform(l=(s-t),m=s+t) %>%
      transform(.,phi=(l/max(grid$l))^2)

cholesky_T <- diag(5)
D <- diag(5)
for (row.i in 1:nrow(grid)) {cholesky_T[grid$s[row.i],grid$t[row.i]] <- grid$phi[row.i]}

Omega <- t(cholesky_T)%*%solve(D)%*%cholesky_T
Sigma <- solve(Omega)

Omega
Sigma











library(MASS)

m <- 30
grid <- expand.grid(s=(1:m),t=(1:m)) %>%
      subset(.,s>t) %>%
      transform(l=(s-t),m=s+t) %>%
      transform(.,phi=(l/max(grid$l))^2)


Sigma <- matrix(0.7,nrow=m,ncol=m) + diag(rep(0.3),m)
Omega <- solve(Sigma)


C <- t(chol(Sigma))
D <- diag(diag(C))
L <- C%*%solve(D)
T_mod <- solve(L)



matplot(Y_chol <- L %*% t(mvrnorm(n=100,mu=rep(0,m),Sigma=I(D^2))),col="pink",type="l")
matlines(t(mvrnorm(n=100,mu=rep(0,m),Sigma=Sigma)),col="blue")

grid <- expand.grid(t=1:m,s=1:m) %>%
      transform(l=t-s,m=s+t)

persp(grid$s, grid$t, z=diag(rep(1,nrow(T_mod)))-T_mod, col = color[facetcol], phi = 30, theta = -30)



      data.frame(grid,
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

      
      ggrid <-  data.frame(grid,
                           phi=as.vector(diag(rep(1,nrow(T_mod)))-T_mod))
      ggrid$phi[ggrid$l<=0] <- 0
      plot(subset(ggrid, m==31)$l,subset(ggrid, m==31)$phi,type="l",
           xlab="l",ylab=expression(phi))
      for (grid.m in 16:46){
            lines(subset(ggrid, m==grid.m)$l,subset(ggrid, m==grid.m)$phi,
                  col=grid.m)      
      }
      plot(subset(ggrid, l==1)$m,subset(ggrid, l==1)$phi,type="l",
           xlab="l",ylab=expression(phi))
      for (grid.l in 2:31){
            lines(subset(ggrid, l==grid.l)$m,subset(ggrid, l==grid.l)$phi,
                  col=grid.l)      
      }
      
      transform(expand.grid(t=1:30,s=1:30),
                l=t-s,m=s+t)
      
      wireframe(phi~s+t,
                data=ggrid,
                light.source = c(5,20,10),
                pretty=TRUE,
                scales = list(arrows = FALSE),
                drape=FALSE,
                screen=list(z = 80, x = -70, y = 3),
                row.values=transform(expand.grid(t=1:30,s=1:30),
                                           l=t-s,m=s+t)$s,
                col.values=transform(expand.grid(t=1:30,s=1:30),
                                     l=t-s,m=s+t)$t,
                cex=0.15,
                col="grey",
                par.settings = list(axis.line = list(col = "transparent")))
      
      