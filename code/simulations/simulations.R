



library(ggthemes)



test1 <- function(x,z,sx=0.3,sz=0.4) { 
      x <- x*20
      (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
                         0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))
}
n <- 500
old.par <- par(mfrow=c(2,2))
x <- runif(n)/20;z <- runif(n);
xs <- seq(0,1,length=30)/20;zs <- seq(0,1,length=30)
pr <- data.frame(x=rep(xs,30),z=rep(zs,rep(30,30)))
truth <- matrix(test1(pr$x,pr$z),30,30)
f <- test1(x,z)
y <- f + rnorm(n)*0.2
b1 <- gam(y~s(x,z))
persp(xs,zs,truth);title("truth")
vis.gam(b1);title("t.p.r.s")
b2 <- gam(y~te(x,z))
vis.gam(b2);title("tensor product")
b3 <- gam(y~ ti(x) + ti(z) + ti(x,z))
vis.gam(b3);title("tensor anova")

## now illustrate partial ANOVA decomp...
vis.gam(b3);title("full anova")
b4 <- gam(y~ ti(x) + ti(x,z,mc=c(0,1))) ## note z constrained!
vis.gam(b4);title("partial anova")
plot(b4)

par(old.par)

















#######################################################################################################
#######################################################################################################


library(gamair)
library(doBy)
data("coast")
data("mack")
head(mack)


mack$log.na <- log(mack$net.area)
mack$t.bd <- (mack$b.depth)^.25
b <- gam(egg.count ~ offset(log.na) + s(lon,lat) + s(lon,lat,by=t.bd)+
               s(lon,lat,by=I(t.bd^2)),
         data=mack,family=Tweedie(p=1.1,link=log),method="ML")
for (i in 1:3) { plot(b,select=i);lines(coast)}




#######################################################################################################
#######################################################################################################




m <- 50
grid <- expand.grid(t1.index=1:m,
                    t2.index=1:m)%>% transform(.,t1=t1.index/m,
                                                t2=t2.index/m)%>%
      orderBy(~t1.index,.)%>%
      subset(.,t1>t2)%>%
      transform(.,l=(t1.index-t2.index)/max(t1.index-t2.index),
                m=0.5*(t1.index+t2.index)/max(0.5*(t1.index+t2.index)))
      

phi <- function(t1.index,t2.index,m){
      if(t1.index-t2.index==1){
            garp <- (2*((t1.index)/m)^2)-0.5
      }
      if(t1.index-t2.index!=1){
            garp <- 0
      }
      garp
}

sapply(1:nrow(grid),function(row.i){
      phi(grid$t1.index[row.i],grid$t2.index[row.i],m=m)
})


indices_of_nonzeros <- as.matrix(expand.grid(t1=t,t2=t) %>% subset(.,(t1-t2)==1))
nonzero_phis <- (2*((2:length(t))/length(t))^2)-0.5
phis <- as.vector(rep(0,sum(lower.tri(T_mat))))
T_mat <- diag(rep(1,length(t)))
T_mat[indices_of_nonzeros] <- -nonzero_phis
phis <- -T_mat[lower.tri(T_mat)]

y <- solve(T_mat)%*%as.vector(rnorm(length(t),mean=0,sd=1))
true_Sigma <- solve(T_mat)%*%T_mat

#M1.index <- grid[,1] %>% unique
#M2.index <- grid[,2] %>% unique
l.index <- grid[,3] %>% unique
m.index <- grid[,4] %>% unique







# X <- sapply(2:length(t),function(i){
#             c(y[1:(i-1)],rep(0,length(t)-i))
#       }) %>% t
X <- matrix(data=0,nrow=length(t),ncol=nrow(grid))
X[2,1] <- y[1]
column_index <- 1
for(i in 3:length(y)){
      X[i,((column_index+1):(column_index+i-1))] <- y[1:(i-1)]
      column_index <- column_index+ i-1
}
X <- X[-1,]

#################################################################################

#B1 <- bsplbase(M1.index, c(0,1,20,3))
#B2 <- bsplbase(M2.index, c(0,1,20,3))

Bl <- bsplbase(grid$l, c(0,1,40,3))
Bm <- bsplbase(grid$m, c(0,1,40,3))


#U <- X%*%B1
U <- X%*%Bl

lambda <- 100
lambdaBar <- 1
lambdaRidge <- 0

d <- 2
D <- diff(diag(ncol(Bl)), diff = d)
#D <- diag(ncol(Bl))
P <- kronecker(diag(rep(1,ncol(Bm))), t(D)%*%D )
dBar <- 3
Dbar <- diff(diag(ncol(Bm)), diff = dBar)
Pbar <- kronecker(t(Dbar)%*%Dbar, diag(rep(1,ncol(Bl))) )
Pridge <- diag(rep(1,ncol(Pbar)))
n1 <- ncol(Bl)
n2 <- ncol(Bm)	

Bl. <- kronecker(Bl, t(rep(1, n2)))
Bm. <- kronecker(t(rep(1, n1)), Bm)

B. <- Bl. * Bm.
Q. <- X %*% B.

M <- solve((t(B.) %*% B. + lambda*P + lambdaBar*Pbar + lambdaRidge*Pridge))
a <-  M%*%t(Q.) %*% as.vector(y)
H <- Q.%*%M%*%t(Q.)
ED <- H %>% diag %>% sum

a <- solve(t(Q.)%*%Q. + lambda*P + lambdaBar*Pbar + lambdaRidge*Pridge,t(Q.)%*%y[-1])
y_hat <- as.vector(Q.%*%a)
phi_hat <- as.vector(B.%*%a)

T_hat <- matrix(data=0,nrow=length(y),ncol=length(y))
diag(T_hat) <- rep(1,length(y))
T_hat[as.matrix(grid[,1:2])] <- -phi_hat 

df <- data.frame(e=as.vector(T_mat - T_hat),expand.grid(l=1:nrow(T_mat),m=1:ncol(T_mat)))
df %>% ggplot(.,aes(x=l,y=m)) + geom_point(aes(colour=e^2)) + theme_minimal() + scale_color_continuous_tableau("Green")




df <- data.frame(type=c(rep("observed",length(y)-1),rep("predicted",length(y)-1)),
           value=c(y[-1],y_hat),
           t=rep((2:m)/m,2))
df %>% ggplot(.,aes(x=t,y=value)) + geom_point(aes(colour=type)) + theme_wsj() 
data.frame(t=2:m,e=y[-1]-y_hat) %>% ggplot(.,aes(x=t,y=e)) + geom_point() + theme_wsj() + scale_color_wsj()

      Omega_hat <- T_hat%*%t(T_hat)
      
entropy_loss <- sum(diag(Omega_hat%*%true_Sigma)) - log(det(Omega_hat%*%true_Sigma)) - m
      
      
#################################################################################
#################################################################################





