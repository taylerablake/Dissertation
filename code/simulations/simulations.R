



ti(x,z)



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





grid <- expand.grid(t1=(1:20)/20,t2=(1:20)/20)%>%
      subset(.,t1>t2)
t <- 1:20

indices_of_nonzeros <- as.matrix(expand.grid(t1=t,t2=t) %>% subset(.,(t1-t2)==1))
nonzero_phis <- (2*((2:length(t))/length(t))^2)-0.5
T_mat <- diag(rep(1,length(t)))
T_mat[indices_of_nonzeros] <- -nonzero_phis
#sigma <- log(t/10+2)
y <- solve(T_mat)%*%as.vector(rnorm(length(t),mean=0,sd=log(2.5)))

true_Sigma <- solve(T_mat)%*%T_mat

M1.index <- grid[,1] %>% unique
M2.index <- grid[,2] %>% unique


X <- sapply(2:length(t),function(i){
            c(y[1:(i-1)],rep(0,length(t)-i))
      }) %>% t
X <- matrix(data=0,nrow=length(t),ncol=nrow(grid))
X[2,1] <- y[1]
column_index <- 1
for(i in 3:length(y)){
      X[i,((column_index+1):(column_index+i-1))] <- y[1:(i-1)]
      column_index <- column_index+ i-1 
}
X

#################################################################################

B1 <- bsplbase(M1.index, c(0,1,20,3))
B2 <- bsplbase(M2.index, c(0,1,20,3))
U <- X%*%B1

lambda <- 1000
lambdaBar <- 10000000

d <- 2
#D <- diff(diag(ncol(B1)), diff = d)
D <- diag(ncol(B1))
P <- kronecker(diag(rep(1,ncol(B2))), t(D)%*%D )
dBar <- 2
Dbar <- diff(diag(ncol(B2)), diff = dBar)
Pbar <- kronecker(t(Dbar)%*%Dbar, diag(rep(1,ncol(B1))) )

n1 <- ncol(B1)
n2 <- ncol(B2)	
# Compute tensor products for estimated alpha surface
#B1. <- kronecker(B1, t(rep(1, n2)))
U. <- kronecker(U, t(rep(1, n2)))
B2. <- kronecker(t(rep(1, n1)), B2)
rm(B1)
rm(B2)
rm(oM1)
rm(oM2)
#B. <- B1. * B2.
Q. <- U. * B2.
dim(U.)
dim(B2.)

a <- solve(t(Q.)%*%Q. + lambda*P + lambdaBar*Pbar,t(Q.)%*%y[-1])




#################################################################################
#################################################################################





