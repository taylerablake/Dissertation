
library(splines)
library(magrittr)
library(ggplot2)
library(plyr)
library(dplyr)

M <- 30

grid <- expand.grid(t=1:M,s=1:M) %>% subset(.,t>s) %>%
      transform(.,l=t-s,
                m=(t+s)/2)

bPars <- rbind(c(1/29,1,50,3,100,2),
               c(1/59,1,60,3,100,3))



Bl <- bsplbase((grid$l/(M-1)),bPars[1,])$base
Bm <- bsplbase((grid$m/((2*M)-1)),bPars[2,])$base
lKnots <- bsplbase((grid$l/(M-1)),bPars[1,])$knots
mKnots <- bsplbase((grid$m/((2*M)-1)),bPars[2,])$knots

knotGrid <- expand.grid(m=mKnots,l=lKnots)[,2:1]
keep <- (knotGrid$m >= 0.5) & (knotGrid$m < max(knotGrid$l)-0.5*knotGrid$l) |
            (knotGrid$m < 0.5) & (knotGrid$m > min(knotGrid$l)+0.5*knotGrid$l)
knotGrid <- transform(knotGrid,keep=factor(keep))
ggplot(knotGrid,aes(x=l,y=m)) + geom_point(aes(colour=keep))
which(knotGrid$keep==TRUE)
m <- 50
grid <- expand.grid(t=1:m,s=1:m)%>%orderBy(~t,.)%>%
      subset(.,t>s)%>%
      transform(.,l=t-s,m=t+s)

rownames(grid) <- 1:nrow(grid)
phi <- function(t1.index,t2.index,m){
      if(t1.index-t2.index==1){
            garp <- (2*((t1.index)/m)^2)-0.5
      }
      if(t1.index-t2.index!=1){
            garp <- 0
      }
      garp
}

true_phi <- sapply(1:nrow(grid),function(row.i){
      phi(grid$t[row.i],grid$s[row.i],m=m)
})

grid <- transform(grid,true_phi=true_phi)
rm(true_phi)


indices_of_nonzeros <- as.matrix(expand.grid(t=(1:m),s=(1:m)) %>% subset(.,(t-s)==1))
nonzero_phis <- (2*((2:m)/m)^2)-0.5

T_mat <- diag(rep(1,m))
phis <- as.vector(rep(0,sum(lower.tri(T_mat))))
T_mat[indices_of_nonzeros] <- -nonzero_phis
phis <- -T_mat[lower.tri(T_mat)]


N <- 1
y <- solve(T_mat)%*%matrix(data=rnorm(N*m,mean=0,sd=1),nrow=m,ncol=N)
true_Sigma <- solve(T_mat)%*%T_mat



l. <- as.vector(outer(rep(1,length(unique(grid$m))),(unique(grid$l))))
Bl <- bsplbase(l./max(grid$l), Pars[1,  ])
knots.l <- seq(Pars[1,1] - Pars[1,4] * ((Pars[1,2] - Pars[1,1])/Pars[1,3]),
               Pars[1,2] + Pars[1,4] * ((Pars[1,2] - Pars[1,1])/Pars[1,3]),
               by =((Pars[1,2] - Pars[1,1])/Pars[1,3]))

Bl <- matrix(rep(knots.l,3),nrow=3,ncol=length(knots.l),byrow=TRUE)


m. <- as.vector(outer(unique(grid$m),rep(1,length(unique(grid$l)))))
Bm <- bsplbase(m./max(grid$m), Pars[2,  ])
knots.m <- seq(Pars[2,1] - Pars[2,4] * ((Pars[2,2] - Pars[2,1])/Pars[2,3]),
               Pars[2,2] + Pars[2,4] * ((Pars[2,2] - Pars[2,1])/Pars[2,3]),
               by =((Pars[2,2] - Pars[2,1])/Pars[2,3]))
Bm <- matrix(rep(knots.m,3),nrow=3,ncol=length(knots.m),byrow=TRUE)






n1 <- ncol(Bl)
n2 <- ncol(Bm)  # Compute tensor products

Bl. <- kronecker(Bl, t(rep(1, n2)))
Bm. <- kronecker(t(rep(1, n1)), Bm)


discard_rows <- left_join(data.frame(l=l.,m=m.,frame="big"),
                          data.frame(grid[,c("l","m")],frame="little")
                          ,by=c("l","m"))$frame.y %>% is.na %>% which

B. <- Bl.*Bm.

knot.grid <- expand.grid(l.knot=knots.l,m.knot=knots.m)
knot.grid <- knot.grid %>% transform(.,t=(l.knot+m.knot)/2,
                                     s=(m.knot-l.knot)/2)
plot(knot.grid$s[knot.grid$s <knot.grid$t],knot.grid$t[knot.grid$s <knot.grid$t],xlab="s",ylab="t")


eps <- max(Pars[1,4] * ((Pars[1,2] - Pars[1,1])/Pars[1,3]),Pars[2,4] * ((Pars[2,2] - Pars[2,1])/Pars[2,3]))
knot.grid <- subset(knot.grid,(s > -eps) & (s < (1+eps)) & (s < t) &(t > -eps) & (t < (1+eps)))
plot(knot.grid$s,knot.grid$t,xlab="s",ylab="t")




reference_grid <- data.frame(l=l.,m=m.)[-discard_rows,]
match_row <- function(row,df){
      match.ind <- match(df[,1],row[1])
      for(j in 2:ncol(df)){
            match.ind <- match.ind*match(df[,j],row[j])
      }
      which(!is.na(match.ind))
}
x_col_ind <-   sapply(1:nrow(grid),function(row.i){
      match_row(grid[row.i,c("l","m")],reference_grid)
})
grid <- transform(grid,column_index=x_col_ind)
X <- matrix(data=0,nrow=m-1,ncol=nrow(grid))
for(row.i in 1:nrow(grid)){
      X[grid$t[row.i]-1,grid$column_index[row.i]] <- as.vector(y)[grid$s[row.i]] 
}
U. <- X%*%B.[-discard_rows,]


d1 <- Pars[1, 6]
#d1 <- 1
D1 <- diag(n1)
if(d1 != 0) {
      for(j in 1:d1) {
            D1 <- diff(D1)
      }
}
lambda1 <- Pars[1, 5]
P1 <- sqrt(lambda1) * kronecker(D1, diag(n2))
d2 <- Pars[2, 6]
#d2 <- 1
D2 <- diag(n2)
if(d2 != 0) {
      for(j in 1:d2) {
            D2 <- diff(D2)
      }
}
lambda2 <- Pars[2, 5]
P2 <- sqrt(lambda2) * kronecker(diag(n1), D2)
Pen <- rbind(P1, P2)


###########################################################################

z1 <- rep(0, n2 * (n1 - d1))
z2 <- rep(0, n1 * (n2 - d2))

n.col <- ncol(U.)
nix <- c(z1, z2)
nix.ridge <- rep(0, n1 * n2)
coef.est <- rep(1, ncol(U.))

mu <- rep(mean(y[-1]), length(y[-1]))
it <- 0
repeat {
      if(it == 0) {
            eta <- mu
      }
      
      it <- it + 1
      if(it > 25)
            break
      
      mu <- eta
      h.prime <- 1
      w <- rep(1, length(y[-1]))
      u <- (y[-1] - mu)/h.prime + eta
      
      
      library(inline)
      library(Rcpp)
      library(RcppArmadillo)
      rcpp_inc <- 'using namespace Rcpp;
      using namespace arma;
      '
      
      ## Rcpp code to calculate the matrix pseudo inverse 
       src <- ' mat m1 = as<mat>(m1in);
       //return(wrap(m1.i()));
       
       mat m2 = inv(m1);
       return(wrap(m2));
       '
       fn <- cxxfunction(signature(m1in="numeric"), src, plugin='RcppArmadillo', rcpp_inc)
       
       m3 <- matrix(runif(16), nr=4)
       res <- fn(m3)
       test <- solve(m3)
       all.equal(test, res)
      
      
      f <- lsfit(rbind(U.,Pen), c(u, nix), wt = c(w, nix + 1) *
                       c(w, (nix + 1)), intercept = F)
      coef.old <- coef.est
      coef.est <- as.vector(f$coef)
      d.coef <- max(abs((coef.est - coef.old)/coef.old))
      if(d.coef < 1e-008)
            break
      print(c(it, d.coef))
      eta <- U. %*% coef.est
}
if(it > 24) {
      warning(paste("parameter estimates did NOT converge in 25 iterations"
      ))
}


H <- hat(f$qr, intercept = F)[1:(m-1)]
trace <- eff.dim <- sum(H)

cv <- press.mu <- press.e <- var.c <- NULL
dev <- sum(f$residuals[1:(m-1)]^2)
dispersion.parm <- dev/((m-1) - trace)
press.e <- f$residuals[1:(m-1)]/(1 - H)
cv <- sqrt(sum((press.e)^2)/(m-1))
press.mu <- y[-1] - press.e

aic <- dev + 2 * trace
w.aug <- c(w, (c(z1, z2) + 1))

Phi_hat <- B. %*% f$coefficients
pl <- length(unique(grid$l))
pm <- length(unique(grid$m))

Phi_hatm <- matrix(Phi_hat, pl, pm, byrow = T)
i.1 <- 1:length(unique(grid$l))/max(grid$l)
in1 <- 1:length(M1.index)
i.2 <- M2.index
in2 <- 1:length(M2.index)
if(length(M1.index) > 100) {
      i.1 <- round(seq(from = min(M1.index), to = max(
            M1.index), length = 100))
      in1 <- round(seq(from = 1, to = p1, length = 50))
}
if(length(M2.index) > 100) {
      i.2 <- round(seq(from = min(M2.index), to = max(
            M2.index), length = 100))
      in2 <- round(seq(from = 1, to = p2, length = 50))
}
persp(M2.index[in2], (M1.index[in1]), t(A.hatm[in1, in2]), xlab
      = y.lab, ylab = x.lab, zlab = z.lab)




###########################################################################