
setwd(file.path("/Users","taylerblake","Documents","Dissertation"))
source(file.path(getwd(),"code","fnc","aux","bsplbase.R"))
source(file.path(getwd(),"code","fnc","draw_psplines.R"))

library(plyr)
library(dplyr)
library(rlist)
library(ggplot2)
library(mgcv)
library(tidyr)
library(splines)
library(reshape2)
library(systemfit)
require(graphics)
library(magrittr)



sourceDir <- function(path, trace = TRUE, ...) {
      for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
            if(trace) cat(nm,":")
            source(file.path(path, nm), ...)
            if(trace) cat("\n")
      }
}

sourceDir(file.path(getwd(),"code","fnc","aux"))
sourceDir(file.path(getwd(),"code","fnc"))


bspline<- function(x,xl,xr,ndx,bdeg){
      dx<- (xr-xl) /ndx
      knots<- seq(xl- bdeg*dx,xr+bdeg*dx,by=dx) 
      B<- spline.des(knots,x,bdeg+1,0*x)$design
      B
}



############################################################################################
############################################################################################

png(filename = file.path(getwd(),"Dissertation TeX","img","pspline_pord2_xsmall_lambda.png"))
draw_psplines(nseg=20,bdeg=3,pord=2,lla=-2,nobs=100)
dev.off()
png(filename = file.path(getwd(),"Dissertation TeX","img","pspline_pord2_small_lambda.png"))
draw_psplines(nseg=20,bdeg=3,pord=2,lla=0, nobs=100)
dev.off()

png(filename = file.path(getwd(),"Dissertation TeX","img","pspline_pord2_medium_lambda.png"))
draw_psplines(nseg=20,bdeg=3,pord=2,lla=2.3,nobs=100)
dev.off()

png(filename = file.path(getwd(),"Dissertation TeX","img","pspline_pord2_large_lambda.png"))
draw_psplines(nseg=20,bdeg=3,pord=2,lla=6,nobs=100)
dev.off()

############################################################################################

x = runif(10)
xg = seq(0, 1, length = 500)
set.seed(123)
y = 1.2 + sin(5  * x) + rnorm(10) * 0.2


B <- bbase(x,  xl = 0, xr = 1, nseg = 57, deg = 3)
nb = ncol(B)
D = diff(diag(nb), diff = 2)

lambda <<- 10 ^ 2
P = lambda * t(D) %*% D
a <- solve(t(B) %*% B + P, t(B) %*% y)
a <- as.vector(a)

cols=rainbow(nb)
Bg <<- bbase(xg, xl = 0, xr = 1, nseg = 57, deg = 3)
A = diag(a)
z = Bg %*% a

png(filename = file.path(getwd(),"Dissertation TeX","img","pspline_10obs_60_basis_functions.png"))
plot(x, y,ylim=c(0,2.5),pch="+",xlab="",ylab="") 
matlines(xg, Bg %*% A, type = 'l', lty = 1, lwd = 2, col= cols,
         xlab = '', ylab = '')
lines(xg, z, col = 'grey',lwd=2)        
dev.off()


##=================================================================================================
##=================================================================================================
##=================================================================================================


x = runif(10)
xg = seq(0, 1, length = 500)
set.seed(123)
y = 1.2 + sin(5  * x) + rnorm(10) * 0.2
f <- function(t){1.2+ sin(5*t)}

lambda <- as.list(exp(seq(-1,10,length.out=50)))
B <- bbase(x,  xl = 0, xr = 1, nseg = 57, deg = 3)
D = diff(diag(ncol(B)), diff = 2)

A <- lapply(lambda,function(l){
            P = l * t(D) %*% D
            a <- solve(t(B) %*% B + P, t(B) %*% y)
            a <- as.vector(a)
            a
      })
Bg <- bbase(xg,  xl = 0, xr = 1, nseg = 57, deg = 3)
beta0hat <- lapply(A,function(a){
      Bg %*% a
})
Beta0 <- matrix(data=unlist(beta0hat),nrow=500,ncol=length(lambda),byrow = FALSE)

png(filename = file.path(getwd(),"Dissertation TeX","img","PS_penalty_section_figure_3.png"))
plot(x, y,pch="+",ylim=c(-0.1,2.6),
     #main=expression(paste(hat(beta)[0](t)," with penalty ||",D[2] ,alpha,"|","|"^2)),
     xlab="",ylab="") 
matlines(xg, Beta0, type = 'l', lty = 1, col=terrain.colors(nb,alpha=0.6),
         xlab = '', ylab = '')
lines(xg,f(xg),col="red")
points(x, y,pch="+")
dev.off()




B <- bspline(x,  xl = 0, xr = 1, ndx = 57, bdeg = 4)
Bg <- bspline(xg,  xl = 0, xr = 1, ndx = 57, bdeg = 4)
D = diff(diag(ncol(B)), diff = 3)
A <- lapply(lambda,function(l){
      P = l * t(D) %*% D
      a <- solve(t(B) %*% B + P, t(B) %*% y)
      a <- as.vector(a)
      a
})

beta0hat <- lapply(A,function(a){
      Bg %*% a
})
Beta0 <- matrix(data=unlist(beta0hat),nrow=500,ncol=length(lambda),byrow = FALSE)

png(filename = file.path(getwd(),"Dissertation TeX","img","PS_penalty_section_figure_4.png"))
plot(x, y,pch="+",ylim=c(0.1,2.5),
     #main=expression(paste(hat(beta)[0](t)," with penalty ||",D[3] ,alpha,"|","|"^2)),
     xlab="",ylab="") 
matlines(xg, Beta0, type = 'l', lty = 1, col=terrain.colors(nb,alpha=0.6),
         xlab = '', ylab = '')
points(x, y,pch="+")
lines(xg,f(xg),col="red")
dev.off()









##=================================================================================================
##=================================================================================================
##=================================================================================================


## Demonstrate limiting behaviour of \hat{beta} as lambda -> infinity
## for various orders of the difference penalty

set.seed(123)
#f <- function(x) 2 * sin(pi * x*2)
f <- function(x){1.2 + 2*sin(5  * x)}
x <- runif(100) %>% sort
xg = seq(0, 1, length = 500)
e <- rnorm(length(x), 0, 0.2)
y <- f(x) + e


lambda <- as.list(10^(seq(-3.5,10,length.out=30)))
B <- bbase(x,  xl = 0, xr = 1, nseg = 67, deg = 3)

for(d in 0:3){
      if(d>0){
            D = diff(diag(ncol(B)), diff = d)      
      }
      if(d==0){
            D = diag(ncol(B))      
      }
      
      A <- lapply(lambda,function(l){
             P = l * t(D) %*% D
             a <- solve(t(B) %*% B + P, t(B) %*% y)
             a <- as.vector(a)
             a
      })
      Bg <- bbase(xg,  xl = 0, xr = 1, nseg = 67, deg = 3)
       beta0hat <- lapply(A,function(a){
             Bg %*% a
       })
       Beta0 <- matrix(data=unlist(beta0hat),nrow=500,ncol=length(lambda),byrow = FALSE)
 
       fileName <- file.path(getwd(),"Dissertation TeX","img",paste0("PS_penalty_section_figure_6_order_",d ,".png"))
       png(filename = fileName)
       plot(x, y,pch="+",
            #main=expression(paste(hat(beta)[0](t)," with penalty ||",D[3] ,alpha,"|","|"^2)),
            xlab="",ylab="",
            ylim=c(min(Beta0[1:(nrow(Beta0)-10),])-0.04,max(Beta0+0.2))) 
       
       matlines(xg, Beta0, type = 'l', lty = 1, col=terrain.colors(nb,alpha=0.6),
                xlab = '', ylab = '')
       points(x, y,pch="+",col="grey")
       lines(xg,f(xg),col="red")
       dev.off()}
 






n <- 100  
x <- runif(n)
f <-  function(x){
      1.6 * x^3 * (3 * (1 - x))^2 + 2 * (x)^4
}
e <- rnorm(n, 0, 0.3)
y <- f(x) + e


B <- bbase(x,  xl = 0, xr = 1, nseg = 67, deg = 3)

for(d in 0:3){
      if(d>0){
            D = diff(diag(ncol(B)), diff = d)      
      }
      if(d==0){
            D = diag(ncol(B))      
      }
      
      A <- lapply(lambda,function(l){
            P = l * t(D) %*% D
            a <- solve(t(B) %*% B + P, t(B) %*% y)
            a <- as.vector(a)
            a
      })
      Bg <- bbase(xg,  xl = 0, xr = 1, nseg = 67, deg = 3)
      beta0hat <- lapply(A,function(a){
            Bg %*% a
      })
      Beta0 <- matrix(data=unlist(beta0hat),nrow=500,ncol=length(lambda),byrow = FALSE)
      
      fileName <- file.path(getwd(),"Dissertation TeX","img",paste0("PS_penalty_section_figure_6_order_",d ,".png"))
      png(filename = fileName)
      plot(x, y,pch="+",type="n",
           #main=expression(paste(hat(beta)[0](t)," with penalty ||",D[3] ,alpha,"|","|"^2)),
           xlab="",ylab="")#,
           #ylim=c(min(Beta0[1:(nrow(Beta0)-10),])-0.04,max(Beta0+0.2))) 
      points(seq(0,1,length.out=50),f(seq(0,1,length.out=50)),col="red",pch=19,cex=0.7)
      matlines(xg, Beta0, type = 'l', lty = 1, col=terrain.colors(nb,alpha=0.6),
               xlab = '', ylab = '')
      points(x, y,pch="+",col="grey")
#      lines(xg,f(xg),col="red")
      dev.off()}





















##=================================================================================================
##=================================================================================================
##=================================================================================================


############################################################################################
############################################################################################





f <- function(x) 2 * sin(pi * x*2)
x <- runif(100) %>% sort
e <- rnorm(length(x), 0, 0.2)
y <- f(x)*x + e


dat <- data.frame(y=y,x1=x,f=f(x),x2=x)
pDat <- data.frame(regressor=seq(0,1,length.out=200),
                   f=f(seq(0,1,length.out=200))) %>% transform(.,signal=regressor*f)
pDat_melt <- melt(pDat,id.vars = c("regressor"))
names(pDat_melt)[grepl("variable",names(pDat_melt))] <- "mean_component"

labls <- c(expression(paste(beta,"(t)")),
           expression(paste("t",beta,"(t)")))
p <- ggplot(data=dat,aes(x=x1,y=y)) + geom_point(shape="+") + ylab("") + xlab("")
p <- p + geom_line(data=pDat_melt,aes(x=regressor,y=value,group=mean_component,colour=mean_component)) + scale_color_ptol("",labels=labls) + theme_minimal()
p
ggsave(filename = file.path(getwd(),"Dissertation TeX","img","PS_VCM_section_figure_1.png"))


lambda <- as.list(exp(seq(-6,6,length.out=40)))
B <- bbase(x,  xl = 0, xr = 1, nseg = 76, deg = 4)
D = diff(diag(ncol(B)), diff = 3)
Q <- diag(x)%*%B

A <- lapply(lambda,function(l){
      P = l * t(D) %*% D
      a <- solve(t(Q) %*% Q + P, t(Q) %*% y)
      a <- as.vector(a)
      a
})
Bg <- bbase(xg,  xl = 0, xr = 1, nseg = 76, deg = 4)
beta1hat <- lapply(A,function(a){
      as.vector(Bg %*% a)
})
Beta1 <- matrix(data=unlist(beta1hat),nrow=nrow(Bg),ncol=length(beta1hat),byrow=FALSE)

png(filename = file.path(getwd(),"Dissertation TeX","img","PS_VCM_section_figure_2.png"))
plot(x, y,pch="+",ylim=c(-2.3,2.6),xlab="",type="n",
     #main=expression(paste(hat(beta)[1](t)," with penalty ||",D[3] ,alpha,"|","|"^2)),
     ylab="") 
matlines(xg, Beta1, type = 'l', lty = 1, col=terrain.colors(nb,alpha=0.6),
         xlab = '', ylab = '')
#points(x, y,pch="+")
lines(xg,f(xg),col="blue",lwd=1.5)
dev.off()

png(filename = file.path(getwd(),"Dissertation TeX","img","PS_VCM_section_figure_3.png"))
plot(x, y,pch="+",ylim=c(-1.8,.75),xlab="",type="n",
     #main=expression(paste(hat(beta)[1](t)," with penalty ||",D[3] ,alpha,"|","|"^2)),
     ylab="") 
matlines(xg, diag(xg)%*%Beta1, type = 'l', lty = 1, col=terrain.colors(nb,alpha=0.5),
         xlab = '', ylab = '')
points(x, y,pch="+",cex=0.7)
lines(xg,xg*f(xg),col="red",lwd=1.5)
dev.off()






B <- bbase(x,  xl = 0, xr = 1, nseg = 57, deg = 3)
D = diff(diag(ncol(B)), diff = 2)
Q <- diag(x)%*%B

A <- lapply(lambda,function(l){
      P = l * t(D) %*% D
      a <- solve(t(Q) %*% Q + P, t(Q) %*% y)
      a <- as.vector(a)
      a
})
Bg <- bbase(xg,  xl = 0, xr = 1, nseg = 57, deg = 3)
beta1hat <- lapply(A,function(a){
      as.vector(Bg %*% a)
})
Beta1 <- matrix(data=unlist(beta1hat),nrow=nrow(Bg),ncol=length(beta1hat),byrow=FALSE)

png(filename = file.path(getwd(),"Dissertation TeX","img","PS_penalty_section_figure_8.png"))
plot(x, y,pch="+",ylim=c(-2.4,2.4),xlab="",ylab="",
     #main=expression(paste(hat(beta)[1](t)," with penalty ||",D[3] ,alpha,"|","|"^2)),
     type="n") 
matlines(xg, Beta1, type = 'l', lty = 1, col=terrain.colors(nb,alpha=0.5),
         xlab = '', ylab = '')
#points(x, y,pch="+")
lines(xg,f(xg),col="blue",lwd=1.5)
dev.off()


png(filename = file.path(getwd(),"Dissertation TeX","img","PS_penalty_section_figure_9.png"))
plot(x, y,pch="+",ylim=c(-1.9,0.9),xlab="",ylab="",
     #main=expression(paste(hat(beta)[1](t)," with penalty ||",D[3] ,alpha,"|","|"^2)),
     type="n") 
matlines(xg, diag(xg)%*%Beta1, type = 'l', lty = 1, col=terrain.colors(nb,alpha=0.6),
         xlab = '', ylab = '')
#points(x, y,pch="+")
lines(xg,xg*f(xg),col="red")
dev.off()
















##=================================================================================================
##=================================================================================================
##=================================================================================================


############################################################################################
############################################################################################






f <- function(x){1.2 + 2*sin(5  * x)}
x <- runif(100) %>% sort
xg = seq(0, 1, length = 500)
e <- rnorm(length(x), 0, 0.2)
y <- f(x) + e


lambda <- as.list(c(10^(seq(-10,10,length.out=30))))
B <- bbase(x,  xl = 0, xr = 1, nseg = 67, deg = 3)

for(d in 0:3){
      if(d>0){
            D = diff(diag(ncol(B)), diff = d)      
      }
      if(d==0){
            D = diag(ncol(B))      
      }
      
      A <- lapply(lambda,function(l){
            P = l * t(D) %*% D
            M <- solve(t(B) %*% B + P)
            a <- M%*%t(B) %*% y
            a <- as.vector(a)
            H <- B%*%M%*%t(B)
            list(H=H,a=a)
      })
      Bg <- bbase(xg,  xl = 0, xr = 1, nseg = 67, deg = 3)
      beta0hat <- lapply(A,function(l){
            Bg %*% l$a
      })
      ED <- lapply(A,function(l){
            sum(diag(l$H))
      })
      LOOCV <- lapply(A,function(l){
            ((y-l$H%*%y)*(1/diag(l$H)))^2 %>% sum %>% sqrt
      })
      
      Beta0 <- matrix(data=unlist(beta0hat),nrow=500,ncol=length(lambda),byrow = FALSE)
      
      
      fileName <- file.path(getwd(),"Dissertation TeX","img",paste0("PS_penalty_section_figure_6_order_",d ,".png"))
      png(filename = fileName)
      plot(x, y,pch="+",
           #main=expression(paste(hat(beta)[0](t)," with penalty ||",D[3] ,alpha,"|","|"^2)),
           xlab="",ylab="",
           ylim=c(min(Beta0[1:(nrow(Beta0)-10),])-0.04,max(Beta0+0.2))) 
      
      matlines(xg, Beta0, type = 'l', lty = 1, col=terrain.colors(nb,alpha=0.6),
               xlab = '', ylab = '')
      points(x, y,pch="+",col="grey")
      lines(xg,f(xg),col="red")
      dev.off()}

      plot(log(unlist(lambda)),unlist(ED),type="l")







##################################################################################################
##################################################################################################

##################################################################################################
##################################################################################################
      












p1 <- p2 <- 100
M1.index <- M2.index <- seq(0,1,length.out=100)

B1 <- bsplbase(M1.index, c(0,1,60,3))
B2 <- bsplbase(M2.index, c(0,1,60,3))
n1 <- ncol(B1)
n2 <- ncol(B2)	# Compute tensor products for estimated alpha surface
Gamma1 <-( 0.5*cos(2*pi*seq(0,1,length.out=(ncol(B1))))) %o% (0.5*sin(2*pi*seq(0,1,length.out=(ncol(B1)))))
Gamma2 <- (0.9*((1-seq(0,1,length.out=(ncol(B1)))))^3)%o% ((2*(seq(0,1,length.out=(ncol(B1))))^2)+(1.7*((1-seq(0,1,length.out=(ncol(B1)))))^3))
eSmoothRow <- matrix(data=loess.smooth(x=seq(0,1,length.out = ncol(B1)^2),
                                       y=rnorm(ncol(B1)^2,sd=1),span=0.13,
                                       degree=1,
                                       evaluation = ncol(B1)^2)$y,
                     nrow=ncol(B1),
                     ncol=ncol(B1),
                     byrow = TRUE)
eSmoothCol <- matrix(data=loess.smooth(x=seq(0,1,length.out = ncol(B1)^2),
                                       y=rnorm(ncol(B1)^2,sd=1),span=0.13,
                                       degree=1,
                                       evaluation = ncol(B1)^2)$y,
                     nrow=ncol(B1),
                     ncol=ncol(B1),
                     byrow = FALSE)
Y <-  B1%*%(Gamma1+Gamma2+eSmoothRow+eSmoothCol)%*%t(B2)



png(filename = file.path(getwd(),"Dissertation TeX","img",paste0("2D_penalty_section_figure_1",".png")))
persp(M2.index, M1.index, z=Y,
      xlab = "", ylab = "", zlab = "",theta = 45, phi = 15,
      shade = 0.65,col="lightblue",
      box=TRUE,d=6,
      ticktype = "detailed",
      nticks=5,lty=2,border = NA,
      cex.axis=0.5,ltheta=120)
dev.off()

      list_of_alpha <- list()
      list_of_hat_mat <- list()


      y <- as.vector(Y)
      B1 <- bsplbase(M1.index, c(0,1,40,3))
      B2 <- bsplbase(M2.index, c(0,1,40,3))
      lambda <- 1000
      lambdaBar <- 10000000
      
      d <- 0
      #D <- diff(diag(ncol(B1)), diff = d)
      D <- diag(ncol(B1))
      P <- kronecker(diag(rep(1,ncol(B2))), t(D)%*%D )
      dBar <- 1
      Dbar <- diff(diag(ncol(B2)), diff = dBar)
      Pbar <- kronecker(t(Dbar)%*%Dbar, diag(rep(1,ncol(B1))) )

oM1 <- outer(rep(1, length(M2.index)),M1.index)
B1 <- bsplbase(as.vector(oM1), c(0,1,40,3))
oM2 <- outer(M2.index, rep(1, length(M1.index)))
B2 <- bsplbase(as.vector(oM2), c(0,1,40,3))
n1 <- ncol(B1)
n2 <- ncol(B2)	# Compute tensor products for estimated alpha surface
B1. <- kronecker(B1, t(rep(1, n2)))
B2. <- kronecker(t(rep(1, n1)), B2)
rm(B1)
rm(B2)
rm(oM1)
rm(oM2)
B. <- B1. * B2.

M <- solve((t(B.) %*% B. + lambda*P + lambdaBar*Pbar))
rm(P)
rm(Pbar)

a <-  M%*%t(B.) %*% as.vector(y)
H <- B.%*%M%*%t(B.)
rm(B.)
rm(M)
ED <- H %>% diag %>% sum %>% round(.,2)


list_of_alpha[[length(list_of_alpha)+1]] <- a
list_of_hat_mat[[length(list_of_hat_mat)+1]] <- H
muHat <- B1%*%matrix(data=a,nrow=43,ncol=43,byrow=FALSE)%*%t(B2)


fileName <- file.path(getwd(),"Dissertation TeX","img",paste0("2D_penalty_section_figure_2_colOrder",dBar,"_rowOrder_",d,"_4.png"))
png(file=fileName)
persp(M2.index, M1.index, z=muHat,
      xlab = "", ylab = "", zlab = "",theta = 40, phi = 15,
      shade = 0.65,col="lightblue",
      box=TRUE,d=6,
      ticktype = "detailed",
      nticks=5,lty=2,border = NA,
      cex.axis=0.5,ltheta=120,
      main=bquote(paste("ED(",lambda,") = ",.(ED))) )
dev.off()













##=================================================================================================
##=================================================================================================
##=================================================================================================


############################################################################################
############################################################################################



############################################################################################

source(file.path(getwd(),"code","fnc","pspline_signal_fit_1d.R"))  
library(ggthemes)  
n <- 50  
x1 <- seq(0,1,length.out=n)
f <-  function(x){
      0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10 
}
e <- rnorm(n, 0, 0.3)
# A continuous `by' variable example.... 
y <- f(x1)*x1 + e
dat <- data.frame(y=y,x1=x1,f=f(x1),x2=x1)
pDat <- data.frame(regressor=seq(0,1,length.out=200),
                   f=f(seq(0,1,length.out=200))) %>% transform(.,signal=regressor*f)
pDat_melt <- melt(pDat,id.vars = c("regressor"))
names(pDat_melt)[grepl("variable",names(pDat_melt))] <- "mean_component"

p <- ggplot(data=dat,aes(x=x1,y=y)) + geom_point() + ylab("") + xlab("x")
p <- p + geom_line(data=pDat_melt,aes(x=regressor,y=value,group=mean_component,colour=mean_component)) + scale_color_ptol("",labels=c("f(x)","xf(x)")) + theme_minimal()
p








n.seg <- 50
x.l <- min(dat$x1)
x.r <- max(dat$x1)
x.max <- x.r + 0.01 * (x.r - x.l)
x.min <- x.l - 0.01 * (x.r - x.l)
d.x <- (x.max - x.min)/n.seg
ps.knots <- seq(x.min - 2 * d.x, x.max + 2 * d.x, by = d.x)


ps.signal.fit <- gam(y ~ s(x1,
                           bs="ps",
                           by=x2,
                           k=(length(ps.knots)-6),
                           m=c(4,4),
                           sp=100,
                           id=1),
                     data=dat[-c(1:3,(nrow(dat)-2):nrow(dat)),],
                     knots = list(x1=ps.knots))


df <- data.frame(x=dat$x1,
           y=dat$y,
           fitted=ps.signal.fit$fitted.values)
ggplot(melt(df,id.vars="x"),aes(x,value)) + geom_point(aes(colour=variable)) + theme_calc() + scale_color_calc()
ggplot(df,aes(x=x,y=y)) + geom_point(colour=calc_pal()(2)[1]) + theme_calc() + geom_line(aes(x=x,y=fitted),colour=calc_pal()(2)[2])

B <- spline.des(ps.signal.fit$smooth[[1]]$knots, seq(0,1,length.out=200), ps.signal.fit$smooth[[1]]$p.order[1]+1)$design
fitted.f <- B%*%ps.signal.fit$coefficients
df <- data.frame(x=seq(0,1,length.out = 200),
                 true.f=f(seq(0,1,length.out = 200)),
                 fitted.f=fitted.f) %>% melt(.,id.vars="x")
ggplot(df,aes(x=x,y=value,group=variable)) + geom_line(aes(colour=variable)) + theme_calc() + scale_color_calc()





ll <- lapply(as.list(exp(seq(-4,5,length.out = 20))),function(l){
                  ps.signal.fit <- gam(y ~ s(x1,bs="ps", by=x2,
                                             sp=l,
                                             k=(length(ps.knots)-(2*3)),
                                             m=c(4,3)),
                                       data=dat[-c(1:3,(nrow(dat)-2):nrow(dat)),],
                                       knots = list(x1=ps.knots))
            
                  llist <- list()
                  llist$sp <- l
                  llist$fitted.values <- ps.signal.fit$fitted.values
                  llist$predicted.grid <- predict.gam(ps.signal.fit,
                                                newdata = data.frame(x1=seq(0,1,length.out=nrow(B)),x2=seq(0,1,length.out=nrow(B))))
                  llist$coefficients <- ps.signal.fit$coefficients
                  llist$fitted.f <- B%*%ps.signal.fit$coefficients 
                  return(llist)
                  })

f.df <- ldply(ll,function(l){return(data.frame(lambda=l$sp,
                                               x=seq(0,1,length.out=length(l$predicted.grid)),
                                               f.hat=as.vector(l$fitted.f)))})
ind <- expand.grid(x=1:length(ll[[1]]$predicted.grid),ind=1:length(ll))$ind
f.df <- cbind.data.frame(f.df,ind)

p <- ggplot(f.df,aes(x=x,y=f.hat)) + geom_line(aes(colour=ind,group=lambda)) + scale_colour_gradient_tableau("Green") + guides(colour="none")
p + geom_line(data=data.frame(x=seq(0,1,length.out=200),true.f=f(seq(0,1,length.out=200)),lambda=0),aes(x=x,y=true.f,group=lambda),colour="red",inherit.aes = FALSE)































