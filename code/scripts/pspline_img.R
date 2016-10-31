
setwd(file.path("/Users","taylerblake","Documents","Dissertation"))
source(file.path(getwd(),"code","fnc","aux","bsplbase.R"))
source(file.path(getwd(),"code","fnc","draw_psplines.R"))

library(plyr)
library(dplyr)
library(rlist)
library(ggplot2)
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




############################################################################################
############################################################################################





f <- function(x) 2 * sin(pi * x*2)
x <- runif(100) %>% sort
e <- rnorm(length(x), 0, 0.2)
y <- f(x) + e

ggplot(data = data.frame(x=x,y=y),aes(x,y)) + geom_point()

  ll <- c(seq(-3,16,length.out=150))
  fit_list <- lapply(as.list(ll),function(l){
                  fit <- pspline.fit( y, x, ps.intervals = 100, 
                                      degree = 3, order = 3,
                                      ridge.adj = 0.0001,lambda=exp(l))
                  fit
              })

  
  fits <- data.frame(x=rep(seq(0,1,length.out=100),length(fit_list)),
             yhat=lapply(fit_list,function(l){spline.des(l$knots,seq(l$knots[2*l$order+1],l$knots[length(l$knots)-(2*l$order+1)],length.out=100), l$degree + 1)$design%*%l$coef}) %>% list.rbind(),
             lambda=expand.grid(1:100,exp(ll))[,2],
             index=expand.grid(1:100,1:length(ll))[,2]) 
  
  diagnostics <- data.frame(aic=lapply(fit_list,function(l){l$aic})%>% unlist,
                            ed=lapply(fit_list,function(l){l$eff.df})%>% unlist,
                            lambda=lapply(fit_list,function(l){l$lambda})%>% unlist,
                            index=1:length(fit_list))
  ggplot(diagnostics,aes(x=log(lambda),y=aic)) + geom_line(colour="blue") + xlab(expression(log(lambda)))        
  ggplot(diagnostics,aes(x=log(lambda),y=ed)) + geom_line(colour="blue") + xlab(expression(log(lambda))) + ylab("tr(H)")
  
  F.n <- ecdf(diagnostics$aic)
  ind <- which(F.n(diagnostics$aic) %in% seq(0.05,1,by=0.05))
  ind <- c(ind,sapply(seq(0.05,1,by=0.05)[!(seq(0.05,1,by=0.05)%in%F.n(diagnostics$aic))],
                      function(q){which.min(abs(F.n(diagnostics$aic)-q))}))
  
  
  p <-  ggplot(subset(fits,lambda%in%diagnostics$lambda[ind]),aes(x=x,y=yhat,group=lambda)) + geom_line(aes(colour=index)) + scale_color_gradientn(colours = heat.colors(16))
  p <- p + geom_point(data=data.frame(x=x,y=y),aes(x,y),colour="grey") + guides(colour="none")
  #true.f <- data.frame(x=seq(min(x),max(x),length.out=200),f=f(seq(min(x),max(x),length.out=200)))
  
  p + geom_line(data=subset(fits,lambda%in%diagnostics$lambda[which.min(diagnostics$aic)]),aes(x=x,y=yhat),size=1.2,colour="blue")


############################################################################################

Y_conts_by_var <- gamSim(eg=3,n=100)
ggplot(data=Y_conts_by_var,aes(x=x1,y=y)) + geom_point()
ggplot(data=Y_conts_by_var,aes(x=x2,y=f)) + geom_point()
ggplot(data=Y_conts_by_var,aes(x=1:nrow(Y_conts_by_var),y=x1*f)) + geom_point()


x1 <- sort(runif(100, 0, 1))
f <-  0.2 * x1^11 * (10 * (1 - x1))^6 + 
      10 * (10 * x1)^3 * (1 - x1)^10 
e <- rnorm(n, 0, 1)
# A continuous `by' variable example.... 
y <- f*x1 + e
Y_conts_by_var <- data.frame(y=y,x1=x1,f=f)

p <- ggplot(data=Y_conts_by_var,aes(x=x1,y=f)) + geom_line(colour="pink") + ylab("") + xlab("")
p <- p + geom_line(aes(x=x1,y=x1*f),colour="blue")
p <- p +  geom_point(aes(x=x1,y=y),colour="grey")
p



