

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
library(stringr)


for (fileName in list.files(file.path(getwd(),"code","fnc","aux"))){
     source(file.path(getwd(),"code","fnc","aux",fileName)) 
}

for (fileName in list.files(file.path(getwd(),"code","fnc"))){
      if(str_detect(fileName,".r") || str_detect(fileName,".R")){
            source(file.path(getwd(),"code","fnc",fileName))      
      }
}




bspline<- function(x,xl,xr,ndx,bdeg){
      dx<- (xr-xl) /ndx
      knots<- seq(xl- bdeg*dx,xr+bdeg*dx,by=dx) 
      B<- spline.des(knots,x,bdeg+1,0*x)$design
      B
}









#################################################################################################
#################################################################################################
## Demonstrate the most basic mechanics of the 2d penalty, showing what a 
## strong row penalty, weak column penalty
## weak row penalty, strong column penalty
## strong row penalty, strong column penalty
## weak row penalty, weak column penalty
##looks like on a sparse set of tensor basis function




p1 <- p2 <- 100
M1.index <- M2.index <- seq(0,1,length.out=100)

B1 <- bsplbase(M1.index, c(0,1,11,3))
B2 <- bsplbase(M2.index, c(0,1,11,3))
ind <- as.matrix(expand.grid(seq(4,10,by=3),seq(4,10,by=3)))

alpha <- matrix(data=0,nrow=ncol(B1),ncol=ncol(B1))
alpha[ind] <- c( 1,    1.5,    1,
                1.8,    1,    0.5,
                 1,     2,     3)

B1 <- bsplbase(M1.index, c(0,1,10,4))
B2 <- bsplbase(M2.index, c(0,1,10,4))
ind <- as.matrix(expand.grid(seq(4,10,by=3),seq(4,10,by=3)))

alpha <- matrix(data=0,nrow=ncol(B1),ncol=ncol(B1))

##############################################################
## strong column penalty


alpha[ind] <- c( 0.5,   1.8,    0.8,
                 1,    1.2,      0.8,
                 1.5,     0.6,     0.8)

fileName <- file.path(getwd(),"Dissertation TeX","img",
                      "model selection",
                      "effective dimension",
                      "2d_PS_section_figure1_big_col_lambda.png")
png(filename = fileName)

persp(M2.index, M1.index, z=B1%*%alpha%*%t(B2),
      xlab = "", ylab = "", zlab = "",theta = 150, phi = 25,
      shade = 0.65,col="lightblue",
      box=TRUE,d=2,
      axes=FALSE,
      lty=2,border = NA,
      cex.axis=0.5,ltheta=40)
dev.off()

alpha[ind] <- c( 1,     .6,     .2,
                 0.2,    0.55,      .9,
                 0.4,   0.3,    0.2)


fileName <- file.path(getwd(),"Dissertation TeX","img",
                      "model selection",
                      "effective dimension",
                      "2d_PS_section_figure1_big_row_lambda.png")
png(filename = fileName)
persp(M2.index, M1.index, z=B1%*%alpha%*%t(B2),
      xlab = "", ylab = "", zlab = "",theta = 150, phi = 25,
      shade = 0.65,col="lightblue",
      box=TRUE,d=2,
      axes=FALSE,
      lty=2,border = NA,
      cex.axis=0.5,ltheta=40)
dev.off()








y <- as.vector(B1%*%alpha%*%t(B2))


startTime <- Sys.time()
psSmooth <- psp2dG(Data=cbind(as.matrix(expand.grid(M1.index,M2.index)[,c(2,1)]),y),
       Pars=matrix(data=c(0,1,30,3,100,2, 0,1,30,3,1000,2),
                   nrow=2,ncol=6,byrow=TRUE),
       ridge.adj = 0.001,
       x.lab = "X",
       y.lab = "Y",
       z.lab = "Z",
       persp.plot = F,
       image.plot = T,
       se = F, 
       family = "gaussian",
       link = "default")
endTime <- Sys.time()
endTime-startTime
names(psSmooth)


persp(M2.index, M1.index,
      z=matrix(data=psSmooth$mu,nrow=100,ncol=100,byrow=FALSE),
      xlab = "", ylab = "", zlab = "",theta = 320, phi = 30,
      shade = 0.65,col="lightblue",
      box=TRUE,d=2,
      axes=FALSE,
      lty=2,border = NA,
      cex.axis=0.5,ltheta=40)






#################################################################################################
#################################################################################################




p1 <- p2 <- 60
M1.index <- M2.index <- seq(0,1,length.out= p1)

B1 <- bsplbase(M1.index, c(0,1,50,3))
B2 <- bsplbase(M2.index, c(0,1,50,3))
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
B1 <- bsplbase(M1.index, c(0,1,50,3))
B2 <- bsplbase(M2.index, c(0,1,50,3))
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
B1 <- bsplbase(as.vector(oM1), c(0,1,6,3))
oM2 <- outer(M2.index, rep(1, length(M1.index)))
B2 <- bsplbase(as.vector(oM2), c(0,1,6,3))
n1 <- ncol(B1)
n2 <- ncol(B2)

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



