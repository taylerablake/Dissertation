setwd(file.path("/Users","taylerblake","Documents","Dissertation"))
source(file.path(getwd(),"code","fnc","aux","bsplbase.R"))
source(file.path(getwd(),"code","Marx Pspline Course",'bases.r'))


library(plyr)
library(dplyr)
library(rlist)
library(ggplot2)
library(tidyr)
library(splines)
library(reshape2)
library(systemfit)
require(graphics)
library(latex2exp)
library(Matrix)

rowDiffPenaltyProjection <- function(row_diff_order,n_knots_row,n_knots_col){
      kronecker((t(diff(diag(n_knots_row), diff = row_diff_order))%*%diff(diag(n_knots_row), diff = row_diff_order)),
                diag(n_knots_col))
}

colDiffPenaltyProjection <- function(col_diff_order,n_knots_row,n_knots_col){
      kronecker(diag(n_knots_row),
                (t(diff(diag(n_knots_col), diff = col_diff_order))%*%diff(diag(n_knots_col), diff = col_diff_order)))
}

p1 <- p2 <- 200
M1.index <- M2.index <- seq(0,1,length.out=200)

oM1 <- outer(rep(1, p2),M1.index)
B1 <- bsplbase(as.vector(oM1), c(0,1,10,2))
oM2 <- outer(M2.index, rep(1, p1))
B2 <- bsplbase(as.vector(oM1), c(0,1,10,2))
n1 <- ncol(B1)
n2 <- ncol(B2)	# Compute tensor products for estimated alpha surface
B1. <- kronecker(B1, t(rep(1, n2)))
B2. <- kronecker(t(rep(1, n1)), B2)
B. <- B1. * B2.
ind <- c(43,43+3,43+6,ncol(B.)-43+1,ncol(B.)-(43+3)+1,ncol(B.)-(43+6)+1,82,85,88)
gamma <- rep(0,ncol(B.))
gamma[ind] <- c(2,8,4,
                8,4,2,
                8,2,4)
gamma[ind] <- 1 
gamma.mat <- matrix(data=rep(gamma,nrow(B.)),ncol=ncol(B.),nrow=nrow(B.),byrow = TRUE)

#P1 <- rowDiffPenaltyProjection(row_diff_order = 1,n_knots_row = n1,n_knots_col = n2)
#P2 <- colDiffPenaltyProjection(col_diff_order = 1,n_knots_row = n1,n_knots_col = n2)
P1 <- rowDiffPenaltyProjection(row_diff_order = 1,n_knots_row = sqrt(ncol(M)),n_knots_col = sqrt(ncol(M)))
P2 <- colDiffPenaltyProjection(col_diff_order = 1,n_knots_row = sqrt(ncol(M)),n_knots_col = sqrt(ncol(M)))

lambda1 <- 0.00001      
lambda2 <- 6
y <- M%*% gamma
y <- rowSums(M)
gamma.hat <- solve(t(M)%*%M + lambda1*t(P1)%*%P1 + lambda2*t(P2)%*%P2)%*%t(M)%*% y
y.hat <- M%*% gamma.hat

G <- data.frame(expand.grid(x=M1.index,tilde.x=M2.index),z=rowSums(B.))
G <- data.frame(expand.grid(x=M1.index,tilde.x=M2.index),z=y.hat)

z <- matrix(rowSums(gamma.mat*B.),nrow=length(M1.index),ncol=length(M1.index))
z <- matrix(rowSums(B.[,ind]),nrow=length(M1.index),ncol=length(M1.index))
persp(M2.index, M1.index,z,
      xlab = "", ylab = "", zlab = "",theta = 45, phi = 30,
      shade = 0.55,#col=color[facetcol],
      box=TRUE,d=6,
      ticktype = "detailed",
      nticks=5,lty=2,border = NA,
      cex.axis=0.5,ltheta=120)

wireframe(z~x*tilde.x,data=G,
          lty=3,
          lwd=0.3,
          ylab="x",
          xlab=expression(tilde(x)),
          zlab="",
          screen=list(z=44,x=-60),
          aspect = c(68/87, 0.8),
          light.source = c(10,0,10),
          #pretty=TRUE,
          scales = list(arrows = FALSE),
          drape=TRUE,
          colorkey=FALSE,
          par.settings = list(axis.line = list(col = "transparent")))































p1 <- p2 <- 200
M1.index <- M2.index <- seq(0,1,length.out=200)

oM1 <- outer(rep(1, p2),M1.index)
B1 <- bsplbase(as.vector(oM1), c(0,1,10,3))
oM2 <- outer(M2.index, rep(1, p1))
B2 <- bsplbase(as.vector(oM2), c(0,1,10,3))
n1 <- ncol(B1)
n2 <- ncol(B2)	# Compute tensor products for estimated alpha surface
B1. <- kronecker(B1, t(rep(1, n2)))
B2. <- kronecker(t(rep(1, n1)), B2)
B. <- B1. * B2.


# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("lightpink","dodgerblue") )
# Generate the desired number of colors from this palette
nbcol <- 200
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
ind <- c(43,43+3,43+6,ncol(B.)-43+1,ncol(B.)-(43+3)+1,ncol(B.)-(43+6)+1,82,85,88)
gamma.mat
#z <- matrix(rowSums(B.[,ind]),nrow=length(M1.index),ncol=length(M1.index))
#z <- rowSums(B.[,ind])
#nrz <- nrow(z)
#ncz <- ncol(z)

#zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices

#M <- data.frame(expand.grid(x=M1.index,tilde.x=M2.index),z=rowSums(B.[,ind]))
M <- data.frame(expand.grid(x=M1.index,tilde.x=M2.index),z=rowSums(B.[,ind]))
newcols <- colorRampPalette(c("grey70", "grey10"))
wireframe(z~x*tilde.x,data=M,
          lty=3,
          lwd=0.3,
          ylab="x",
          xlab=expression(tilde(x)),
          zlab="",
          screen=list(z=44,x=-65),
          aspect = c(68/87, 0.8),
          light.source = c(10,0,10),
          #pretty=TRUE,
          scales = list(arrows = FALSE),
          col.regions = newcols(100),
          drape=TRUE,
          colorkey=FALSE,
          par.settings = list(axis.line = list(col = "transparent")))




gamma <- rep(0,ncol(B.))
gamma[ind] <- 1 
gamma.mat <- matrix(data=rep(gamma,nrow(B.)),ncol=ncol(B.),nrow=nrow(B.),byrow = TRUE)
M <- data.frame(expand.grid(x=M1.index,tilde.x=M2.index),z=rowSums(gamma.mat*B.))
wireframe(z~x*tilde.x,data=M,
          lty=3,
          lwd=0.3,
          ylab="x",
          xlab=expression(tilde(x)),
          zlab="",
          screen=list(z=44,x=-65),
          aspect = c(68/87, 0.8),
          light.source = c(10,0,10),
          #pretty=TRUE,
          scales = list(arrows = FALSE),
          col.regions = newcols(100),
          drape=TRUE,
          colorkey=FALSE,
          par.settings = list(axis.line = list(col = "transparent")))

gamma[ind] <- 1:9
gamma[ind] <- c(2,2,2,2,4,6,2,4,6)/6
gamma.mat <- matrix(data=rep(gamma,nrow(B.)),ncol=ncol(B.),nrow=nrow(B.),byrow = TRUE)
M <- data.frame(expand.grid(x=M1.index,tilde.x=M2.index),z=rowSums(gamma.mat*B.))
wireframe(z~x*tilde.x,data=M,
          lty=3,
          lwd=0.3,
          ylab="x",
          xlab=expression(tilde(x)),
          zlab="",
          screen=list(z=44,x=-65),
          aspect = c(68/87, 0.8),
          light.source = c(10,0,10),
          #pretty=TRUE,
          scales = list(arrows = FALSE),
          col.regions = newcols(100),
          drape=TRUE,
          colorkey=FALSE,
          par.settings = list(axis.line = list(col = "transparent")))
