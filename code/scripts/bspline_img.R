


setwd(file.path("/Users","taylerblake","Documents","Dissertation"))
source(file.path(getwd(),"code","fnc","aux","bsplbase.R"))

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


bSpline <- spline.des(knots=seq(0,1,by=.05), x=seq(0,1,length.out = 200), ord = 4,outer.ok = TRUE)
bS <- data.frame(expand.grid(x=seq(0,1,length.out = 200),j=1:ncol(bSpline$design)),
                 B=matrix(bSpline$design,nrow=200*ncol(bSpline$design),ncol=1),
                 q=3)
      bS <- merge(bS,data.frame(j=1:ncol(bSpline$design),knot=paste0("x",1:ncol(bSpline$design))),by="j")

p <- ggplot(subset(bS,((j == 2 | j %in% (13:16))&q==3)),
            aes(x=x,y=B,group=j)) + ylab("") + xlab("") 
p <- p + scale_x_continuous(breaks = bSpline$knots[c(2:6,13:20)],
                            labels = paste0("x",c(c(1:5,12:19))))
p <- p + theme(axis.text.y = element_blank(),
               axis.text.x = element_text(angle = 45,vjust=0.1),
               axis.ticks.y = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"))
p + geom_line(aes(x=x,y=B,group=j)) 
ggsave(filename = file.path(getwd(),"Dissertation TeX","img","uni_cubic_bsplines.png"),
       width = 7.25,height = 4,units = "in")





bSpline <- spline.des(knots=seq(0,1,by=.05), x=seq(0,1,length.out = 200), ord = 2,outer.ok = TRUE)
bS <- data.frame(expand.grid(x=seq(0,1,length.out = 200),
                                      j=1:ncol(bSpline$design)),
                          B=matrix(bSpline$design,nrow=200*ncol(bSpline$design),ncol=1),
                          q=1)

p <- ggplot(subset(bS,  j %in% c(3,14:17)&q==1 ),
            aes(x=x,y=B,group=j)) + ylab("") + xlab("") 
p <- p + scale_x_continuous(breaks = bSpline$knots[c(3:5,14:19)],
                            labels = paste0("x",c(c(1:3,14:19))))
p <- p + theme(axis.text.y = element_blank(),
               axis.text.x = element_text(angle = 45,vjust=0.1),
               axis.ticks.y = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"))
p + geom_line(aes(x=x,y=B,group=j))# + facet_wrap(~ q,nrow=2,labeller = label_both)
ggsave(filename = file.path(getwd(),"Dissertation TeX","img","uni_linear_bsplines.png"),
       width = 7.25,height = 4,units = "in")











############################################################################################
############################################################################################








# Simulate data
n = 100
x = seq(0, 1, length = n)
set.seed(123)
y = 1.2 + sin(5  * x) + rnorm(n) * 0.2
order = 2
lla = 1

# Preparations
xg = seq(0, 1, length = 500)
nseg = 20
bdeg = 3
pord = 2









############################################################################################
############################################################################################


p1 <- p2 <- 200
M1.index <- M2.index <- seq(0,1,length.out=200)

oM1 <- outer(rep(1, p2),M1.index)
B1 <- bsplbase(as.vector(oM1), c(0,1,4,3))


oM2 <- outer(M2.index, rep(1, p1))
B2 <- bsplbase(as.vector(oM2), c(0,1,4,3))
n1 <- ncol(B1)
n2 <- ncol(B2)	# Compute tensor products for estimated alpha surface
B1. <- kronecker(B1, t(rep(1, n2)))
B2. <- kronecker(t(rep(1, n1)), B2)
B. <- B1. * B2.
dim(B.)

# Create a function interpolating colors in the range of specified colors

jet.colors <- colorRampPalette( c("deepskyblue2","green") )
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
z <- matrix(B.[,25],nrow=length(M1.index),ncol=length(M1.index))

nrz <- nrow(z)
ncz <- ncol(z)

zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]


# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)
par(bg="white")
b <- bsplbase(seq(0,1,length.out=200),c(0,1,4,3))[,4]

png(filename = file.path(getwd(),"Dissertation TeX","img","bicubic_basis_function.png"))
persp(M2.index, M1.index, z,
      xlab = "", ylab = "", zlab = "",theta = 45, phi = 9,
      shade = 0.55,col=color[facetcol],
      box=TRUE,d=6,
      ticktype = "detailed",
      nticks=5,lty=2,border = NA,
      cex.axis=0.5,ltheta=120,zlim=c(0,max(b+0.02)))
persp(M2.index, M1.index, z,
      xlab = "", ylab = "", zlab = "",theta = 45, phi = 9,
      shade = 0.55,col=color[facetcol],
      box=TRUE,d=6,
      ticktype = "detailed",
      nticks=5,lty=2,border = NA,
      cex.axis=0.5,ltheta=120,zlim=c(0,max(b+0.02)))-> res

for (i in seq(0.1,0.6,by=0.1)){
      lines (trans3d(x=M1.index, y = 1, z = i, pmat = res),lty=3,col="gray84")
      lines (trans3d(x=0, y = M1.index, z = i, pmat = res),lty=3,col="gray84")
}
for(i in seq(0.2,1,by=0.2)){
      lines (trans3d(x=0, y = i, z = seq(0,max(b)+0.02,by=0.01), pmat = res),lty=3,col="gray84")
      lines (trans3d(x=i, y = 1, z = seq(0,max(b)+0.02,by=0.01), pmat = res),lty=3,col="gray84")
}
lines(trans3d(x=0, y = M2.index, z = b, pmat = res))
lines(trans3d(x=M1.index, y = 1, z = b, pmat = res))
#dev.off()

############################################################################################
############################################################################################
png(filename = file.path(getwd(),"Dissertation TeX","img","bicubic_bspline_contour.png"))
contour(x=M1.index, y=M2.index, z, nlev = 8, method="edge", xlim=c(0,1),ylim=c(0,1))
for(i in seq(0.2,0.8,by=0.2)){
      abline(v=i,lty=3,col="lightgray")
      abline(h=i,lty=3,col="lightgray")
}
dev.off()
############################################################################################
############################################################################################
## make a 3d surface plot of a handful of cubic b-spline functions

B <- bsplbase(seq(0,1,by=0.01),bpars=c(0,1,12,3)) 
plot(seq(0,1,by=0.01),B[,12],xlab="x",ylab="B(x)",type="l")
B. <- data.frame(x=rep(seq(0,1,by=0.01),ncol(B)),b=matrix(B,nrow=ncol(B)*nrow(B),byrow = FALSE),knot=expand.grid(seq(0,1,by=0.01),n=1:ncol(B))$n)
B. %>% ggplot(.,aes(x=x,y=b,group=knot)) + geom_line(aes(colour=factor(knot))) + guides(colour="none")


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
#z <- matrix(rowSums(B.[,ind]),nrow=length(M1.index),ncol=length(M1.index))
#z <- rowSums(B.[,ind])
#nrz <- nrow(z)
#ncz <- ncol(z)

#zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices

#M <- data.frame(expand.grid(x=M1.index,tilde.x=M2.index),z=rowSums(B.[,ind]))
M <- data.frame(expand.grid(x=M1.index,tilde.x=M2.index),z=rowSums(B.[,ind]))
newcols <- colorRampPalette(c("grey70", "grey10"))

par(bg="white")
png(filename = file.path(getwd(),"Dissertation TeX","img","sparse_bicubic_basis.png"))
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
dev.off()






gamma <- rep(0,ncol(B.))
gamma[ind] <- c(1,1,1,
                1,2,3,
                1,2,3)
gamma.mat <- matrix(data=rep(gamma,nrow(B.)),nrow=nrow(B.),ncol=ncol(B.),byrow=TRUE)


#M <- data.frame(expand.grid(x=M1.index,tilde.x=M2.index),z=rowSums(B.[,ind]))
M <- data.frame(expand.grid(x=M1.index,tilde.x=M2.index),z=rowSums(gamma.mat*B.))

newcols <- colorRampPalette(c("grey50", "grey20"))
png(filename = file.path(getwd(),"Dissertation TeX","img","large_row_penalty.png"))
wireframe(z~x*tilde.x,data=M,
          lty=3,
          lwd=0.2,
          ylab="x",
          xlab=expression(tilde(x)),
          zlab="",
          screen=list(z=30,x=-55),
          #aspect = c(68/87, 0.9),
          light.source = c(10,0,10),
          pretty=TRUE,
          #scales = list(arrows = FALSE),
          col.regions = newcols(100),
          drape=TRUE,
          colorkey=FALSE,
          par.settings = list(axis.line = list(col = "transparent")))
dev.off()


gamma[ind] <- c(1,3,1,
                1,1,3,
                2,2,1)
gamma.mat <- matrix(data=rep(gamma,nrow(B.)),nrow=nrow(B.),ncol=ncol(B.),byrow=TRUE)
M <- data.frame(expand.grid(x=M1.index,tilde.x=M2.index),z=rowSums(gamma.mat*B.))
png(filename = file.path(getwd(),"Dissertation TeX","img","large_col_penalty.png"))
wireframe(z~x*tilde.x,data=M,
          lty=3,
          lwd=0.2,
          ylab="x",
          xlab=expression(tilde(x)),
          zlab="",
          screen=list(z=30,x=-55),
          #aspect = c(68/87, 0.9),
          light.source = c(10,0,10),
          pretty=TRUE,
          #scales = list(arrows = FALSE),
          col.regions = newcols(100),
          drape=TRUE,
          colorkey=FALSE,
          par.settings = list(axis.line = list(col = "transparent")))
dev.off()
