
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





f <- function(x){1.2 + 2*sin(5  * x)}
x <- runif(10)
xg = seq(0, 1, length = 500)
e <- rnorm(length(x), 0, 0.6)
y <- f(x) + e



lambda <- as.list(c(10^(seq(-7,10,length.out=50))))
B <- bbase(x,  xl = 0, xr = 1, nseg = 67, deg = 3)

Beta0_list <- list()
ED_list <- list()
LOOCV_list <- list()


for(d in 0:4){
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
            ((y-l$H%*%y)*(1/(1-diag(l$H))))^2 %>% sum %>% sqrt
      })
      
      Beta0 <- matrix(data=unlist(beta0hat),nrow=500,ncol=length(lambda),byrow = FALSE)
      Beta0_list[[d+1]] <- Beta0
      
      fileName <- file.path(getwd(),"Dissertation TeX","img","model selection",paste0("PS_ED_section_figure_1_order_",d ,".png"))
      png(filename = fileName)
      plot(x, y,pch="+",
           #main=expression(paste(hat(beta)[0](t)," with penalty ||",D[3] ,alpha,"|","|"^2)),
           xlab="",ylab="")#,
           #ylim=c(min(Beta0[1:(nrow(Beta0)-10),])-0.04,max(Beta0+0.2))) 
      
      matlines(xg, Beta0, type = 'l', lty = 1, col=terrain.colors(nb,alpha=0.6),
               xlab = '', ylab = '')
      points(x, y,pch="+",col="grey")
      lines(xg,f(xg),col="red")
      dev.off()

      ED_list[[d+1]] <- unlist(ED)
      #LOOCV_list[[d+1]] <- unlist(LOOCV)
      
      }


ED <- data.frame(ed=as.vector(unlist(ED_list)),
                 expand.grid(l=unlist(lambda),d=(0:4)))

p <- ggplot(ED,aes(x=log(l),y=ed,group=d)) + geom_line(aes(colour=factor(d)))
p <- p + guides(colour=guide_legend("d")) 
p <- p + xlab(expression(paste(log(lambda)))) + ylab("ED")
p +theme_minimal() + scale_color_wsj()
ggsave(file.path(getwd(),"Dissertation TeX","img","model selection","PS_ED_section_figure_1.png"))





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


f <- function(x){1.2 + 2*sin(5  * x)}
x <- runif(100)
xg = seq(0, 1, length = 500)
e <- rnorm(length(x), 0, 0.6)
y <- f(x) + e



lambda <- as.list(c(10^(seq(-8,10,length.out=50))))
B <- bbase(x,  xl = 0, xr = 1, nseg = 67, deg = 3)

Beta0_list <- list()
ED_list <- list()
LOOCV_list <- list()


for(d in 0:4){
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
            ((y-l$H%*%y)*(1/(1-diag(l$H))))^2 %>% sum %>% sqrt
      })
      
      Beta0 <- matrix(data=unlist(beta0hat),nrow=500,ncol=length(lambda),byrow = FALSE)
      Beta0_list[[d+1]] <- Beta0
      
      fileName <- file.path(getwd(),"Dissertation TeX","img","model selection",paste0("PS_LOOCV_section_figure_1_order_",d ,".png"))
      png(filename = fileName)
      plot(x, y,pch="+",
           #main=expression(paste(hat(beta)[0](t)," with penalty ||",D[3] ,alpha,"|","|"^2)),
           xlab="",ylab="")#,
      #ylim=c(min(Beta0[1:(nrow(Beta0)-10),])-0.04,max(Beta0+0.2))) 
      
      matlines(xg, Beta0, type = 'l', lty = 1, col=terrain.colors(nb,alpha=0.6),
               xlab = '', ylab = '')
      points(x, y,pch="+",col="grey")
      lines(xg,f(xg),col="red")
      dev.off()
      
      LOOCV_list[[d+1]] <- unlist(LOOCV)
      
}

LOOCV <- lapply(A,function(l){
      ((y-l$H%*%y)*(1/diag(l$H)))^2 %>% sum %>% sqrt
})
LOOCV_list[[d+1]] <- unlist(LOOCV)



LOOCV <- data.frame(LOOCV=as.vector(unlist(LOOCV_list)),
                 expand.grid(l=unlist(lambda),d=(0:4)))

p <- ggplot(subset(LOOCV,d>0),
            aes(x=log(l),y=LOOCV,group=d)) + geom_line() + facet_wrap(~ d,nrow=2,scales="free")
#p <- p + guides(colour=guide_legend("d")) 
p <- p + xlab(expression(paste(log(lambda)))) + ylab("LOOCV")
p +theme_minimal() + scale_color_wsj()
ggsave(file.path(getwd(),"Dissertation TeX","img","model selection","PS_LOOCV_section_figure_1.png"))


ddply(LOOCV,.(d),summarize,
      minCV=min(LOOCV),best_lambda=l[which.min(LOOCV)])

# LOOCV <- LOOCV_list %>% unlist %>% matrix(data=.,nrow=length(lambda),ncol=4,byrow=FALSE)
# plot(log(unlist(lambda)),unlist(LOOCV),type="l",
#      xlab=expression(paste(log(lambda))),
#      ylab=expression(paste(CV(lambda))))
