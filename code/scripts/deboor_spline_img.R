
betaDat <- data.frame(x=rbeta(1000,.7,2))
F_n <- ecdf(betaDat$x)
y_breaks <- c(F_n(0.1),diff(F_n(seq(0.1,1,by=0.1))))*1000
x_breaks <- seq(0.1,1,by=0.1)
density(betaDat$x,bw=0.05,from=0,to=1)%>%plot

H <- betaDat %>% ggplot(.,aes(x)) + geom_histogram(aes(y=..density..),
                                                   binwidth=.1,
                                                   center=0.05,
                                                   fill="white",
                                                   color="grey")
H <- H + theme_minimal() + xlim(-0.3,1.3) + scale_x_continuous(breaks=x_breaks[-c(5:8)],
                                                         labels=c(expression(tau[1]),
                                                                  expression(tau[2]),
                                                                  expression(tau[3]),
                                                                  expression(tau[4]),
                                                                  # expression(tau[5]),
                                                                  # expression(tau[6]),
                                                                  # expression(tau[7]),
                                                                  # expression(tau[8]),
                                                                  expression(tau[n]),
                                                                  expression(tau[n+1])))

H <- H + scale_y_continuous(breaks=y_breaks[-10],labels=c(expression(h[1]),
                                                expression(h[2]),
                                                expression(h[3]),
                                                expression(h[4]),
                                                expression(h[5]),rep(fixed(""),4)))
H <- H + xlab("") + ylab("")
H + geom_line(data=data.frame(x=c(0,seq(0.01,1,length.out=100)),
                               y=c(3.58,dbeta(seq(0.01,0.98,length.out=100)+0.02,0.7,2))),
              aes(x,y),inherit.aes = FALSE) 
ggsave(file=file.path(getwd(),"Dissertation TeX","img","histogram_smoothing.png"))





#############################################################################################


tpf <- function(x,breakpoint,degree){
      max(0,(x-breakpoint)^degree)      
}

4.65-(0.00001/2)
4.65+(0.00001/2)
xi <- c(seq(0.1,3.1,by=1),4.5,4.8,seq(6.1,9.1,by=1),10)
diff(xi)
alpha <- c(1.3,0,0,0,
           1.1/diff(xi)[4],
           -2.2/diff(xi)[5] - 1.1/diff(xi)[4],
           2.2/diff(xi)[5] + 1.1/diff(xi)[6],
           -1.1/diff(xi)[6],
           0,0,0,0)
alpha %>% round(.,3)
broken_line <- function(x,Xi,Alpha,p){
      
      Alpha[1] + sum(Alpha[-1]*sapply(Xi,function(xi_arg){tpf(x,xi_arg,p)}))
}

sapply(xi,function(x_arg){broken_line(x_arg,xi,round(alpha,4),1)})
broken_line(9.6,xi,round(alpha,5),1)


