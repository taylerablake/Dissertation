
source(file.path(getwd(),"lib","build-grid.R"))


M <- 20
N <- 100
grid <- build_grid(M)

Sigma <- matrix(0.7,nrow=m,ncol=m) + diag(rep(0.3),m)
Omega <- solve(Sigma)

C <- t(chol(Sigma))
D <- diag(diag(C))
L <- C%*%solve(D)
T_mod <- solve(L)
D <- diag(rep(diag(C)[-1],N))

grid <- transform(grid,l_s=l/max(grid$l),
                  m_s=m/max(grid$m))






y <- mvrnorm(n=N,mu=rep(0,M),Sigma=Sigma)
y_vec <- as.vector(t(y[,-1]))


X <- matrix(data=0,nrow=N*(M-1),ncol=choose(M,2))
no.skip <- 0
for (t in 2:M){
      X[((0:(N-1))*(M-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- y[,1:(t-1)]
      no.skip <- no.skip + t - 1
}




bPars <- rbind(c(0,1,30,3,100,3),
               c(0,1,30,3,100,3))


Bl <- bsplbase(as.vector(grid$l)/max(grid$l),
               bPars[1,],outer.okay = TRUE)$base
Bm <- bsplbase(grid$m/max(grid$m),
               bPars[2,],outer.okay = TRUE)$base

knots.l <- bsplbase(as.vector(grid$l)/max(grid$l),
                    bPars[1,],outer.okay = TRUE)$knots[1:ncol(Bl)]
knots.m <-bsplbase(as.vector(grid$m)/max(grid$m),
                   bPars[2,],outer.okay = TRUE)$knots[1:ncol(Bm)]


n1 <- ncol(Bl)
n2 <- ncol(Bm)
knot_grid <- data.frame(expand.grid(m=knots.m,l=knots.l),
                        expand.grid(m_index=1:length(knots.m),
                                    l_index=1:length(knots.l)))
knot_grid <- subset(knot_grid, (m > 0.5*l) & (m < 1-(0.5*l))) 
knot_grid <- knot_grid[,c(2,1,4,3)]


B. <- NULL
for(this_pair in 1:nrow(knot_grid)) {
      B. <- cbind(B.,
                  Bl[,match(knot_grid$l[this_pair],knots.l)]*Bm[,match(knot_grid$m[this_pair],knots.m)])
}

nrow_Dlm_m <-sum(table(knot_grid$l_index)[table(knot_grid$l_index)>dm]-dm)
ls_to_fix <- as.numeric(names(table(knot_grid$l_index)[table(knot_grid$l_index)>dm]))


D_list <- list()
for(this_knot in 1:length(ls_to_fix)) {
      D_list[[this_knot]] <- diff(diag(sum(knot_grid$l_index==ls_to_fix[this_knot])),
                                  differences = dm)
}

Dlm_m <- as.list((lapply(D_list,nrow) %>% unlist %>% cumsum %>% c(0,.))[-(length(D_list)+1)]) %>%
      list.zip(before_rows=.,
               D=D_list) %>%
      lapply(.,function(l) {
            rbind(matrix(data=0,
                         nrow=l$before_rows,
                         ncol=ncol(l$D)),
                  l$D,
                  matrix(data=0,
                         nrow=nrow_Dlm_m-(l$before_rows+nrow(l$D)),
                         ncol=ncol(l$D)))
      }) %>%
      list.cbind %>%
      cbind(.,matrix(as.vector(rep(0,
                                   nrow(Dlm_m)*sum(table(knot_grid$l_index)<=dm))),
                     ncol=sum(table(knot_grid$l_index)<=dm),
                     nrow=nrow(Dlm_m)))



B1 <- Bl 
B2 <- Bm
B3 <- grid$m/max(grid$m)*Bl
B4 <- ((grid$m/max(grid$m))^2)*Bl
B5 <- grid$l/max(grid$l)*Bm
B6 <- ((grid$l/max(grid$l))^2)*Bm
B7 <- B.


B <- cbind(as.vector(rep(1,nrow(Bl))),
           B1,B2,B3,B4,B5,B6,B7)
U <- X %*% B





