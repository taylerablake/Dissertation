
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




bPars <- rbind(c(0,1,35,3,100,3),
               c(0,1,35,3,100,3))


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
knot_grid <- transform(knot_grid,
                       m_index=m_index-min(m_index)+1)

B. <- NULL
for(this_pair in 1:nrow(knot_grid)) {
      B. <- cbind(B.,
                  Bl[,match(knot_grid$l[this_pair],knots.l)]*Bm[,match(knot_grid$m[this_pair],knots.m)])
}

nrow_Dlm_m <-sum(table(knot_grid$l_index)[table(knot_grid$l_index)>dm]-dm)
ls_to_fix <- as.numeric(names(table(knot_grid$l_index)[table(knot_grid$l_index)>dm]))


dl <- 3
if(dl>0){
  if(sum(knot_grid$m_index==1)>dl){
    Dlm_l <- matrix(data=0,nrow=sum(knot_grid$m_index==1)-dl,
                    ncol=nrow(knot_grid))
    Dlm_l[,which(knot_grid$m_index==1)] <- diff(diag(sum(knot_grid$m_index==1)),
                                                differences = dl)            
    for(i in 2:length(unique(knot_grid$m_index))){
      d_l <- matrix(data=0,nrow=sum(knot_grid$m_index==i)-dl,
                    ncol=nrow(knot_grid))
      d_l[,which(knot_grid$m_index==i)] <- diff(diag(sum(knot_grid$m_index==i)),
                                                differences = dl)      
      Dlm_l <- rbind(Dlm_l,d_l)
    }
  }
  
  if((sum(knot_grid$m_index==1)<=dl) &(sum(knot_grid$m_index==2)>dl)){
    Dlm_l <- matrix(data=0,nrow=sum(knot_grid$m_index==2)-dl,
                    ncol=nrow(knot_grid))
    Dlm_l[,which(knot_grid$m_index==2)] <- diff(diag(sum(knot_grid$m_index==2)),
                                                differences = dl)            
    for(i in 3:length(unique(knot_grid$m_index))){
      if(sum(knot_grid$m_index==i)>dl){
        d_l <- matrix(data=0,nrow=sum(knot_grid$m_index==i)-dl,
                      ncol=nrow(knot_grid))
        d_l[,which(knot_grid$m_index==i)] <- diff(diag(sum(knot_grid$m_index==i)),
                                                  differences = dl)      
        Dlm_l <- rbind(Dlm_l,d_l)      
      }
    }
  }
}


dm <- 3
if(dm>0){
  if(sum(knot_grid$l_index==1)>dm){
    Dlm_m <- matrix(data=0,nrow=sum(knot_grid$l_index==1)-dm,
                    ncol=nrow(knot_grid))
    Dlm_m[,which(knot_grid$l_index==1)] <- diff(diag(sum(knot_grid$l_index==1)),
                                                differences = dm)            
    for(i in 2:max(knot_grid$l_index)){
      if(sum(knot_grid$l_index==i)>dm){
        d_m <- matrix(data=0,nrow=sum(knot_grid$l_index==i)-dm,
                      ncol=nrow(knot_grid))
        d_m[,which(knot_grid$l_index==i)] <- diff(diag(sum(knot_grid$l_index==i)),
                                                  differences = dm)      
        Dlm_m <- rbind(Dlm_m,d_m)
      }
    }
  }
  
  if((sum(knot_grid$l_index==1)<=dl) &(sum(knot_grid$l_index==2)>dm)){
    Dlm_m <- matrix(data=0,nrow=sum(knot_grid$l_index==2)-dm,
                    ncol=nrow(knot_grid))
    Dlm_m[,which(knot_grid$l_index==2)] <- diff(diag(sum(knot_grid$l_index==2)),
                                                differences = dm)            
    for(i in 3:length(unique(knot_grid$l_index))){
      if(sum(knot_grid$l_index==i)>dm){
        d_m <- matrix(data=0,nrow=sum(knot_grid$l_index==i)-dm,
                      ncol=nrow(knot_grid))
        d_m[,which(knot_grid$l_index==i)] <- diff(diag(sum(knot_grid$l_index==i)),
                                                  differences = dm)      
        Dlm_m <- rbind(Dlm_m,d_m)      
      }
    }
  }
}



B1 <- Bl 
B2 <- Bm
B3 <- grid$m/max(grid$m)*Bl
B4 <- ((grid$m/max(grid$m))^2)*Bl
B5 <- grid$l/max(grid$l)*Bm
B6 <- ((grid$l/max(grid$l))^2)*Bm
B7 <- B.

B_list <- list()
B_list[[1]] <- B1
B_list[[2]] <- B2
B_list[[3]] <- B3
B_list[[4]] <- B4
B_list[[5]] <- B5
B_list[[6]] <- B6
B_list[[7]] <- B7

# B <- cbind(as.vector(rep(1,nrow(Bl))),
#            B1,B2,B3,B4,B5,B6,B7)
B <- cbind(as.vector(rep(1,nrow(Bl))),
           Bl,Bm,B.)

U <- X %*% B

B_list <- list()
B_list[[1]] <- Bl
B_list[[2]] <- Bm
B_list[[3]] <- B.
nrow_Pen <- nrow(D1)+nrow(D2)+nrow(Dlm_l)+nrow(Dlm_m)


D1 <- diff(diag(ncol(Bl)),
           differences=dl)
D2 <- diff(diag(ncol(Bm)),
           differences=dm)
lambda_l <- 8
lambda_m <- 8
lambda_l2 <- 8
lambda_m2 <- 8
Pen <- NULL
Pen <- rbind(Pen,cbind(lambda_l*D1,
                       matrix(data=0,
                              nrow=nrow(D1),
                              ncol=ncol(B)-ncol(D1)-1)))
Pen <- rbind(Pen,cbind(matrix(data=0,
                              nrow=nrow(D1),
                              ncol=ncol(D1)),
             lambda_m*D2,
             matrix(data=0,
                    nrow=nrow(D1),
                    ncol=ncol(B)-1-ncol(D2)-ncol(D1))))
Pen <- rbind(Pen,
             cbind(matrix(data=0,nrow=(nrow(Dlm_l)+nrow(Dlm_m)),
                          ncol=(ncol(D1)+ncol(D2))),
                   rbind(lambda_l2*Dlm_l,
                         lambda_m2*Dlm_m)))
Pen <- cbind(rep(0,nrow(Pen)),Pen)
Pen <- rbind(Pen,
             ridge.adj*cbind(rep(0,ncol(B)-1),diag(ncol(B)-1)))

nix <- as.vector(rep(0,nrow(Pen)))
w <- c(as.vector(1/diag(D)))
y_vec <-  c(as.vector(t(y[,-1]))) 

fit <- lsfit(rbind(U,Pen),
             c(y_vec, nix),
             wt = c(w, nix + 1) *c(w, nix + 1),
             intercept = F)

phi <- diag(M)
phi[lower.tri(phi)] <- B %*% fit$coefficients
my_title <- c(as.expression(bquote(lambda[l] == .(lambda_l))),
              as.expression(bquote(lambda[m] == .(lambda_m))))

wireframe(phi,
          scales = list(arrows = FALSE),
          xlab="",
          ylab="",
          zlab=expression(phi),
          main=my_title,
          zlim=c(-0.3,1),
          screen = list(x = 90, y = -70))


wireframe(t(phi),
          scales = list(arrows = FALSE),
          xlab="",
          ylab="",
          zlab=expression(phi),
          main=my_title,
          zlim=c(-0.3,1))





Bl <- bsplbase(as.vector(grid$l)/max(grid$l),
               bPars[1,],outer.okay = TRUE)$base
Bm <- bsplbase(grid$m/max(grid$m),
               bPars[2,],outer.okay = TRUE)$base

