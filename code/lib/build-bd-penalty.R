

build_bd_pen <- function (lambda_list,
                          basis_list,
                          diff_ord_list) {
  
  
  
  D1 <- diff(diag(ncol(B1)),
             differences=dl)
  D2 <- diff(diag(ncol(Bm)),
             differences=dm)
  D3 <- diff(diag(ncol(B3)),
             differences=dl)
  D4 <- diff(diag(ncol(B4)),
             differences=dl)
  D5 <- diff(diag(ncol(B5)),
             differences=dm)
  D6 <- diff(diag(ncol(B6)),
             differences=dm)
  dl <- 2
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
  D7 <- Dlm_l
  
  
  dm <- 2
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
  D8 <- Dlm_m
  Diff_list <- list()
    Diff_list[[1]] <- D1
    Diff_list[[2]] <- D2
    Diff_list[[3]] <- D3
    Diff_list[[4]] <- D4
    Diff_list[[5]] <- D5
    Diff_list[[6]] <- D6
    Diff_list[[7]] <- list(D7=D7,D8=D8)
    
    
  lam_list <- as.list(c(lambda_l,lambda_m,lambda_l, lambda_l, lambda_m, lambda_m))
  lam_list[[7]] <- list(lam_l=lambda_l,
                        lam_m=lambda_m)

  
  if (length(lambda_list) != length(basis_list)) {
    stop( "the number of penalty parameters must match the number of functional components" )
  }

  l <- list.zip(lambda=lam_list,D=Diff_list)
  Pen_list <- list()
  for (i in 1:length(l)) {
    if(is.list(l[[i]]$D)){
      Pen_list[[i]] <-  list.rbind(lapply(list.zip(lam=l[[i]]$lambda,D=l[[i]]$D),
                                          function(ll) {
                                            ll$lam*ll$D
                                          }))
    }
    if(!is.list(l[[i]]$D)) {
     Pen_list[[i]] <-  l[[i]]$lambda*l[[i]]$D    
    }
  }  
  
  ncol_blkdiag_Pen <- lapply(B_list,ncol) %>% unlist %>% sum %>% add(1)
  col_to_skip <- c(0,cumsum(unlist(lapply(Pen_list,ncol))))
  
  big_Pen <- NULL
  for(blk in 1:length(Pen_list)) {
    if (!is.list(Pen_list[[i]])){
      big_Pen <- rbind(big_Pen,
                       cbind(matrix(data=0,
                                    nrow=nrow(Pen_list[[i]]),
                                    ncol=ncol_blkdiag_Pen-col_to_skip[i]-ncol(Pen_list[[i]])),
                             Pen_list[[i]],
                             matrix(data=0,
                                    nrow=nrow(Pen_list[[i]]),
                                    ncol=ncol_blkdiag_Pen-col_to_skip[i]-ncol(Pen_list[[i]]))))  
    }
  }
  big_Pen <- cbind(rep(0,nrow(big_Pen)),big_Pen)

  
    Pen <- rbind(## P1 ####################################################
               cbind(rep(0,nrow(Dl)),
                     lambda_pair$lam_l1*Dl,
                     matrix(data=0,nrow=nrow(Dl),
                            ncol=((2*ncol(Dm))+ncol(Dl)))),
               ## P2 ####################################################
               cbind(rep(0,nrow(Dl)),
                     matrix(data=0,nrow=nrow(Dm),
                            ncol=ncol(Dl)),
                     lambda_pair$lam_m1*Dm,
                     matrix(data=0,nrow=nrow(Dm),ncol=ncol(Dl)+ncol(Dm))),
               ## P3 ####################################################
               cbind(rep(0,nrow(Dl)),
                     matrix(data=0,nrow=nrow(Dm),
                            ncol=ncol(Dl)+ncol(Dm)),
                     lambda_pair$lam_l2*Dl,
                     matrix(data=0,nrow=nrow(Dm),ncol=ncol(Dm))),
               ## P4 ####################################################
               cbind(rep(0,nrow(Dl)),
                     matrix(data=0,nrow=nrow(Dm),
                            ncol=((2*ncol(Dl))+ncol(Dm))),
                     lambda_pair$lam_m2*Dm))
  
}
