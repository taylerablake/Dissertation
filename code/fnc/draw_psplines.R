

source(file.path(getwd(),"code","Marx Pspline Course","Computer_labs",'bases.r'))

draw_psplines = function(nseg,bdeg,pord,lla,nobs) {

      
      x = seq(0, 1, length = nobs)
      xg = seq(0, 1, length = 500)
      set.seed(123)
      y = 1.2 + sin(5  * x) + rnorm(nobs) * 0.2
      
      
      B <- bbase(x,  xl = 0, xr = 1, nseg = nseg, deg = bdeg)
      nb = ncol(B)
      D = diff(diag(nb), diff = pord)
      
      lambda <<- 10 ^ lla
      P = lambda * t(D) %*% D
      a <- solve(t(B) %*% B + P, t(B) %*% y)
      a <- as.vector(a)
      
      
      cols = hcl(h = seq(60, 240, length = nb), c =90, l = 70)
      Bg <<- bbase(xg, xl = 0, xr = 1, nseg = nseg, deg = bdeg)
      A = diag(a)
      z = Bg %*% a
      plot(x, y,ylim=c(0,2.5),pch="+",col="grey55",xlab="",ylab="") 
      matlines(xg, Bg %*% A, type = 'l', lty = 1, lwd = 2, col= cols,
               xlab = '', ylab = '')
      lines(xg, z, col = 'red')        
      knots = seq(0, 1, length = nseg+1)
      points(knots, 0 * knots, pch = 15, cex = 0.8)
      points(knots, a[-c(1,length(a))], pch=19, col="red", cex=0.8)
      tl =paste('n_basis = ', nb, ', penalty order = ', pord, 
                ', degree = ', bdeg, ', log10(lambda) = ', lla, sep = '')
      #title(tl)
}


