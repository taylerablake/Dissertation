


drawit = function() {
      nb = ncol(B)
      cols = hcl(h = seq(60, 240, length = nb), c =90, l = 70)
      Bg <<- bbase(xg, xl = 0, xr = 1, nseg = nseg, deg = bdeg)
      A = diag(a)
      z = Bg %*% a
      plot(x, y, lty = 1, type = 'l') 
      matlines(xg, Bg %*% A, type = 'l', lty = 1, lwd = 2, col= cols,
               xlab = '', ylab = '', ylim = c(0, 1))
      lines(xg, z, col = 'red', lwd = 3)        
      knots = seq(0, 1, length = nseg + 1)
      points(knots, 0 * knots, pch = 15, cex = 0.8)
      tl =paste('P-splines, n = ', nb, ', order = ', pord, 
                ', degree = ', bdeg, ', log10(lambda) = ', lla, sep = '')
      title(tl)
}
