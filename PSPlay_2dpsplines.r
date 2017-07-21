



# Simulate data
set.seed(1985)
n = 300
x1 = runif(n)
x2 = runif(n)

y = 1.2 + sin(5  * x1) * cos(2*pi*x2) + rnorm(n) * 0.2
order = 2
lla = 1

# Preparations
x1g = seq(0, 1, length = 500)
x2g = seq(0, 1, length = 500)
nseg1 = 8
nseg2 = 8
bdeg = 3
pord1 = 2
pord2 = 2




# row-tensor product (or Box product of two matrices)
rowtens <- function(X, B) {
  one.1 <- matrix(1, 1, ncol(X))
  one.2 <- matrix(1, 1, ncol(B))
  kronecker(X, one.2) * kronecker(one.1, B)
}




drawit = function() {
  #nb = ncol(B)
  #cols = hcl(h = seq(60, 240, length = nb), c =90, l = 70)
  B1g <<- bbase(x1g, xl = 0, xr = 1, nseg = nseg1, deg = bdeg)
  B2g <<- bbase(x2g, xl = 0, xr = 1, nseg = nseg2, deg = bdeg)
  Bg <<- rowtens(B2g,B1g)
  #A = diag(a)
  z = Bg %*% a
  plot(x, y, col="grey",pch="+") 
  matlines(xg, Bg %*% A, type = 'l', lty = 1, lwd = 2, col= cols,
           xlab = '', ylab = '', ylim = c(0, 1))
  lines(xg, z, col = 'red', lwd = 3)        
  knots = seq(0, 1, length = nseg + 1)
  points(knots, 0 * knots, pch = 15, cex = 0.8)
  tl =paste('P-splines, n = ', nb, ', order = ', pord, 
            ', degree = ', bdeg, ', log10(lambda) = ', lla, sep = '')
  title(tl)
}


ps.smooth = function(p){
  nseg <<- floor(p$nseg)
  bdeg <<- p$bdeg
  lla <<- p$lla
  lambda <<- 10 ^ lla
  pord <<- p$ord 
  #  cat(nseg, bdeg, pord, '\n')
  B <<- bbase(x,  xl = 0, xr = 1, nseg = nseg, deg = bdeg)
  nb = ncol(B)
  D = diff(diag(nb), diff = pord)
  P = lambda * t(D) %*% D
  a <<- solve(t(B) %*% B + P, t(B) %*% y)
  a <<- as.vector(a)
  drawit()
  return(p)
}


Bg <<- bbase(xg, xl = 0, xr = 1, nseg = nseg, deg = bdeg)
A = diag(a)
z = Bg %*% a
plot(x, y, col="grey",pch="+") 
matlines(xg, Bg %*% A, type = 'l', lty = 1, lwd = 2, col= cols,
         xlab = '', ylab = '', ylim = c(0, 1))
lines(xg, z, col = 'red', lwd = 3)        
knots = seq(0, 1, length = nseg + 1)
points(knots, 0 * knots, pch = 15, cex = 0.8)
tl =paste('P-splines, n = ', nb, ', order = ', pord, 
          ', degree = ', bdeg, ', log10(lambda) = ', lla, sep = '')
title(tl)

B1 <- bsplbase(as.vector(oM1), Pars[1,  ])
oM2 <- outer(M2.index, rep(1, p1))
B2 <- bsplbase(as.vector(oM2), Pars[2,  ])
n1 <- ncol(B1)
n2 <- ncol(B2)	# Compute tensor products for estimated alpha surface
B1. <- kronecker(B1, t(rep(1, n2)))
B2. <- kronecker(t(rep(1, n1)), B2)
B. <- B1. * B2.


