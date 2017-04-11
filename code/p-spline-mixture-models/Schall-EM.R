
# Choose true function to simulate
if (Type. == "f01") {
    R1 <- sin(2 * pi * (x1 - 0.5))
    R2 <- cos(2 * pi * (x2 - 0.5))
    y0 <- R1 + R2
} else if (Type. == "f02") {
    R <- ((x1 - 0.5)^2 + (x2 - 0.5)^2)
    y0 <- cos(2 * pi * sqrt(R))
} else if (Type. == "f03") {
    R1 <- sin(2 * pi * x1)
    R2 <- cos(2 * pi * x2)
    R12 <- sin(2 * pi * (x2 - x1))
    y0 <- R1 + R2 + R12
} else if (Type. == "f04") {
    y0 <- 2 * x1 * sin(4 * pi * x2)
} else {
    stop("Model not specified")
}
y <- y0 + rnorm(n) * sigma.


start1 <- proc.time()[3]

# Build Mixed Model Bases For nested bases for the interaction terms, create
# new mixed model basis with number of segments as a integer divisor of nseg1,
# and nseg2 for div=1, G1=G1n and G2=G2n

MM1 <- MM.basis(x1, min(x1) - 0.01, max(x1) + 0.01, nseg1, bdeg, pord)  # bases for x1
X1 <- MM1$X
G1 <- MM1$Z
d1 <- MM1$d
B1 <- MM1$B

MM2 <- MM.basis(x2, min(x2) - 0.01, max(x2) + 0.01, nseg2, bdeg, pord)  # bases for x2
X2 <- MM2$X
G2 <- MM2$Z
d2 <- MM2$d
B2 <- MM2$B

MM1n <- MM.basis(x1, min(x1) - 0.01, max(x1) + 0.01, nseg1/div, bdeg, pord)  # Nested bases for x1
G1n <- MM1n$Z
d1n <- MM1n$d
B1n <- MM1n$B

MM2n <- MM.basis(x2, min(x2) - 0.01, max(x2) + 0.01, nseg2/div, bdeg, pord)  # Nested bases for x2
G2n <- MM2n$Z
d2n <- MM2n$d
B2n <- MM2n$B

c1 <- ncol(B1)
c2 <- ncol(B2)

c1n <- ncol(B1n)
c2n <- ncol(B2n)

one1. <- matrix(X1[, 1], ncol = 1)
one2. <- matrix(X2[, 1], ncol = 1)

x1. <- matrix(X1[, 2], ncol = 1)
x2. <- matrix(X2[, 2], ncol = 1)

#####################
X <- rowtens(X2, X1)  # -> Fixed effects
#
d3 <- c(rep(1, c2n - pord) %x% d1n + d2n %x% rep(1, c1n - pord))
Delta1 <- diag(1/sqrt(d1))
Delta2 <- diag(1/sqrt(d2))
Delta3 <- diag(1/sqrt(d3))

# random effects matrices
Z1 <- G1 %*% Delta1  # smooth random comp. fx1
Z2 <- G2 %*% Delta2  # smooth random comp. fx2
Z1x2 <- rowtens(x2., G1n)  # linear:smooth comp. x2:fx1
Z2x1 <- rowtens(G2n, x1.)  # linear:smooth comp. fx2:x1
Z12 <- rowtens(G2n, G1n) %*% Delta3  # smooth interaction  fx1:fx2
#####################

# Random effects matrix
Z <- cbind(Z1, Z2, Z1x2, Z2x1, Z12)

# Compute cross-products

XtX. <- crossprod(X)
XtZ. <- crossprod(X, Z)
ZtZ. <- crossprod(Z)
Xty. <- crossprod(X, y)
Zty. <- crossprod(Z, y)
yty <- sum(y^2)
####################################################################################################

V <- construct.block2(XtX., XtZ., ZtZ.)
u <- c(Xty., Zty.)

# Administrative helpers to locate sub-bases

np <- c(ncol(X), ncol(Z1), ncol(Z2), ncol(Z1x2), ncol(Z2x1), ncol(Z12))

nblocks <- 6  # number of components (1 fixed component + 5 random effects)

idx <- rep(1:nblocks, np)
end1 <- proc.time()[3]
nn <- sum(np)

# Initialize algorithm
la <- c(0, c(rep(1, nblocks)))

# Set lower threshold on variances
thr <- 1e-08
start2 <- proc.time()[3]
# Run Schall's EM algorithm
if (flag.) {
    cat(" It  D-lambdas   fixed     fx1      fx2    fx1:x2    x1:fx2   fx1:fx2\n")
}
for (it in 1:200) {
    
    # Build penalty matrix
    P <- diag(la[idx])
    
    # Solve equations
    Vu <- cbind(V, u)
    Hb <- solve(V + P, Vu)
    HH <- Hb[, 1:nn]
    b <- Hb[, nn + 1]
    
    # Compute effective dimensions and variances
    ed <- tapply(diag(HH), idx, sum)
    ssb <- tapply(b^2, idx, sum)
    tau <- ssb/ed + thr
    ssr <- yty - crossprod(b, 2 * u - V %*% b)
    sig2 <- (ssr/(n - sum(ed)))
    
    # New lambdas and convergence check
    lanew <- c(0, rep(sig2, nblocks - 1))/tau
    dla <- mean(abs(log10(la[2:nblocks]) - log10(lanew[2:nblocks])))
    la <- lanew
    if (flag.) {
        cat(sprintf("%1$3d %2$10.6f", it, dla))
        cat(sprintf("%8.3f", ed), "\n")
    }
    if (dla < 1e-06) 
        break
}
end2 <- proc.time()[3]

M <- cbind(X, Z)
Fit <- M %*% b


## Select spatials points for prediction
ngrid <- 50  # generate a n x n grid....
mx1 <- seq(min(x1), max(x1), length = ngrid)
mx2 <- seq(min(x2), max(x2), length = ngrid)
gx1 <- rep(mx1, ngrid)
gx2 <- rep(mx2, rep(ngrid, ngrid))

M.Grid <- Grid2d(gx1, gx2, nseg1, nseg2, div, bdeg, pord, ngrid, b)

Fit.grid <- matrix(M.Grid %*% b, ngrid, ngrid)

index1 <- 1 * (idx == 2)
index1[2] <- 1
Fit1 <- matrix(M.Grid %*% (index1 * b), ngrid, ngrid)

index2 <- 1 * (idx == 3)
index2[3] <- 1
Fit2 <- matrix(M.Grid %*% (index2 * b), ngrid, ngrid)

index12a <- 1 * (idx == 4)
index12a[4] <- 1
Fit12a <- matrix(M.Grid %*% (index12a * b), ngrid, ngrid)

index12b <- 1 * (idx == 5)
Fit12b <- matrix(M.Grid %*% (index12b * b), ngrid, ngrid)

index12c <- 1 * (idx == 6)
index12c[4] <- 1
Fit12c <- matrix(M.Grid %*% (index12c * b), ngrid, ngrid)

y.hat <- c(Fit)
RSS <- sum((y - Fit)^2)
AIC. <- n * log(RSS/n) + 2 * sum(ed)
MSE. <- mean((y0 - Fit)^2)

cat("----------- PS-ANOVA-Schall -------------\n")
cat("Sample size           ", n, "\n")
cat("sigma                 ", sigma., "\n")
cat("sigma (estimated)     ", round(sig2^0.5, 3), "\n")
if (div == 1) {
    cat("No nested bases \n")
} else {
    cat("Nested bases for interactions, with nseg1=", nseg1, ", nseg2=", nseg2, "and div=", 
        div, "\n")
}
cat("------------------------------------------\n")
cat("Timings:\nconstructing matrices", end1 - start1, "seconds\n")
cat("Estimation           ", (end2 - start2), "seconds\n")
cat("-----------------------------------------\n")
cat("   AIC", AIC., "\n")
cat("          fixed    fx1       fx2     fx1:x2     x1:fx2    fx1:fx2\n")
cat("    ed", sprintf("%8.3f", ed), "\n")
cat("sum ed", sum(ed), "\n")
cat("MSE", MSE., "\n")
cat("-----------------------------------------\n")
