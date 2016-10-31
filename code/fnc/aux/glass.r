"glass"<-
      function(x, ps.intervals = NULL, lambda = 0, x.signal = NULL, signal.index = 
                     NULL, varying.index = NULL, degree = 3, order = 3, ridge.adj = 1e-005, 
               ridge.inv = 0.0001)
      {
            # NOTE: see Instructions and Examples
            # x: regressor of interest (also see signal.index and varying.index below) 
            # ps.intervals: number of equally spaced B-spline intervals 
            #	(the number of knots is equal to ps.int+2*degree+1
            # lambda: non-negative regularization parameter for difference penalty
            # x.signal: if SIGNAL Component, this is the n X p signal matrix (p >> n possible)
            # signal.index: if SIGNAL Component, this is the axis where B-splines are 
            #	constructed (specify above x as 1:ncol(x.signal))
            # varying.index: if VARYING Component, this is the indexing variable where B-splines 
            #	are constructed (e.g. time)
            # degree: degree of B-spline basis
            # order: order of difference penalty (0 is the ridge penalty)
            # ridge.adj, ridge.inv: small positive numbers to stabilize 
            #	linear dependencies among B-spline bases 
            #
            # Support functions: glass.wam(), gam.slist(), gam.wlist(), plot.glass()
            #
            # Reference: 
            # Eilers, P.H.C. and Marx, B.D. (2002). Generalized Linear Additive Smooth Structures. 
            #	Journal of Compuational and Graphical Statistics, 11(4): 758-783.
            # Marx, B.D. and Eilers, P.H.C. (1998). Direct generalized linear modeling with 
            #	penalized likelihood. CSDA, 28(2): 193-209.	
            # Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines and penalties 
            #	(with comments and rejoinder). Statistical Science, 11(2): 89-121.
            #
            #
            # (c) 1998, 2003 Brian D. Marx
            #
            number.knots <- ps.intervals + 2 * degree + 1
            smooth.dummy <- F
            if(!missing(x.signal)) {
                  x.signal <- as.matrix(x.signal)
                  x.index <- signal.index
                  if(length(x) != nrow(x.signal)) {
                        warning(paste(
                              "\nThe x argument needs to have the same length as nrow(x.signal)"
                        ))
                        if(missing(signal.index)) {
                              signal.index <- 1:ncol(x.signal)
                        }
                  }
            }
            if(!missing(varying.index)) {
                  x.index <- varying.index
                  if(length(x) != length(x.index)) {
                        warning(paste(
                              "length of x and varying.index should be equal"
                        ))
                  }
            }
            else if(missing(x.signal) && missing(varying.index)) {
                  x.index <- x
                  smooth.dummy <- T
            }
            if(!missing(x.signal) && !missing(varying.index)) {
                  warning(paste("varying coef models should not be used with signal/curve\nregression"
                  ))
            }
            xl <- min(x.index)
            xr <- max(x.index)
            xmax <- xr + 0.01 * (xr - xl)
            xmin <- xl - 0.01 * (xr - xl)
            dx <- (xmax - xmin)/ps.intervals
            nx <- names(x.index)
            x.index <- as.vector(x.index)
            nax <- is.na(x.index)
            if(nas <- any(nax))
                  x.index <- x[!nax]
            sorder <- degree + 1
            if(!missing(ps.intervals)) {
                  nAknots <- ps.intervals - 1
                  if(nAknots < 1) {
                        nAknots <- 1
                        warning(paste("ps.intervals was too small; have used 2"
                        ))
                  }
                  if(nAknots > 0) {
                        Aknots <- seq(from = xmin - degree * dx, to = xmax + 
                                            degree * dx, by = dx)
                  }
                  else knots <- NULL
            }
            sbasis <- spline.des(Aknots, x.index, sorder, 0 * x.index)$design
            if(!missing(x.signal) && missing(varying.index)) {
                  signal.basis <- sbasis
                  basis <- x.signal %*% sbasis
            }
            if(!missing(varying.index) && missing(x.signal)) {
                  varying.basis <- sbasis
                  basis <- x * sbasis
            }
            if(missing(x.signal) && missing(varying.index)) {
                  basis <- sbasis
                  sbasis <- NULL
            }
            n.col <- ncol(basis)
            if(nas) {
                  nmat <- matrix(NA, length(nax), n.col)
                  nmat[!nax,  ] <- basis
                  basis <- nmat
            }
            dimnames(basis) <- list(1:nrow(basis), 1:n.col)
            if((order - n.col + 1) > 0) {
                  order <- n.col - 1
                  warning(paste("order was too large; have used ", n.col - 1))
            }
            if(lambda < 0) {
                  lambda <- 0
                  warning(paste("lambda was negative; have used ", lambda))
            }
            if(lambda > 10000) {
                  lambda <- 10000
                  warning(paste("lambda was >10000; for stability have used", 
                                lambda))
            }
            aug <- diag(n.col)
            if(order != 0) {
                  for(tt in 1:order) {
                        aug <- diff(aug)
                  }
            }
            pen.aug <- sqrt(lambda) * aug
            attr(basis, "knots") <- Aknots
            if(ridge.adj == 0) {
                  attr(basis, "pen.augment") <- pen.aug
            }
            if(ridge.adj != 0) {
                  attr(basis, "pen.augment") <- rbind(pen.aug, sqrt(ridge.adj) * 
                                                            diag(n.col))
            }
            attr(basis, "lambda") <- lambda
            attr(basis, "order") <- order
            signal.dummy <- varying.dummy <- F
            attr(basis, "signal.basis") <- attr(basis, "varying.basis") <- NULL
            if(!missing(x.signal)) {
                  attr(basis, "signal.basis") <- signal.basis
                  signal.dummy <- T
            }
            attr(basis, "signal.dummy") <- signal.dummy
            if(!missing(varying.index)) {
                  attr(basis, "varying.basis") <- varying.basis
                  varying.dummy <- T
            }
            attr(basis, "smooth.dummy") <- smooth.dummy
            attr(basis, "varying.dummy") <- varying.dummy	
            #	attr(basis, "ridge.gam") <- ridge.gam
            attr(basis, "ridge.inv") <- ridge.inv
            attr(basis, "ps.int") <- ps.intervals
            attr(basis, "degree") <- degree
            basis
      }

