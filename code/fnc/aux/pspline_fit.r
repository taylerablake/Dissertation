"pspline.fit"<-
      function(response, x.var, ps.intervals = 8, wts = NULL, degree = 3, order = 3, 
               link = "default", family = "gaussian", m.binomial = NULL, r.gamma = 
                     NULL, lambda = 0, x.predicted = NULL, ridge.adj = 0.0001)
      {
            # Function pspline.fit: univariate smoother using P-splines.
            # Input: x.var= explanatory variable on abcissae.
            # Input: response= response variable.
            # Input: family=gaussian, binomial, poisson, Gamma distribution.
            # Input: wts= vector of weights; default is vector of ones.
            # Input: m.binomial=vector of binomial trials. Default is 1 vector.
            # Input: r.gamma=vector of gamma shape parameters. Default is 1 vector.
            # Input: link= link function (identity, log, sqrt, logit, probit, cloglog, loglog, recipical).
            # Input: ps.intervals= number of intervals for B-splines. Default=8.
            # Input: degree= degree of B-splines. Default=3.
            # Input: order= order of difference penalty. Default=3.
            # Input: lambda= smoothness regulalizing parameter ( >= 0). Default=0.
            # Input: x.predicted=a list of x variables for prediction and twice stderr limits.
            # Result: a scatterplot of (response, x.var) with smoothed fit and se bands.
            # Output: A list: including AIC= deviance + 2*trace(Hat), dispers.parm, etc.
            #
            # Reference: Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines and
            #            penalties (with comments and rejoinder). Statistical Science, 11(2): 89-121.
            #
            #
            # Support functions: pspline.checker(), pspline.fitter(), pspline.predictor()
            #
            #
            # (c) 1995 Paul Eilers & Brian Marx
            #
            y <- response
            x <- x.var
            if(missing(wts)) {
                  wts <- rep(1, length(y))
            }
            parms <- pspline.checker(family, link, degree, order, ps.intervals, 
                                     lambda, ridge.adj, wts)
            family <- parms$family
            link <- parms$link
            q <- parms$degree
            d <- parms$order
            ridge.adj <- parms$ridge.adj
            lambda <- parms$lambda
            ndx <- parms$ps.intervals
            wts <- parms$wts
            if(missing(m.binomial)) {
                  m.binomial <- rep(1, length(y))
            }
            if(missing(r.gamma)) {
                  r.gamma <- rep(1, length(y))
            }
            n <- length(y)
            xl <- min(x)
            xr <- max(x)
            xmax <- xr + 0.01 * (xr - xl)
            xmin <- xl - 0.01 * (xr - xl)
            dx <- (xmax - xmin)/ndx
            knots <- seq(xmin - q * dx, xmax + q * dx, by = dx)
            b <- spline.des(knots, x, q + 1, 0 * x)$design
            n.col <- ncol(b)
            if(d < 0) {
                  d <- min(3, (n.col - 1))
                  warning(paste("penalty order cannot be negative: have used", d)
                  )
            }
            if((d - n.col + 1) > 0) {
                  d <- n.col - 1
                  warning(paste("penalty order was too large: have used", d))
            }
            if(ridge.adj > 0) {
                  nix.ridge <- rep(0, n.col)
                  p.ridge <- sqrt(ridge.adj) * diag(rep(1, n.col))
            }
            p <- diag(n.col)
            if(d != 0) {
                  for(j in 1:d) {
                        p <- diff(p)
                  }
            }
            p <- sqrt(lambda) * p
            nix <- rep(0, n.col - d)
            b <- as.matrix(b)
            ps.fit <- pspline.fitter(family, link, n.col, m.binomial, r.gamma, y, b,
                                     p, p.ridge, nix, nix.ridge, ridge.adj, wts)
            mu <- ps.fit$mu
            coef <- ps.fit$coef
            w <- ps.fit$w
            e <- 1e-009
            h <- hat(ps.fit$f$qr, intercept = F)[1:n]
            trace <- sum(h) - 1
            if(family == "binomial") {
                  dev <- 2 * sum((y + e) * log((y + e)/(mu + e)) + (m.binomial - 
                                                                          y + e) * log((m.binomial - y + e)/(m.binomial - mu + e)
                                                                          ))
                  dispersion.parm <- 1
            }
            if(family == "poisson") {
                  dev <- 2 * sum(y * log(y + e) - y - y * log(mu) + mu)
                  dispersion.parm <- 1
            }
            if(family == "Gamma") {
                  dev <- -2 * sum(r.gamma * (log((y + e)/mu) - ((y - mu)/mu)))
                  ave.dev <- dev/n
                  dispersion.parm <- (ave.dev * (6 + ave.dev))/(6 + 2 * ave.dev)
            }
            if(family == "gaussian") {
                  dev <- sum(ps.fit$f$residuals^2)
                  dispersion.parm <- dev/(n - trace)
            }
            aic <- dev + 2 * trace
            x.seq <- seq(xl, xr, length = 50)
            b.seq <- spline.des(knots, x.seq, q + 1, 0 * x.seq)$design
            w.aug <- c(w, (nix + 1))
            yhat <- b.seq %*% as.vector(ps.fit$coef)
            half.meat <- sqrt(c(w)) * b
            meat <- t(half.meat) %*% half.meat
            if(ridge.adj > 0) {
                  bread <- solve(meat + t(p) %*% p + t(p.ridge) %*% p.ridge)
            }
            if(ridge.adj == 0) {
                  bread <- solve(meat + t(p) %*% p)
            }
            half.sw <- half.meat %*% bread
            var.beta <- t(half.sw) %*% half.sw
            var.yhat <- b.seq %*% var.beta %*% t(b.seq)
            stdev.yhat <- as.vector(sqrt(diag(var.yhat)))
            stdev.yhat <- sqrt(dispersion.parm) * stdev.yhat
            pivot <- 2 * stdev.yhat
            upper <- yhat + pivot
            lower <- yhat - pivot
            summary.yhat <- cbind(lower, yhat, upper)
            if(link == "logit") {
                  summary.yhat <- 1/(1 + exp( - summary.yhat))
            }
            if(link == "probit") {
                  summary.yhat <- apply(summary.yhat, c(1, 2), pnorm)
            }
            if(link == "cloglog") {
                  summary.yhat <- (1 - exp( - exp(summary.yhat)))
            }
            if(link == "loglog") {
                  summary.yhat <- exp( - exp( - summary.yhat))
            }
            if(link == "sqrt") {
                  summary.yhat <- summary.yhat^2
            }
            if(link == "log") {
                  summary.yhat <- exp(summary.yhat)
            }
            if(link == "recipical") {
                  summary.yhat <- 1/(summary.yhat)
            }
            if(family == "binomial" && mean(m.binomial) != 1) {
                  matplot(x.seq, summary.yhat, type = "l", lty = c(2, 1, 2), xlab
                          = "regressor", ylab = "estimated mean", main = 
                                "P-spline fit with twice std error bands")
            }
            if(mean(m.binomial) == 1) {
                  matplot(x.seq, summary.yhat, type = "l", lty = c(2, 1, 2), xlab
                          = "regressor", ylab = "estimated mean", ylim = c(min(
                                min(y), min(summary.yhat[, 1])), max(max(y), max(
                                      summary.yhat[, 3]))), main = 
                                "P-spline fit with twice std error bands")
                  matpoints(x, y, type = "p", pch = "O")
            }
            ps.predict <- NULL
            if(!missing(x.predicted)) {
                  ps.predict <- pspline.predictor(x.predicted, knots, link, coef, 
                                                  q, var.beta, dispersion.parm)
            }
            llist <- list()
            llist$family <- family
            llist$link <- link
            llist$ps.intervals <- ndx
            llist$knots <- knots
            llist$order <- d
            llist$degree <- q
            llist$lambda <- lambda
            llist$aic <- aic
            llist$fits <- spline.des(knots, x, q + 1, 0 * x.seq)$design%*%as.vector(ps.fit$coef)
            llist$deviance <- dev
            llist$eff.df <- trace
            llist$df.resid <- n - trace
            llist$dispersion.param <- dispersion.parm
            llist$summary.predicted <- ps.predict$summary.pred
            llist$coef <- coef
            llist
      }



