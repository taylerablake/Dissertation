"pspline.fitter"<-
      function(family, link, n.col, m.binomial, r.gamma, y, b, p, p.ridge, nix, 
               nix.ridge, ridge.adj, wts, ...)
      {
            coef.est <- rep(1, ncol(b))
            if(family == "binomial") {
                  mu <- (y + 0.5 * m.binomial)/2
            }
            if(family == "Gamma" || family == "poisson") {
                  mu <- (y + 3)
            }
            if(family == "gaussian") {
                  mu <- rep(mean(y), length(y))
            }
            it <- 0
            repeat {
                  if(it == 0) {
                        if(link == "identity") {
                              eta <- mu
                        }
                        if(link == "log") {
                              eta <- log(mu)
                        }
                        if(link == "sqrt") {
                              eta <- sqrt(mu)
                        }
                        if(link == "logit") {
                              eta <- log(mu/(m.binomial - mu))
                        }
                        if(link == "recipical") {
                              eta <- 1/mu
                        }
                        if(link == "probit") {
                              eta <- qnorm(mu/m.binomial)
                        }
                        if(link == "cloglog") {
                              eta <- log( - log(1 - mu/m.binomial))
                        }
                        if(link == "loglog") {
                              eta <-  - log( - log(mu/m.binomial))
                        }
                  }
                  it <- it + 1
                  if(it > 25)
                        break
                  if(link == "identity") {
                        mu <- eta
                        h.prime <- 1
                  }
                  if(link == "log") {
                        mu <- exp(eta)
                        h.prime <- mu
                  }
                  if(link == "sqrt") {
                        mu <- eta^2
                        h.prime <- 2 * eta
                  }
                  if(link == "logit") {
                        mu <- m.binomial/(1 + exp( - eta))
                        h.prime <- mu * (1 - mu/m.binomial)
                  }
                  if(link == "recipical") {
                        mu <- 1/eta
                        h.prime <-  - (mu^2)
                  }
                  if(link == "probit") {
                        mu <- m.binomial * pnorm(eta)
                        h.prime <- m.binomial * dnorm(eta)
                  }
                  if(link == "cloglog") {
                        mu <- m.binomial * (1 - exp( - exp(eta)))
                        h.prime <- (m.binomial) * exp(eta) * exp( - exp(eta))
                  }
                  if(link == "loglog") {
                        mu <- m.binomial * exp( - exp( - eta))
                        h.prime <- m.binomial * exp( - eta) * exp( - exp( - eta
                        ))
                  }
                  if(family == "gaussian") {
                        w <- rep(1, length(y))
                  }
                  if(family == "poisson") {
                        w <- h.prime^2/mu
                  }
                  if(family == "binomial") {
                        w <- h.prime^2/(mu * (1 - mu/m.binomial))
                  }
                  if(family == "Gamma") {
                        w <- (r.gamma * h.prime^2)/mu^2
                  }
                  u <- (y - mu)/h.prime + eta
                  if(ridge.adj > 0) {
                        f <- lsfit(rbind(b, p, p.ridge), c(u, nix, nix.ridge), 
                                   wt = c(wts, nix + 1, nix.ridge + 1) * c(w, (nix +
                                                                                     1), (nix.ridge + 1)), intercept = F)
                  }
                  if(ridge.adj == 0) {
                        f <- lsfit(rbind(b, p), c(u, nix), wt = c(wts, nix + 1) *
                                         c(w, (nix + 1)), intercept = F)
                  }
                  coef.old <- coef.est
                  coef.est <- as.vector(f$coef)
                  d.coef <- max(abs((coef.est - coef.old)/coef.old))
                  if(d.coef < 1e-008)
                        break
                  print(c(it, d.coef))
                  eta <- b %*% coef.est
            }
            if(it > 24) {
                  warning(paste("parameter estimates did NOT converge in 25 iterations"
                  ))
            }
            llist <- list(coef = coef.est, mu = mu, f = f, w = w * wts)
            return(llist)
      }