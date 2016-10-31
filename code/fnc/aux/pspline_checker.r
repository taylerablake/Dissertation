"pspline.checker"<-
      function(family, link, degree, order, ps.intervals, lambda, ridge.adj, wts)
      {
            if(link == "default" && family == "gaussian") {
                  link <- "identity"
            }
            if(link == "default" && family == "poisson") {
                  link <- "log"
            }
            if(link == "default" && family == "binomial") {
                  link <- "logit"
            }
            if(link == "default" && family == "Gamma") {
                  link <- "log"
            }
            if(family != "binomial" && family != "gaussian" && family != "poisson" && 
               family != "Gamma") {
                  warning(paste("Improper FAMILY option. Choose: gaussian, poisson, binomial or Gamma"
                  ))
            }
            if((family == "binomial") && (link != "logit" && link != "probit" && 
                                          link != "cloglog" && link != "loglog")) {
                  warning(paste("Improper LINK option with family=binomial. Choose: logit, probit, loglog, cloglog"
                  ))
            }
            if((family == "Gamma") && (link != "log" && link != "recipical" && link !=
                                       "identity")) {
                  warning(paste("Improper LINK option with family=Gamma. Choose: recipical, log, identity"
                  ))
            }
            if((family == "poisson") && (link != "log" && link != "sqrt" && link != 
                                         "identity")) {
                  warning(paste("Improper LINK option with family=poisson. Choose: log, sqrt, identity"
                  ))
            }
            if((family == "gaussian") && (link != "identity")) {
                  warning(paste("Improper LINK option with family=gaussian. Choose: identity"
                  ))
            }
            if(degree < 0) {
                  degree <- 1
                  warning(paste("degree must be non-neg integer: have used 1"))
            }
            if(order < 0) {
                  order <- 0
                  warning(paste("order must be non-neg integer: have used 0"))
            }
            if(ps.intervals < 2) {
                  ps.intervals <- 2
                  warning(paste("ps.intervals must be positive integer, > 1: have used 2"
                  ))
            }
            if(lambda < 0) {
                  lambda <- 0
                  warning(paste("lambda cannot be negative: have used 0"))
            }
            if(ridge.adj < 0) {
                  ridge.adj <- 0
                  warning(paste("ridge.adj cannot be negative: have used 0"))
            }
            if(min(wts) < 0) {
                  warning(paste("At least one weight entry is negative"))
            }
            llist <- list(family = family, link = link, degree = degree, order = 
                                order, ps.intervals = ps.intervals, lambda = lambda, ridge.adj
                          = ridge.adj, wts = wts)
            return(llist)
      }