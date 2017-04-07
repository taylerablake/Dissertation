pspline2d.checker <- function(family, link, degree1, degree2, order1, order2, ps.intervals1, 
               ps.intervals2, lambda1, lambda2, ridge.adj, wts){
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
            if(degree1 < 0) {
                  degree1 <- 1
                  warning(paste("degree1 must be non-neg integer: have used 1"))
            }
            if(order1 < 0) {
                  order1 <- 0
                  warning(paste("order1 must be non-neg integer: have used 0"))
            }
            if(ps.intervals1 < 2) {
                  ps.intervals1 <- 2
                  warning(paste("ps.intervals1 must be positive integer, > 1: have used 2"
                  ))
            }
            if(lambda1 < 0) {
                  lambda1 <- 0
                  warning(paste("lambda1 cannot be negative: have used 0"))
            }
            if(degree2 < 0) {
                  degree2 <- 1
                  warning(paste("degree2 must be non-neg integer: have used 1"))
            }
            if(order2 < 0) {
                  order2 <- 0
                  warning(paste("order2 must be non-neg integer: have used 0"))
            }
            if(ps.intervals2 < 2) {
                  ps.intervals2 <- 2
                  warning(paste("ps.intervals2 must be positive integer, > 1: have used 2"
                  ))
            }
            if(lambda2 < 0) {
                  lambda2 <- 0
                  warning(paste("lambda2 cannot be negative: have used 0"))
            }
            if(ridge.adj < 0) {
                  ridge.adj <- 0
                  warning(paste("ridge.adj cannot be negative: have used 0"))
            }
            if(min(wts) < 0) {
                  warning(paste("At least one weight entry is negative"))
            }
            llist <- list(family = family, link = link, degree1 = degree1, order1
                          = order1, ps.intervals1 = ps.intervals1, lambda1 = lambda1, 
                          degree2 = degree2, order2 = order2, ps.intervals2 = 
                                ps.intervals2, lambda2 = lambda2, ridge.adj = ridge.adj, wts = 
                                wts)
            return(llist)
      }