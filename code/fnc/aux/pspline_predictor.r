"pspline.predictor"<-
      function(x.predicted, knots, link, coef, q, var.beta, dispersion.parm, ...)
      {
            b.pred <- spline.des(knots, x.predicted, q + 1, 0 * x.predicted, outer.ok = TRUE)$design
            eta.pred <- b.pred %*% as.vector(coef)
            b.pred <- as.matrix(b.pred)
            if(length(x.predicted) > 1) {
                  var.pred <- (b.pred) %*% var.beta %*% t(b.pred)
            }
            if(length(x.predicted) == 1) {
                  var.pred <- t(b.pred) %*% var.beta %*% (b.pred)
            }
            stdev.pred <- as.vector(sqrt(diag(var.pred)))
            stdev.pred <- sqrt(dispersion.parm) * stdev.pred
            pivot <- as.vector(2 * stdev.pred)
            upper <- eta.pred + pivot
            lower <- eta.pred - pivot
            summary.pred <- cbind(lower, eta.pred, upper)
            if(link == "logit") {
                  summary.pred <- 1/(1 + exp( - summary.pred))
            }
            if(link == "probit") {
                  summary.pred <- apply(summary.pred, c(1, 2), pnorm)
            }
            if(link == "cloglog") {
                  summary.pred <- (1 - exp( - exp(summary.pred)))
            }
            if(link == "loglog") {
                  summary.pred <- exp( - exp( - summary.pred))
            }
            if(link == "sqrt") {
                  summary.pred <- summary.pred^2
            }
            if(link == "log") {
                  summary.pred <- exp(summary.pred)
            }
            if(link == "recipical") {
                  summary.pred <- summary.predd <- 1/(summary.pred)
                  summary.pred <- summary.predd[, 3:1]
            }
            summary.pred <- as.matrix(summary.pred)
            dimnames(summary.pred) <- list(NULL, c("-2std_Lower", "Predicted", 
                                                   "+2std_Upper"))
            llist <- list(summary.pred = summary.pred)
            return(llist)
      }