#--- function: boundSearch
# OBJECTIVE: conduct the search function for the upper and lower bounds for the delta confidence interval
# INPUTS:     data        : (list) List object of the form returned by genNet
#             alpha       : (scalar) Value between 0 and 1. Represents the test size for the desired $(1 - \alpha)$ CI. Default is 0.05
#             B.sampDist  : (integer) Number of simualted datasets to use in the sampling distributions. Default is 100
#             B.optIters  : (integer) The maximum number of iterations nloptr can use to gather information on the location of the upper and lower bounds. Default is 100 
#             plot        : (boolean) Indicator of whether or not to return an illustrative plot with the MLE and bounds.  Default is FALSE
# OUTPUTS:    lb          : estimate of lower bound for delta CI
#             ub          : estimate of upper bound for delta CI
#             deltaMLE    : delta_{\hat{\theta}} based on MLE for theta

boundSearch <- function(data, alpha = 0.05, B.sampDist = 100, B.optIters = 100, plot = FALSE){
  dat     <- data[[2]]
  network <- data[[3]]
  
  #-- initial values for theta
  times   <- sort(unique(dat$t))
  start   <- optim(par    = sapply(c(0.1, 0.1, 0.01), logit),
                    fn     = likMW,
                    method = "BFGS",
                    dat    = dat,
                    net    = network)
  optStart  <- inv.logit(start$par)
  deltaMLE  <- optStart[2] - optStart[1]

  #--- using nloptr
  optHistory <<- data.frame(p0 = NA, p1 = NA, ep = NA, pval = NA, V = NA, delta = NA)     # initialize the results data frame


  evalLB <- nloptr(x0 = optStart,
                   eval_f = function(x) optFun(x, deltaMLE, "lower", dat, network, B.sampDist, prt = TRUE)["V"],
                   lb = rep(1e-7, 3), ub = rep( 1 - 1e-7, 3),
                   opts = list("algorithm" = "NLOPT_LN_NEWUOA_BOUND", "maxeval" = B.optIters, print_level = 0))

  evalUB <- nloptr(x0 = optStart,
                   eval_f = function(x) optFun(x, deltaMLE, "upper", dat, network, B.sampDist, prt = TRUE)["V"],
                   lb = rep(1e-7, 3), ub = rep( 1 - 1e-7, 3),
                   opts = list("algorithm" = "NLOPT_LN_NEWUOA_BOUND", "maxeval" = B.optIters, print_level = 0))
  
  
  histo  <- optHistory[-1, ]                                                            # removing the initialized values in this data frame
  
  #--- penalized spline
  modUp     <- scam(pval ~ s(delta, k = 5, bs = "mpd"), data = histo[histo$delta >= deltaMLE, ]) # monotonic decreasing
  modLo     <- scam(pval ~ s(delta, k = 5, bs = "mpi"), data = histo[histo$delta <= deltaMLE, ]) # monotonic increasing
  
  findUp    <- function(x) {(predict.scam(modUp, data.frame(delta = x)) - alpha)^2}
  findLo    <- function(x) {(predict.scam(modLo, data.frame(delta = x)) - alpha)^2}
  
  lb        <- optimize(findLo, interval = c(-1, deltaMLE))$minimum
  ub        <- optimize(findUp, interval = c(deltaMLE, 1))$minimum
  
  if(plot == TRUE){
    fnUp      <- function(x) predict.scam(modUp, data.frame(delta = x))
    fnLo      <- function(x) predict.scam(modLo, data.frame(delta = x))
  
    illPlot <- ggplot(histo, aes(delta, pval)) + 
      geom_point() + 
      theme_classic() +
      labs(x = expression(delta[theta^"*"]), y = expression(rho)) +
      geom_hline(yintercept = alpha, color = "gray", linetype = "dashed") +
      stat_function(fun = fnUp, xlim = c(deltaMLE, 1), colour = "red", size = 1) + 
      stat_function(fun = fnLo, xlim = c(-1, deltaMLE), colour = "blue", size = 1) +
      scale_x_continuous(limits = c(-1, 1), 
                         breaks = c(-1, lb, deltaMLE, ub, 1), 
                         labels = c("-1"     = "-1", 
                                    lb       = round(lb, 3), 
                                    deltaMLE = round(deltaMLE, 3), 
                                    ub       = round(ub, 3),
                                    "1"      = "1")) +
      scale_y_continuous(limits = c(0, 1),
                         breaks = c(alpha, 1)) +
      theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 13))
    
    print(illPlot)
  }
  

  #-- return ub, lb, and flag
  return(data.frame(lb, ub, deltaMLE))
  
  
}