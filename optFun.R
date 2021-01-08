#--- function: optFun
# OBJECTIVE: to calculate the objective function for determining the upper or lower bound of the 
#            confidence interval for delta with the internal print option for global access
# INPUTS:     x0          : (vector) initial values for theta, with arguments p0, p1, epsilon -- in that order 
#             deltaMLE    : (scalar) the delta_{\hat{theta}} value found via MLE for theta
#             whichBd     : (string in ["lower", "upper"]) indicates lower or upper bound
#             dat         : (data.frame) data frame with id, cluster, treatment, time, and outcome information. Of the form of genNet()[[2]]
#             net         : (array [n, n, tau+1]) network information; 1 for contact between i,j, 0 otherwise. Of the form of genNet()[[3]]
#             B.sampDist  : (integer) Number of simualted datasets to use in the sampling distributions. Default is 100
#             prt         : (boolean) option to print the itnernal results for use outside of the algorithm. boundSearch relies on this being TRUE. Default is TRUE
# OUTPUTS:    val         : value of the objective function form the ci function
#             delta       : the value of delta corresponding to the input theta
optFun <- function(x0, deltaMLE, whichBd, dat, net, B.sampDist = 100, prt = TRUE){
  hyp <- data.frame(p0 = x0[1], p1 = x0[2], ep = x0[3]) 
  val <- ci(dat,
            net,
            hyp,
            stat5,
            B.sampDist,
            alpha   = 0.05,
            sided   = "one",
            target  = "delta",
            whichBd,
            deltaMLE)
  
  delta   <- hyp$p1 - hyp$p0
  retVal  <- c(x0, pval = val[1], V = -val[2], delta = delta) # for finding the minimum with the optimizer
  
  if(prt == TRUE){
    optHistory <<- rbind(optHistory, retVal)
    
  }
  return(c(val, delta))
}