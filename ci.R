#--- function:  ci
# OBJECTIVE: assess p-value for a test and given hypothesis for assessment of a CI
# INPUTS:     dat         : (data.frame) data frame with id, cluster, treatment, time, and outcome information. Of the form of genNet()[[2]]
#             net         : (array [n, n, tau+1]) network information; 1 for contact between i,j, 0 otherwise. Of the form of genNet()[[3]]
#             hyp         : (data.frame) Data frame with values for p0, p1, ep 
#             test.stat   : (function) Function that takes in data and returns a test statistic
#             B.sampDist  : (integer) Number of simualted datasets to use in the sampling distributions. Default is 100
#             alpha       : (scalar) Value between 0 and 1. Represents the test size for the desired $(1 - \alpha)$ CI. Default is 0.05
#             sided       : (string in ["one", "two"]) Indication of whether extremely large values ("one") or extremely large and small values ("two") are used for the test. Default is "two", but "one" is used for stat5
#             target      : (boolean) indicator of whether the results are to be used for a 3D CI for theta (FALSE) or a 1D interval for delta (TRUE) 
#             whichBd     : (string in ["lower", "upper"]) indicates lower or upper bound
#             deltaMLE    : (scalar) the delta_{\hat{theta}} value found via MLE for theta
# OUTPUTS:    pval        : if target == TRUE, pvalue of test with parameter combination hyp
#             objVal      : if target == TRUE, then the output is the hypothesized delta value times an indicator of whether 
#                           the param combo is in the confidence region 
#             results     : if target == FALSE, data frame with parameter combination, pvalue and indicator of whether the combo is in the CI
ci <- function(dat, net, hyp, test.stat, B.sampDist=100, alpha, sided="two", target, whichBd, deltaMLE){
  # number of param combinations we are assessing in this run
  until     <- ifelse(length(dim(hyp)) == 0, 1, nrow(hyp))               
  results   <- hyp

  if(target == "theta"){
    results[["pval"]]   <- t( sapply(1:until, function(i) para.boot(dat, net, hyp, test.stat, B.sampDist, sided) ) )
    results[["inC"]]    <- (results$pval > alpha)
    
    return(results)
  } 
  else if(target == "delta"){

    if(between(hyp$p1, 0, 1)){
      results[["pval"]]   <- para.boot(dat, net, hyp, test.stat, B.sampDist, sided)     
    } else{
      # coding for implausible combos resulting from the switch to searching on delta
      print("implausible p1 value considered")
      results[["pval"]]   <- -100                                                  
    }
    
    # sign for objective function (because it maximizes)
    beta                <- recode(whichBd, "lower" = -1, "upper" = 1)                
    gamma               <- 0.0001
    objVal              <- beta*((hyp$p1 - hyp$p0) - deltaMLE) / ((results$pval - alpha)^2 + gamma)
    
    return(c(pval = results$pval, V = objVal))
  }
}
