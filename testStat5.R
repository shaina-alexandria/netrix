#--- function: stat5
# OBJECTIVE: compute test statistic 5 (T_5 from manuscript) for a given dataset
# INPUTS:   data    : (data.frame) data frame with id, cluster, treatment, time, and outcome information. Of the form of genNet()[[2]]
#           network : (array [n, n, tau+1]) network information; 1 for contact between i,j, 0 otherwise. Of the form of genNet()[[3]]
#           hyp     : (data.frame) Data frame with values for p0, p1, ep 
#           start   : (vector) initial values for theta for MLE search
# OUTPUT:   stat    : the value of the test statistic

stat5 <- function(data, network, hyp, start = c(0.1, 0.1, 0.01), ...){

  # a list of all the time points to be utilized
  times         <- sort(unique(data$t))
  
  # checking whether parameters will be able to be identified. If not, skip nloptr
  identifiable  <- (var(data$s[data$t == times[2]], na.rm = TRUE) > 0)               
  
  # throwing an error if the start values are not finite
  if(all(is.finite(sapply(start, logit))) == FALSE){browser()}
  
  # defining the test stat to be zero when unidentifiable
  if(identifiable == 0){
    stat <- 0                                                     
  }
  else{
    
    # using optim to find the MLEs
    ests <- optim(par    = sapply(start, logit),
                  fn     = likMW,
                  method = "BFGS",
                  dat    = data,
                  net    = network)
    
    # reconverting estimates from the logit scale
    stat0           <- as.data.frame(t(sapply(ests$par, inv.logit)))
    
    # naming the columns
    colnames(stat0) <- c("p0","p1","ep")
    
    # calculating the test stat based on the hypothesis AND MLE estimates from observed data
    stat  <- (stat0$p0 - hyp$p0)^2 + (stat0$p1 - hyp$p1)^2 + (stat0$ep - hyp$ep)^2
  }
  
  return(stat)
}
