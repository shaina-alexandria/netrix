#--- function: perm.txt
# OBJECTIVE:
# INPUTS:     data      : (data.frame) data frame with id, cluster, treatment, time, and outcome information. Of the form of genNet()[[2]]
#             net       : (array [n, n, tau+1]) network information; 1 for contact between i,j, 0 otherwise. Of the form of genNet()[[3]]
#             test.stat : (function) Function that takes in data and returns a test statistic
# OUTPUTS:    ace.hat   : the value of the test statistic with the permuted dataset

perm.txt <- function(data, net., test.stat.){
  m         <- length(unique(data$cluster[data$txt == 1]))
  whoTxt    <- sample(unique(data$cluster), m)
  data$txt  <- (data$cluster %in% whoTxt)
  ace.hat   <- test.stat.(data, net.)

  return(ace.hat)
}

#--- function: perm.test
# OBJECTIVE:
# INPUTS:   data        : (list) List object of the form returned by genNet
#           test.stat   : (function) Function that takes in data and returns a test statistic
#           n.sims      : (integer) The number of permutations used to generate a sampling distribution for the test statistic. Default is 1000
# OUTPUTS:  pval        : the p-value from the RI test
#           stat        : the value of the test statistic for the observed data

perm.test <- function(data, test.stat, n.sims = 1000){
  data.b      <- data[[2]]
  net         <- data[[3]]

  ace.hat.obs <- test.stat(data.b, net)
  samp.dist   <- replicate(n.sims, perm.txt(data.b, net, test.stat) ) 

  dist.cent   <- samp.dist - mean(samp.dist, na.rm = T)
  obs.cent    <- ace.hat.obs - mean(samp.dist, na.rm = T)
  pval        <- sum(abs(dist.cent) >= abs(obs.cent))/n.sims 

  return(c(pval = pval, stat = obs.cent))
}






