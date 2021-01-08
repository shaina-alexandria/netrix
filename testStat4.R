#--- function: stat4
# OBJECTIVE: compute test statistic 4 (T_4 from manuscript) for a given dataset
# INPUTS:   data    : (data.frame) data frame with id, cluster, treatment, time, and outcome information. Of the form of genNet()[[2]]
#           network : (array [n, n, tau+1]) network information; 1 for contact between i,j, 0 otherwise. Of the form of genNet()[[3]]
# OUTPUT:   stat    : the value of the test statistic


stat4 <- function(data, network, ...){

  ace.hat <- function(y, txt, x = 1){
    mean(y[txt == 1 & x == 1], na.rm = T) - mean(y[txt == 0 & x == 1], na.rm = T)
  }
  
  data.y  <- subset(data, subset = (y.prop != -9))
  stat    <- with(data.y, ace.hat(y.prop, txt, s) )

  return(stat)
}
