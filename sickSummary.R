#--- function: sickSummary
# OBJECTIVE: calculates summary information of sickness required by the test statistics
# INPUTS:   ids       : (vector) participant ids
#           sick1     : (vector) Indicator of individuals with the outcome at time point t. Will be a vector of -9s when t+1 = 0
#           sick2     : (vector) Indicator of individuals with the outcome at time point t+1. Will be a vector of -9s when t = tau
#           txt       : (vector) Indicator of whether an individual is in the treatment group
#           network1  : (matrix [n, n]) network information for time point t; 1 for contact between i,j; 0 otherwise.
# OUTPUTS:  res       : data frame with information on id, y.bin, y.prop, and s (same as input sick1)

sickSummary <- function(ids, sick1, sick2, txt, network1){
    # calculating information to be stored in time "1"
  stopifnot(is.vector(sick1) == TRUE)
  
  # y.prop: uses network from time t  and illness from time t+1
  net1        <- as.matrix(network1)
  # (i,j) is 1 if i and j are friends at time 1 and person j is sick at t+1, 
  #          0 if i and j aren't friends OR i and j are friends and j is not sick, 
  #          -9 if information is missing
  
  # total number of friends
  n.f1      <- rowSums(net1) 
  
  # total number of friends from week 1 sick at time 2
  n.f.sick1 <- colSums(sick2*net1)             
  
  propSick1 <- ifelse(n.f1 > 0, n.f.sick1/n.f1, -9)     
  anySick1  <- as.numeric(ifelse(propSick1 == -9, -9, ifelse(propSick1 == 0, 0, 1)))
  
  res       <- data.frame(id     = ids, 
                          y.bin  = anySick1, 
                          y.prop = propSick1, 
                          s      = sick1)
  
  return(res)
}