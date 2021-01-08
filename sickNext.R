#--- function: sickNext 
# OBJECTIVE: calculate probabilities of sickness at time 1 for the network and generate sickness indicator 
# INPUTS:     trans.prob  : (vector) transmission probability for each individual in the network
#             net         : (matrix [n, n]) network information for time point t; 1 for contact between i,j; 0 otherwise.
#             sick0       : (vector) Indicator of individuals with the outcome at time point t
#             ep          : (scalar) Number between 0 and 1. Probability of getting sick from outside of the observed network
#             exclude     : (boolean) An indicator of whether or not immunity is incurred after experiencing the outcome. Default is TRUE
#             sickAlready : (vector) indicator of whether a person has already been sick in the past (works with exclude)
# OUTPUTS:    sick1:      : vector of outcomes at next time point 

sickNext <- function(trans.prob, net, sick0, ep, exclude = TRUE, sickAlready){
  inner <- 1 - (net * (sick0 * trans.prob)) # by row here, so next by col
  r     <- 1 - ((1 - ep) * colProds(inner))
  sick  <- rbinom(length(r), size = 1, prob = r)
  
  if(exclude == TRUE){
    sick1 <- ifelse(sickAlready == 1, 0, sick)              # letting people be immune
  } else{
    sick1 <- sick                                           # not letting people be immune
  }

  return(sick1)
}
