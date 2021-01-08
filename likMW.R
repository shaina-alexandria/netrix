#--- function: likMW
# OBJECTIVE: to calculate the loglikelihood for multiple weeks given observed data and parameter guesses and return the negative loglikelihood
# INPUTS:   par   : (vector) with arguments p0, p1, epsilon -- in that order 
#           dat   : (data.frame) data frame with id, cluster, treatment, time, and outcome information. Of the form of genNet()[[2]]
#           net   : (array [n, n, tau+1]) network information; 1 for contact between i,j, 0 otherwise. Of the form of genNet()[[3]]
# OUTPUTS:  neglik: negative log likelihood 
likMW <- function(par, dat, net){

  # a list of all the time points to be utilized
  times <- sort(unique(dat$t))
  
  # an initial choice for the parameters
  par1  <- sapply(par, function(x) (inv.logit(x) + 1e-7)*(1 - 2e-7))
  
  # initializing the vector of likelihood contributions
  l.t   <- rep(NA, length(times) - 1)
  
  # calculate the likelihood contribution from each time point alone and sum at the end
  for(t in times[2:length(times)]){
    
    # defining the time invariate txt A
    A    <- dat$txt[dat$t == (t-1)]
    
    # defining S0 for a specific time pair
    S0   <- dat$s[dat$t == (t-1)]
    
    # defining S1 for a specific time pair
    S1   <- dat$s[dat$t == t]
    
    # defining E for a specific time pair (since arrays don't have zero indices, t represent network data for time t-1)
    E    <- net[ , ,t]
    
    # calculating the inner part of the log in the likelihood
    q    <- ( A*par1[2] + (1 - A)*par1[1] ) * S0
    
    # turning the vector into a matrix so I know how it's multiplying
    q1   <- matrix(q, nrow = length(q), ncol = length(q))   
    
    # the resultant total risk from contacts
    Q    <- colProds(1 - (E * q1))
    
    # the likelihood conribution for time t
    l.t[t]  <- sum(S1 * log(1 - ((1 - par1[3]) * Q))) + sum((1 - S1) * log((1 - par1[3]) * Q))
  }
  
  # summing to get the full likelihood contribution
  l      <- sum(l.t)
  
  # returning the negative log likelihood or a lower bound
  negLik <- max(-1e8, -1*l)
  
  if(is.finite(negLik) == FALSE){browser()}
  return(negLik)
}