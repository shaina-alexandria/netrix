# OBJECTIVE: likelihood function used for test statistic calculation
# INPUTS:  par : parameter set (null) against which to compare the observed data
#          A   : vector of treatment assignments in the observed data
#          S0  : vector of illness status at time 0
#          S1  : vector of illness status at time 1
#          E   : observed network (matrix) of contact information
# OUTPUTS: negLik: negative likelihood of binomial distribution given the observed data

lik <- function(par, A, S0, S1, E){
  
  par1 <- sapply(par, inv.logit)
  q    <- ( A*par1[2] + (1-A)*par1[1] ) * S0
  q1   <- matrix(q, nrow = length(q), ncol = length(q))          # turning the vector into a matrix so I know how it's multiplying
  Q    <- colProds(1 - (E * q1))
  l    <- sum(S1 * log(1 - ((1 - par1[3]) * Q))) + sum((1 - S1) * log((1 - par1[3]) * Q))
  
  negLik <- -1*l
  return(negLik)
}