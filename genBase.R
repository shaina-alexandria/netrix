#--- function: genBase
# OBJECTIVE: generate baseline sickness and calculate trans.prob
# INPUTS:     p0          : (scalar) Number between 0 and 1. Per-contact probability of transmission for an individual in the control group 
#             p1          : (scalar) Number between 0 and 1. Per-contact probability of transmission for an individual in the treatment group
#             ep          : (scalar) Number between 0 and 1. Probability of getting sick from outside of the observed network
#             N           : (integer) Number of individuals in the network
#             h           : (integer) Number of clusters in the network 
#             m           : (integer) Number of clusters to be assignmed to the treatment group (m < h)
#             s           : (integer) Number of individuals with the outcome at baseline (time 0)
#             m.friends   : (scalar) Mean number of friends individuals have within the network (used for ERGM model selection only) 
#             t.pts       : (integer) Number of time points, including baseline (tau + 1)
#             immunity    : (boolean) An indicator of whether or not immunity is incurred after experiencing the outcome. Default is FALSE
# OUTPUTS:    dat.1       : dataframe with id, cluster, txt, and t information for all t.pts time points. To be built on by genNet
#             sick0       : indicator of outcome at baseline 
#             trans.prob  : transmission probability for each individual in the network
genBase <- function(p0, p1, ep, N, h, m, s, m.friends, t.pts, immunity){
  id      <- 1:N
  hall.0  <- sample(rep(1:h, ceiling(N/h)))
  hall    <- hall.0[1:N]
  whoTxt  <- sample(1:h, m)
  txt     <- sapply(hall, function(x) ifelse(any(x == whoTxt), 1, 0))
  whoSick <- sample(1:N, s)
  sick0   <- sapply(id, function(x) ifelse(any(x == whoSick), 1, 0))
  
  trans.prob0 <- p0*(1 - txt) + p1*txt
  # making sure there are no negative trans.probs
  trans.prob  <- sapply(trans.prob0, function(x) max(x,0))          
  
  dat.1   <- data.frame(id      = rep(id, t.pts), 
                        cluster = rep(hall, t.pts), 
                        txt     = rep(txt, t.pts), 
                        t       = rep(0:(t.pts-1), each = N)
                        )                                          
  return(list(dat.1, sick0, trans.prob))
}
