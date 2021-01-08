#--- function: genNet
# OBJECTIVE: to generate a data set with id, hall number, treatment, number of friends, sick at time 0, binary outcome y, 
#            and proportional outcome y to pass through to the permutation test using a scale-free model
# INPUTS:     p0          : (scalar) Number between 0 and 1. Per-contact probability of transmission for an individual in the control group 
#             p1          : (scalar) Number between 0 and 1. Per-contact probability of transmission for an individual in the treatment group
#             ep          : (scalar) Number between 0 and 1. Probability of getting sick from outside of the observed network
#             netgenmodel : (string in ["scale-free", "ergm"]) model used to generate network. Default is "scale-free"
#             N           : (integer) Number of individuals in the network
#             h           : (integer) Number of clusters in the network 
#             m           : (integer) Number of clusters to be assignmed to the treatment group (m < h)
#             s           : (integer) Number of individuals with the outcome at baseline (time 0)
#             m.friends   : (scalar) Mean number of friends individuals have within the network (used for ERGM model selection only) 
#             tau         : (integer) Number of time points, subsequent to baseline. Default is 9 (like in eX-FLU)
#             immunity    : (boolean) An indicator of whether or not immunity is incurred after experiencing the outcome. Default is FALSE
# OUTPUTS:   list with 3 elements:
#            1. netgenmodel : input string returned
#            2. data        : dataframe with information on id, cluster, txt, s (indicator of outcome), 
#                             y.prop (the proportion of friends who got sick in the week after the ego was sick),
#                             y.bin (indicator of if any friends got sick in the week after the ego was sick)
#            3. net         : network for tau+1 time points and N individuals

genNet <- function(p0, p1, ep=0.01, netgenmodel = "scale-free", N=504, h=112, m=56, s=25, m.friends=4, tau=9, immunity=FALSE){
  # t.points is the total number of time points, including baseline
  t.pts     <- tau + 1
   
  #=== generating baseline information
  dat.b     <- genBase(p0, p1, ep, N, h, m, s, m.friends, t.pts, immunity)

  # initialize network for all time points
  net       <- array(0, dim = c(N, N, t.pts))
  
  if(netgenmodel == "scale-free"){
    #=== creating the network (based on scale-free network model)
    # initial set of nodes
    m0                <- 2
    
    # empty initial network
    net0              <- matrix(0, nrow = N, ncol = N)
    
    # first step: all m0 existing nodes are connected to node m0+1
    net0[1:m0, m0+1]  <- 1
    net0[m0+1, 1:m0]  <- 1
    
    # dynamically add edges with each new node from m0+2 to N
    for(k in (m0 + 2):N){
      # calculate degree of each node so far
      degree    <- rowSums(net0)[1:(k-1)]
      
      # probability of connecting to node i based on degree
      nodeProb  <- degree / sum(degree)
      
      # make m0 connections for node k based on nodeProb 
      whichEdge <- sample(1:(k-1), m0, prob = nodeProb)
      
      # add to network
      net0[whichEdge, k]  <- 1
      net0[k, whichEdge]  <- 1
      
    }
    
    # populating baseline network
    net[ , ,1]  <- net0
  }
  
  if(netgenmodel == "ergm"){
    #=== creating the network -- ergm
    network.0   <- rbinom(N*(N-1)/2, 1, m.friends/N)
    network.1   <- matrix(0, nrow = N, ncol = N)
    network.1[lower.tri(network.1)] <- network.0
    network     <- t(network.1)
    network[lower.tri(network)]     <- network.0   
    net[ , ,1]  <- network
    
  }
  
  if(!(netgenmodel %in% c("scale-free", "ergm"))){message("model must be 'scale-free' or 'ergm'"); break}
  
  if (t.pts > 1){
    
    for (i in 2:t.pts){
      # network at previous time point
      ni          <- net[ , ,(i-1)]
      
      # total number of connections
      totF        <- sum(ni)
      
      # sampling 1% of contacts to be dropped
      drops       <- sample(which(ni[lower.tri(ni)] == 1), round(0.1*totF))
      
      # sampling 1% of non-contacts to be added
      adds        <- sample(which(ni[lower.tri(ni)] == 0), round(0.1*totF))
      
      # unlisted matrix values
      newF        <- ni[lower.tri(ni)]
      
      # alterations
      newF[drops] <- 0
      newF[adds]  <- 1
      
      # reformatting as symmetric matrix
      ni[lower.tri(ni)] <- newF
      ni                <- t(ni)
      ni[lower.tri(ni)] <- newF
      net[ , ,i]        <- ni  
      
      if( isSymmetric(net[ , ,i]) == FALSE ) warning("network matrix not symmetric")
    }
  }

  #=== generating sickness at time t
  data  <- genSick(dat.b[[1]], dat.b[[2]], ep, immunity, net, dat.b[[3]])

  return(list(netgenmodel, data, net))
}
