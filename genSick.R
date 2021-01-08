#--- function: genSick
# OBJECTIVE: generate sickness trajectory from given info (can be used in paraboot also)
# INPUTS:     dat         : (data.frame) data frame with id, cluster, treatment, time, and outcome information. Of the form of genNet()[[2]]
#             sick0       : (integer) Number of individuals with the outcome at time point t
#             ep          : (scalar) Number between 0 and 1. Probability of getting sick from outside of the observed network
#             immunity    : (boolean) An indicator of whether or not immunity is incurred after experiencing the outcome. Default is FALSE
#             net         : (array [n, n, tau+1]) network information; 1 for contact between i,j, 0 otherwise. Of the form of genNet()[[3]]
#             trans.prob  : (vector) transmission probability for each individual in the network
# OUTPUTS:    data:       : data.frame built on input dat, with additional information on outcome at next time point t+1, y.prop, and y.bin

genSick <- function(dat, sick0, ep, immunity, net, trans.prob){
  N       <- length(unique(dat$id))
  t.pts   <- length(unique(dat$t))
  sick.t  <- matrix(-9, nrow = N, ncol = (t.pts)) 

  #=== filling in sickness information sequentially
  sick.t[ ,1] <- sick0
  
  if (immunity == TRUE){
    for (i in 2:t.pts){
      sick.already  <- rowSums(sick.t, na.rm = TRUE)
      sick.t[ ,i]   <- sickNext(trans.prob, 
                                net[ , ,(i-1)], 
                                sick.t[ ,(i-1)], 
                                ep, 
                                sickAlready = sick.already)
    }
  }
  if (immunity == FALSE){
    for (i in 2:t.pts){ 
      null.count    <- nrow(sick.t)
      sick.t[ ,i]   <- sickNext(trans.prob, 
                                net[ , ,(i-1)], 
                                sick.t[ ,(i-1)], 
                                ep, 
                                exclude     = FALSE, 
                                sickAlready = null.count)
    }

  }

  # putting buffer weeks for the boundaries
  sick.jnk  <- cbind(rep(-9, nrow(sick.t)), sick.t, rep(-9, nrow(sick.t))) 
  net.jnk   <- abind(net[ , ,1]*(-9), net, along = 3)
  times     <- sort(unique(dat$t))
  
  summ  <- do.call(rbind, lapply(2:(t.pts+1), function(x) { 
    summ.t  <- sickSummary(dat$id[dat$t == times[1]],  
                           sick.jnk[ ,x], 
                           sick.jnk[ ,(x+1)], 
                           dat$txt[dat$t == times[1]], 
                           net.jnk[ , ,x])
    res     <- cbind(summ.t, t = rep(times[x - 1], nrow(summ.t)) )
    
 
    return(res)
  })) 

  dat2 <- as.data.frame(summ)
  dat   <- merge(dat, dat2, all = T, sort = F)

  dat$s[dat$t == times[2]] <- sick.t[ ,t.pts]
  
  data <- with(dat, data.frame(id, cluster, txt, t, s, y.prop, y.bin))

  
  return(data)
}