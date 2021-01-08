#--- function: para.boot.comp1
# OBJECTIVE: reassign treatment, recalculate transmission probs, and regenerate data
# INPUTS:   data0       : (data.frame) data frame with id, cluster, treatment, time, and outcome information. Of the form of genNet()[[2]]
#           net         : (array [n, n, tau+1]) network information; 1 for contact between i,j, 0 otherwise. Of the form of genNet()[[3]]
#           out.i       : (data.frame) Data frame with values for p0, p1, ep 
#           test.stat   : (function) Function that takes in data and returns a test statistic
#           immunity    : (boolean) An indicator of whether or not immunity is incurred after experiencing the outcome. Default is FALSE
# OUTPUTS:  t           : value of test statistic for generated data

para.boot.comp1 <- function(dataO, net, out.i, test.stat, immunity=FALSE){

  #== rerandomize treatment
  clust           <- unique(dataO$cluster)                                            
  m               <- length(unique(dataO$cluster[dataO$txt == 1]))                     # number of unique treated clusters
  whoTxt          <- sample(clust, m)                                               # reassigning treatment
  txt.r           <- sapply(dataO$cluster, function(x) ifelse(any(x == whoTxt), 1, 0))
  data.n          <- dataO
  data.n[["txt"]] <- txt.r                                                          # copy of dataset but with a rerandomized treatment regime
  dat.1           <- data.n[ ,c("id", "t", "cluster", "txt")]
  times           <- sort(unique(dataO$t))

  #== using treatment information to generate new trans probs
  trans.prob0   <- out.i$p0 + (out.i$p1 - out.i$p0)*dat.1$txt[dat.1$t == times[1]]
  trans.prob    <- sapply(trans.prob0, function(x) max(x, 0))                        # making sure there are no negative trans.probs

  #== regenerate time series
  sick0   <- dataO$s[dataO$t == times[1]] 
  data    <- genSick(dat.1, sick0, out.i$ep, immunity, net, trans.prob)

  #== calc t from new data (for sampling dist)
  t   <- test.stat(data, net, out.i)[[1]]

  return(t)
}

#--- function: para.boot
# OBJECTIVE: calculate pvalue under alternative hyp described by out.i using para.boot.comp1
# INPUTS:   data0       : (data.frame) data frame with id, cluster, treatment, time, and outcome information. Of the form of genNet()[[2]]
#           net         : (array [n, n, tau+1]) network information; 1 for contact between i,j, 0 otherwise. Of the form of genNet()[[3]]
#           out.i       : (data.frame) Data frame with values for p0, p1, ep 
#           test.stat   : (function) Function that takes in data and returns a test statistic
#           B           : (integer) Number of simualted datasets to use in the sampling distributions. Default is 100
#           sided       : (string in ["one", "two"]) Indication of whether extremely large values ("one") or extremely large and small values ("two") are used for the test. Default is "two", but "one" is used for stat5
#           display     : (boolean) Indicator of whether to provide summaries, observed test stats, and plots for trobuleshooting. Default is FALSE
#           immunity    : (boolean) An indicator of whether or not immunity is incurred after experiencing the outcome. Default is FALSE
# OUTPUTS:  pval        : vector of pvalues for test 

para.boot <- function(dataO, net, out.i, test.stat, B = 100, sided = "two", display = FALSE, immunity = FALSE){
  boot.stats  <- replicate(B, para.boot.comp1(dataO, net, out.i, test.stat, immunity) )   # sampling distribution
  obsT        <- test.stat(dataO, net, out.i)[[1]]                              # observed test stat

  #== calc pval
  if(length(obsT) == 1){
    
    if(sided == "two"){
      boot.center <- mean(boot.stats)
      extreme     <- abs(obsT - boot.center) <= abs(boot.stats - boot.center)
    }
    if(sided == "one"){
      extreme     <- (obsT <= boot.stats)
    }
    
    pval  <- mean(extreme)
    
    if(display == TRUE){
      print(summary(boot.stats))
      print(paste("Tstat=",obsT))
    
      if(var(boot.stats) > 0){
        hist(boot.stats)
        abline(v = obsT, lwd = 2, col = "olivedrab")
      }
    }
  } 
  else{
    
    if(sided == "two"){
      boot.center <- rowMeans(boot.stats)
      extreme     <- sapply(1:length(obsT), 
                            function(i) (abs(obsT[i] - boot.center[i]) <= abs(boot.stats[i,] - boot.center[i])) 
      ) # dimensions flip
    }
    if(sided == "one"){
      extreme     <- sapply(1:length(obsT), 
                            function(i) (obsT[i] <= boot.stats[i,]) 
      ) # dimensions flip
    }

    pval  <- colMeans(extreme)
  }

  return(pval)
}
