#--- function: stat1
# OBJECTIVE: compute test statistic 1 (T_1 from manuscript) for a given dataset
# INPUTS:   data    : (data.frame) data frame with id, cluster, treatment, time, and outcome information. Of the form of genNet()[[2]]
#           network : (array [n, n, tau+1]) network information; 1 for contact between i,j, 0 otherwise. Of the form of genNet()[[3]]
# OUTPUT:   stat    : the value of the test statistic

stat1 <- function(data, network, ...){
  
  #----- recalculating the x.u based on new treatment assignments
  times     <- sort(unique(data$t))
  
  for(t in times[-1]){
    net0      <- as.matrix(network[ , ,t])
    txt       <- data$txt[data$t == (t-1) ]
    n.f0      <- rowSums(pmax(net0, 0))                                    # total number of friends
    x.u       <- ifelse(n.f0 != 0, colSums(txt*pmax(net0, 0))/n.f0 , -9)   # proportion of friends treated
  
    data$x.u[data$t == t] <- x.u
  }
  
  #----- test statistic
  data.y  <- subset(data, subset = (x.u != -9 & n.f0 > 0 & t != times[1]))
  mod3    <- glm(s ~ x.u, data = data.y, family = binomial)
  #mod3    <- glm(x.u ~ txt + s + txt*s, data = data.y, family = binomial)
  
  #print(summary(mod3)$coef)
  #print(with(data[data$x.u != -9 & n.f0 > 0 & t != times[1], ], table(s, round(x.u, 1))))
  # 
  # tab <- with(data[data$x.u != -9 & n.f0 > 0 & t != times[1], ], table(s, round(x.u, 1)))
  # 
  # print(round(tab[2, ]/colSums(tab), 3) )

  
  stat    <- summary(mod3)$coef["x.u","z value"]

  
  return(stat)
}
