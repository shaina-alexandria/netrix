#--- function: stat2
# OBJECTIVE: compute test statistic 2 (T_2 from manuscript) for a given dataset
# INPUTS:   data    : (data.frame) data frame with id, cluster, treatment, time, and outcome information. Of the form of genNet()[[2]]
#           network : (array [n, n, tau+1]) network information; 1 for contact between i,j, 0 otherwise. Of the form of genNet()[[3]]
# OUTPUT:   stat    : the value of the test statistic


stat2 <- function(data, network, ...){

  #----- recalculating x.c based on new treatment assignments
  times     <- sort(unique(data$t))
  
  for(t in times[-1]){
    net0      <- as.matrix(network[ , ,t])
    sick0     <- data$s[data$t == (t-1) ]
    txt       <- data$txt[data$t == (t-1) ]
    n.f.sick0 <- colSums(pmax(sick0*net0, 0))                      # total number of sick friends at 0 time pt
    x.c       <- ifelse(n.f.sick0 != 0, colSums(sick0*txt*pmax(net0, 0))/n.f.sick0, -9)  # proportion of sick friends treated
    
    data$x.c[data$t == t] <- x.c
  }
  
  #----- test statistic
  data.y    <- subset(data, subset = (x.c != -9 & n.f.sick0 > 0 & t != times[1]))
  mod4      <- glm(s ~ x.c, data = data.y, family = binomial)
  # 
  # data.y2    <- subset(data, subset = (n.f.sick0 > 0 & t != times[1]))
  # mod4.1    <- glm(s ~ n.f + x.c, data = data.y2, family = binomial)
  
  stat4.0   <- is.na(mod4$coef["x.c"])
  # print(summary(mod4.1)$coef)
  # print(with(data[data$x.c != -9 & n.f.sick0 > 0 & t != times[1], ], table(s, t)) )
  # tab <- with(data[data$x.c != -9 & n.f.sick0 > 0 & t != times[1], ], table(s, round(x.c, 1)))
  # print(round(tab[2, ]/colSums(tab), 3) )
  

  stat      <- ifelse(stat4.0 == FALSE, summary(mod4)$coef["x.c","z value"], NA) # divide by se? sure

  return(as.numeric(stat))
}
