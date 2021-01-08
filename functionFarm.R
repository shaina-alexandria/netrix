#=========================================================#
#         Function Farm for Data Generation in netRIx     #
#                                                         #
#         Last Updated: 1-8-21                            #
#         Maintainer: Shaina Alexandria                   #
#         Contact: shaina.alexandria22@gmail.com          #
#=========================================================#

# PURPOSE:  sourcing all the functions needed to run a parametric bootstrap for CIs and Fisher's RI tests

#=== sourcing all the reliant functions 

library(matrixStats)
library(abind)
library(nloptr)
library(tidyr)
library(dplyr)
library(Matrix)
library(gtools)
library(scam)
library(ggplot2)

#===== Functions for inference =====#

# OBJECTIVE: likelihood function used for test statistic calculation
source("likMW.R")

# OBJECTIVE: bootstrapped sickness information given a network and baseline information
source("paraboot.R")

# OBJECTIVE: obtain a confidence interval with Bayesian Optimization
source("ci.R")
source("optFun.R")
source("boundSearch.R")

# OBJECTIVE: the test statistic functions
source("testStat1.R")
source("testStat2.R")
source("testStat3.R")
source("testStat4.R")
source("testStat5.R")

# OBJECTIVE: conducting the permutation test
source("permtest.R")

#===== Functions for data generation =====#

# OBJECTIVE: generate baseline sickness and calculate trans.prob
source("genBase.R")

# OBJECTIVE: calculate probabilities of sickness at time 1 for the network and generate sickness indicator 
source("sickNext.R")

# OBJECTIVE: same as sick.friend.summ.gen1, but includes 2 weeks of network info for accurate calculation of stuff
source("sickSummary.R")

# OBJECTIVE: generate sickness trajectory from given info (can be used in paraboot also)
source("genSick.R")

# OBJECTIVE: to generate a data set with id, hall number, treatment, number of friends, sick at time 0, binary outcome y, 
#            and proportional outcome y to pass through to the permutation test
#            using a scale-free model
source("genNet.R")




