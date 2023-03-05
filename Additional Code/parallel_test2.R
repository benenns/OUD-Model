rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)
library(rBeta2009)
library(parallel)
library(foreach)
library(doParallel)
library(tidyr)

# Set number of cores
#n_cores <- detectCores()
n_cores <- 10
registerDoParallel(n_cores)

foreach(i=4:1, .combine='c') %dopar% {
  Sys.sleep(3 * i)
  i
}

foreach(i=4:1, .combine='c', .inorder=FALSE) %dopar% {
  Sys.sleep(3 * i)
  i
}