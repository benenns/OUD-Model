rm(list = ls()) # to clean the workspace

library(parallel)
library(foreach)
library(doParallel)

type <- if (exists("mcfork", mode="function")) "FORK" else "PSOCK"
cores <- getOption("mc.cores", detectCores())
cl <- makeCluster((cores/2), type = type)
print(cl)
registerDoParallel(cl)
getDoParRegistered()

# Generate data
data <- 1:1e9
data_list <- list("1" = data,
                  "2" = data,
                  "3" = data,
                  "4" = data)
# Single core
time_benchmark <- system.time(
  lapply(data_list, mean)
)

time_foreach <- system.time({
  r <- foreach::foreach(i = 1:length(data_list),
                        .combine = rbind) %dopar% {
                          mean(data_list[[i]])
                        }
})
time_foreach[3]
# Stop cluster to free up resources
stopCluster(cl) # Function doesn't work
getDoParRegistered()

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

unregister_dopar()
getDoParRegistered()
