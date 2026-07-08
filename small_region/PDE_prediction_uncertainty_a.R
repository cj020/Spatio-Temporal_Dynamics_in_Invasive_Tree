# 11/11/25 prediction uncertainty
# 11/16/25 prediction as sample from ZIBinom for 1000 MCMC samples

library(fields)
library(RColorBrewer)
library(mvtnorm)
# library(gridExtra)
library(fBasics)
library(coda)
library(tidyverse)
library(egg)
# library(patchwork)
library(ggplot2)
# library(cowplot)
# library(grid)

# setwd("/home/yawenguan/jiaqichen/Region88/paper/real data/prediction uncertainty")
source("/home/yawenguan/jiaqichen/Rcode/raster and homogenization functions v6.70.R")

dat <- read.csv("/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/dat_8520_small_area_300_by_300_train_and_test.csv")

# number of MCMC samples
N_sample <- 1000
# number of prediction for each MCMC sample
# N_pred <- 1000

MCMC.result <- read.csv("/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/MCMC_result_20000_real_data_300_by_300_zibinom_homo.csv")
N_iteration <- nrow(MCMC.result)
N_burnin <- 15002

set.seed(12345)
MCMC_sample <- sample(x = N_burnin:N_iteration, size = N_sample, replace = FALSE)

# number of grids N by N
N <- 300
# rescale epsilon
EPSILON <- 15
# number of grids after rescaling
M <- N/EPSILON
# Length of grid
LoG <- 30
# time scale of 1 year
TSCALE <- 1
# number of years
TIME <- 35
# number of years as training data
TRAIN <- 30
# Time steps within 1 year
STEP <- 5
# delta t: length of time of 1 step
DELTA_T <- TSCALE/STEP
# number of covariates
K <- 2
# size in binomial
SIZE <- 100

DIFFUSION <- 1
GROWTH <- 1
GAMMA <- 1

rlt_table <- matrix(dat$trees, ncol = 1)

# edge.homog
dat <- dat %>% mutate(edge.homog = ifelse(s1 <= EPSILON|
                                            s1 >= (N-EPSILON+1)|
                                            s2 <= EPSILON|
                                            s2 >= (N-EPSILON+1), 1, 0))

# original raster and rescaling raster
onematrix.NbyN <- matrix(rep(1, N*N), nrow = N, ncol = N)
origin.rst <- raster(onematrix.NbyN, xmn = 0, xmx = 1, ymn = 0, ymx = 1, 
                     crs = NA)
rescaling.rst <- aggregate(origin.rst, fact = EPSILON, na.rm = TRUE, 
                           fun = myharmonicmeanfun)

neighbormatrix <- neighborhood(rescaling.rst)

# covariates
covariates <- dat %>% filter(year == 1985) %>% 
  dplyr::select("distance.str","DEM_30m", "s1", "s2")
dist.to.str <- covariates[,"distance.str"]
dem30m <- covariates[,"DEM_30m"]
s10 <- covariates[,"s1"]
s20 <- covariates[,"s2"]

# edge85.homog
edge85 <- as.vector(dat[dat$year==1985, "edge"])
edge85.homog <- as.vector(dat[dat$year==1985, "edge.homog"])

# stream85 & homogenization
stream85 <- as.vector(dat[dat$year==1985, "Stream"])
stream.rst <- origin.rst
stream.rst[] <- stream85
stream.bar <- aggregate(stream.rst, fact = EPSILON, na.rm = TRUE, 
                        fun = mean)
stream.bar.vec <- stream.bar[]
stream.bar.vec <- ifelse(stream.bar.vec > 0, 1, 0)

# grid that involving homogenized stream
stream.bar.rst <- rescaling.rst
stream.bar.rst[] <- stream.bar.vec
stream.homog.rst <- disaggregate(stream.bar.rst, EPSILON)
stream85.homog <- stream.homog.rst[]
stream.homog <- rep(stream85.homog, TIME+1)

dat <- cbind(dat, stream.homog)

# treein85
treein85 <- as.vector(dat[dat$year==1985, "trees"])

# predction 
for (i in 1:N_sample) {
  beta1.int_est <- MCMC.result[MCMC_sample[i],2]
  beta3.int_est <- MCMC.result[MCMC_sample[i],3]
  betag.int_est <- MCMC.result[MCMC_sample[i],4]
  beta1.dist_est <- MCMC.result[MCMC_sample[i],5]
  beta3.dist_est <- MCMC.result[MCMC_sample[i],6]
  betag.dist_est <- MCMC.result[MCMC_sample[i],7]
  beta1.dm_est <- MCMC.result[MCMC_sample[i],8]
  beta3.dm_est <- MCMC.result[MCMC_sample[i],9]
  betag.dm_est <- MCMC.result[MCMC_sample[i],10]
  
  # true value parameters
  beta1.int <- beta1.int_est
  beta1.dist <- beta1.dist_est
  beta1.dm <- beta1.dm_est
  
  # beta for growth rate
  beta3.int <- beta3.int_est
  beta3.dist <- beta3.dist_est
  beta3.dm <- beta3.dm_est
  
  # beta for gamma
  betag.int <- betag.int_est
  betag.dist <- betag.dist_est
  betag.dm <- betag.dm_est
  
  diffu.rate.pred <- exp(-(beta1.int + 
                             beta1.dist * dist.to.str + 
                             beta1.dm * dem30m))
  
  growth.rate.pred <- exp(-(beta3.int + 
                              beta3.dist * dist.to.str + 
                              beta3.dm * dem30m))
  
  gamma.pred <- 1/(1+exp(-(betag.int + 
                             betag.dist * dist.to.str + 
                             betag.dm * dem30m)))
  
  # diffu.rate.pred[which(stream85.homog==1 | edge85.homog == 1)] <- 0
  # growth.rate.pred[which(stream85.homog==1 | edge85.homog == 1)] <- 0
  growth.rate.pred[which(stream85 ==1 | edge85 == 1)] <- 0
  gamma.pred[which(stream85.homog==1 | edge85.homog == 1)] <- 1
  
  
  # solve PDE numerically with homogenization
  process <- mynumericalPDEhomofun.Rcpp_mat(n = N, epsilon = EPSILON, 
                                            initial.vec = treein85, 
                                            stream.vec = stream.bar.vec,
                                            diffu.vec = diffu.rate.pred, 
                                            growth.vec = growth.rate.pred,
                                            length_of_grid = LoG,
                                            time = TIME, 
                                            step = STEP, 
                                            delta_t = DELTA_T, 
                                            min.value = 0, max.value = 1)
  
  mu.pred <-  as.vector(process$result.mu)
  
  diffu.prediction <- rep(diffu.rate.pred, TIME+1)
  growth.prediction <- rep(growth.rate.pred, TIME+1)
  gamma.prediction <- rep(gamma.pred, TIME+1)
  
  # pred_table <- NULL
  
  pred_vec <- rep(0, N*N*(TIME+1))
  
  for (k in 1:nrow(dat)) {
    
    #the edges and values at 0 are not simulated
    if(dat[k, "edge.homog"] == 1 | dat[k, "stream.homog"] == 1){
      pred_vec[k] <- 0
    } else{
      meansim <- mu.pred[k]
      gamma <- gamma.prediction[k]
      if(dat[k, "year"] == 1985){
        pred_vec[k] <- meansim
      } else{
        if(meansim >= 0 & meansim <= 1){
          pred_vec[k] <- myrZIbinom(n = 1, gamma = gamma, p = meansim, s = SIZE)/SIZE
        }
        else{
          if(meansim > 1)
            pred_vec[k] <- 1
          if(meansim < 0)
            pred_vec[k] <- 0
        }
      }
    }
  }
  
  # for (j in 1:N_pred) {
  #   for (k in 1:nrow(dat)) {
  #     
  #     #the edges and values at 0 are not simulated
  #     if(dat[k, "edge.homog"] == 1 | dat[k, "stream.homog"] == 1){
  #       pred_vec[k] <- 0
  #     } else{
  #       meansim <- mu.pred[k]
  #       gamma <- gamma.prediction[k]
  #       if(dat[k, "time"] == 0){
  #         pred_vec[k] <- meansim
  #       } else{
  #         pred_vec[k] <- myrZIbinom(n = 1, gamma = gamma, p = meansim, s = SIZE)/SIZE
  #       }
  #     }
  #   }
  # }
  
  # pred_table <- cbind(pred_table, pred_vec)
  # 
  # pred_mean <- rowMeans(pred_table)
  # pred_3 <- round(pred_mean, 3)
  
  pred_3 <- round(pred_vec, 3)
  
  rlt_table <- cbind(rlt_table, pred_3)
}

rlt_table <- rlt_table[,-1]
rlt_df <- as.data.frame(rlt_table)

lower <- apply(rlt_table, 1, quantile, 0.025, na.rm = TRUE)
upper <- apply(rlt_table, 1, quantile, 0.975, na.rm = TRUE)

sum_withinCI <- rep(0, nrow(dat))

for (l in 1:nrow(dat)) {
  if(dat$trees[l] >= lower[l] & dat$trees[l] <= upper[l])
    sum_withinCI[l] <- 1
}

dat <- cbind(dat, lower, upper, sum_withinCI)
write.csv(dat, "/home/yawenguan/jiaqichen/Region88/paper/real data/prediction uncertainty/prediction_uncertainty_with_1000_MCMC_result_real_data_300_by_300_a.csv")

sink(file = "/home/yawenguan/jiaqichen/Region88/paper/real data/prediction uncertainty/prediction coverage PDE method a.txt")
cat("prediction coverage in training data: ")
mean(sum_withinCI[1:(N*N*(TRAIN+1))])
cat("prediction coverage in test data: ")
mean(sum_withinCI[(N*N*(TRAIN+1)+1):nrow(rlt_table)])
sink(file = NULL)