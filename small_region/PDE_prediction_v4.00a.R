# 3/26/24 optimization delete advection, death, exp growth and zero.rm
# 4/2/24 training data set
# 4/9/24 logistic betas for gamma, lambda
# 4/24/24 two covariates
# 5/8/24 intercept
# 7/11/24 zero-inflated binomial
# 7/15/24 rate calculation optimizing
# 7/31/24 use source raster and homogenization functions v4.R
# 8/14/24 use source raster and homogenization functions v5.03.R and adjustment
# 9/6/24 use source raster and homogenization functions v6.00.R and adjustment
# 9/6/24 use "DEM_30m" as covariate instead of "DEM_30m"
# 9/10/24 change diffusion function, growth function
# 10/09/24 use standardized dist and dem30m
# 1/28/25 use raster and homogenization functions v6.20.R
# 2/3/25 use raster and homogenization functions v6.30.R
# 2/3/25 stream homogenization
# 2/11/25 add a vector denoting grids involving homogenized stream and adjustment
# 2/12/25 adjustment, if stream.homog==1 or edge.homog ==1, growth.rate.pred=0, gamma.pred=1
# 4/25/25 use functions v6.50.R
# 5/7/25 add homogenized observation
# 8/18/25 use functions v6.70.R
# 8/18/25 N = 300, Epsilon = 15, LoG = 30
# 11/23/25 If mu.pred>1, prediction <- 1. It mu.pred < 0, prediction <- 0.

library(raster)
library(fields)
library(RColorBrewer)
library(rasterVis)
library(mvtnorm)
library(gridExtra)
library(fBasics)
library(coda)
library(tidyverse)
library(egg)

# setwd("/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a")
source("/home/yawenguan/jiaqichen/Rcode/raster and homogenization functions v6.70.R")
dat <- read.csv("/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/dat_8520_small_area_300_by_300_train_and_test.csv")
MCMC.result <- read.csv("/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/MCMC_result_20000_real_data_300_by_300_zibinom_homo.csv")

# N_sp <- 10000
burnin <- 5002:nrow(MCMC.result)
# burnin <- 1

# number of MCMC samples to calculate prediction
N_sample <- 1000

N_iteration <- nrow(MCMC.result)
N_start <- 15002

set.seed(12345)
MCMC_sample <- sample(x = N_start:N_iteration, size = N_sample, replace = FALSE)

# beta1.int_est <- 5
# beta3.int_est <- -2
# betag.int_est <- -1
# beta1.dist_est <- 0
# beta3.dist_est <- 0
# betag.dist_est <- 0
# beta1.dm_est <- 0
# beta3.dm_est <- 0
# betag.dm_est <- 0

beta1.int_est <- mean(MCMC.result[burnin,2])
beta3.int_est <- mean(MCMC.result[burnin,3])
betag.int_est <- mean(MCMC.result[burnin,4])
beta1.dist_est <- mean(MCMC.result[burnin,5])
beta3.dist_est <- mean(MCMC.result[burnin,6])
betag.dist_est <- mean(MCMC.result[burnin,7])
beta1.dm_est <- mean(MCMC.result[burnin,8])
beta3.dm_est <- mean(MCMC.result[burnin,9])
betag.dm_est <- mean(MCMC.result[burnin,10])

sink(file = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/parameters estimation.txt")
cat("beta diffusion intercept: ")
beta1.int_est
cat("beta growth intercept: ")
beta3.int_est
cat("beta gamma intercept: ")
betag.int_est
cat("beta diffusion dist: ")
beta1.dist_est
cat("beta growth dist: ")
beta3.dist_est
cat("beta gamma dist: ")
betag.dist_est
cat("beta diffusion dm: ")
beta1.dm_est
cat("beta growth dm: ")
beta3.dm_est
cat("beta gamma dm: ")
betag.dm_est
sink(file = NULL)

# acceptance rate
sink(file = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/acceptance rate.txt")
cat("beta diffusion intercept: ")
mean(MCMC.result[1:nrow(MCMC.result),11])
cat("beta growth intercept: ")
mean(MCMC.result[1:nrow(MCMC.result),12])
cat("beta gamma intercept: ")
mean(MCMC.result[1:nrow(MCMC.result),13])
cat("beta diffusion dist: ")
mean(MCMC.result[1:nrow(MCMC.result),14])
cat("beta growth dist: ")
mean(MCMC.result[1:nrow(MCMC.result),15])
cat("beta gamma dist: ")
mean(MCMC.result[1:nrow(MCMC.result),16])
cat("beta diffusion dm: ")
mean(MCMC.result[1:nrow(MCMC.result),17])
cat("beta growth dm: ")
mean(MCMC.result[1:nrow(MCMC.result),18])
cat("beta gamma dm: ")
mean(MCMC.result[1:nrow(MCMC.result),19])
sink(file = NULL)

# 95% confidence interval
beta1.int_sort <- sort(MCMC.result[burnin,2])
beta1.int_lb <- beta1.int_sort[floor(length(burnin)*0.025)]
beta1.int_ub <- beta1.int_sort[ceiling(length(burnin)*0.975)]

beta3.int_sort <- sort(MCMC.result[burnin,3])
beta3.int_lb <- beta3.int_sort[floor(length(burnin)*0.025)]
beta3.int_ub <- beta3.int_sort[ceiling(length(burnin)*0.975)]

betag.int_sort <- sort(MCMC.result[burnin,4])
betag.int_lb <- betag.int_sort[floor(length(burnin)*0.025)]
betag.int_ub <- betag.int_sort[ceiling(length(burnin)*0.975)]

beta1.dist_sort <- sort(MCMC.result[burnin,5])
beta1.dist_lb <- beta1.dist_sort[floor(length(burnin)*0.025)]
beta1.dist_ub <- beta1.dist_sort[ceiling(length(burnin)*0.975)]

beta3.dist_sort <- sort(MCMC.result[burnin,6])
beta3.dist_lb <- beta3.dist_sort[floor(length(burnin)*0.025)]
beta3.dist_ub <- beta3.dist_sort[ceiling(length(burnin)*0.975)]

betag.dist_sort <- sort(MCMC.result[burnin,7])
betag.dist_lb <- betag.dist_sort[floor(length(burnin)*0.025)]
betag.dist_ub <- betag.dist_sort[ceiling(length(burnin)*0.975)]

beta1.dm_sort <- sort(MCMC.result[burnin,8])
beta1.dm_lb <- beta1.dm_sort[floor(length(burnin)*0.025)]
beta1.dm_ub <- beta1.dm_sort[ceiling(length(burnin)*0.975)]

beta3.dm_sort <- sort(MCMC.result[burnin,9])
beta3.dm_lb <- beta3.dm_sort[floor(length(burnin)*0.025)]
beta3.dm_ub <- beta3.dm_sort[ceiling(length(burnin)*0.975)]

betag.dm_sort <- sort(MCMC.result[burnin,10])
betag.dm_lb <- betag.dm_sort[floor(length(burnin)*0.025)]
betag.dm_ub <- betag.dm_sort[ceiling(length(burnin)*0.975)]


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

pred_table <- NULL

# index <- c(1:(N*N))
# s1 <- rep(c(1:N), each = N)
# s2 <- rep(c(1:N), N)

# edge.homog
dat <- dat %>% mutate(edge.homog = ifelse(s1 <= EPSILON|
                                      s1 >= (N-EPSILON+1)|
                                      s2 <= EPSILON|
                                      s2 >= (N-EPSILON+1), 1, 0))

# dat <- dat %>% mutate(edge.homog = ifelse(s1 <= EPSILON|
#                                       s1 >= (N-EPSILON+1)|
#                                       s2 <= EPSILON|
#                                       s2 >= (N-EPSILON+1)|
#                                       Stream==1, 1, 0))

# original raster and rescaling raster
onematrix.NbyN <- matrix(rep(1, N*N), nrow = N, ncol = N)
origin.rst <- raster(onematrix.NbyN, xmn = 0, xmx = 1, ymn = 0, ymx = 1, 
                     crs = NA)
rescaling.rst <- aggregate(origin.rst, fact = EPSILON, na.rm = TRUE, 
                           fun = myharmonicmeanfun)

neighbormatrix <- neighborhood(rescaling.rst)

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

# true parameter calculation
covariates <- dat %>% filter(year == 1985) %>% 
  dplyr::select("distance.str","DEM_30m", "s1", "s2")
dist.to.str <- covariates[,"distance.str"]
dem30m <- covariates[,"DEM_30m"]
s10 <- covariates[,"s1"]
s20 <- covariates[,"s2"]

# zerovector <- rep(0, N*N)

# if(DIFFUSION == 1){
diffu.rate.pred <- exp(-(beta1.int + 
                           beta1.dist * dist.to.str + 
                           beta1.dm * dem30m))

# diffu.rate.pred <- exp(beta1.int + 
#                          beta1.dist * dist.to.str + 
#                          beta1.dm * dem30m)/10000
# } else{
#   diffu.rate.pred <- zerovector
# }

# if(GROWTH == 1){
growth.rate.pred <- exp(-(beta3.int + 
                            beta3.dist * dist.to.str + 
                            beta3.dm * dem30m))

# growth.rate.pred <- exp(beta3.int + 
#                           beta3.dist * dist.to.str + 
#                           beta3.dm * dem30m)
# } else{
#   growth.rate.pred <- zerovector
# }

# if(GAMMA == 1){
gamma.pred <- 1/(1+exp(-(betag.int + 
                           betag.dist * dist.to.str + 
                           betag.dm * dem30m)))
# } else{
#   gamma.pred <- zerovector
# }

# edge85.homog
edge85 <- as.vector(dat[dat$year==1985, "edge"])
edge85.homog <- as.vector(dat[dat$year==1985, "edge.homog"])

# stream85
# stream85 <- as.vector(dat[dat$year==1985, "Stream"])

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

# edge85.homog <- dat %>% filter(year==1985) %>% select(edge.homog)

# diffu.rate.pred[which(stream85.homog==1 | edge85.homog == 1)] <- 0
# growth.rate.pred[which(stream85.homog==1 | edge85.homog == 1)] <- 0
growth.rate.pred[which(stream85 ==1 | edge85 == 1)] <- 0
gamma.pred[which(stream85.homog==1 | edge85.homog == 1)] <- 1

# growth.rate.pred[which(stream85==1)] <- 0
# gamma.pred[which(stream85==1)] <- 1

# diffu.rate.pred[which(edge85.homog==1)] <- 0
# growth.rate.pred[which(edge85.homog==1)] <- 0
# gamma.pred[which(edge85.homog==1)] <- 1

# my PDE solution
# mu0 <- dat$trees[dat$year==1985]
# diffu.r <- diffu.rate.pred
# growth.r <- growth.rate.pred

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

# process <- mynumericalPDEhomofun.Rcpp_mat(n = N, epsilon = EPSILON, 
#                                           initial.vec = mu0, 
#                                           diffu.vec = diffu.r, 
#                                           growth.vec = growth.r,
#                                           time = TIME, 
#                                           step = STEP, 
#                                           delta_t = DELTA_T, 
#                                           min.value = 0, max.value = 1)

mu.pred <-  as.vector(process$result.mu)

dat <- cbind(dat, mu.pred)
diffu.prediction <- rep(diffu.rate.pred, TIME+1)
growth.prediction <- rep(growth.rate.pred, TIME+1)
gamma.prediction <- rep(gamma.pred, TIME+1)
dat <- cbind(dat, diffu.prediction, growth.prediction, gamma.prediction)


# prediction 
for (k in 1:N_sample) {
  beta1.int_est <- MCMC.result[MCMC_sample[k],2]
  beta3.int_est <- MCMC.result[MCMC_sample[k],3]
  betag.int_est <- MCMC.result[MCMC_sample[k],4]
  beta1.dist_est <- MCMC.result[MCMC_sample[k],5]
  beta3.dist_est <- MCMC.result[MCMC_sample[k],6]
  betag.dist_est <- MCMC.result[MCMC_sample[k],7]
  beta1.dm_est <- MCMC.result[MCMC_sample[k],8]
  beta3.dm_est <- MCMC.result[MCMC_sample[k],9]
  betag.dm_est <- MCMC.result[MCMC_sample[k],10]
  
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
  
  mu.vec <-  as.vector(process$result.mu)
  
  diffu.vec <- rep(diffu.rate.pred, TIME+1)
  growth.vec <- rep(growth.rate.pred, TIME+1)
  gamma.vec <- rep(gamma.pred, TIME+1)
  
  mean.pred <- mu.vec[which(dat$edge.homog==0 & dat$time!=0)]
  
  if(sum(mean.pred < 0) > 0){
    warning("There is value less than min value.")
  }
  
  if(sum(mean.pred > 1) > 0) {
    warning("There is value greater than max value.")
  }
  
  pred_vec <- rep(0, N*N*(TIME+1))
  
  for (j in 1:nrow(dat)) {
    
    #the edges and values at 0 are not simulated
    if(dat[j, "edge.homog"] == 1 | stream.homog[j] == 1){
      pred_vec[j] <- 0
    } else{
      meansim <- mu.vec[j]
      gamma <- gamma.vec[j]
      if(dat[j, "year"] == 1985){
        pred_vec[j] <- meansim
      } else{
        if(meansim >= 0 & meansim <= 1)
          pred_vec[j] <- myrZIbinom(n = 1, gamma = gamma, p = meansim, s = SIZE)/SIZE
        else{
          if(meansim > 1)
            pred_vec[j] <- 1
          if(meansim < 0)
            pred_vec[j] <- 0
        }
      }
    }
  }
  
  pred_table <- cbind(pred_table, pred_vec)
}

prediction <- rowMeans(pred_table)
dat <- cbind(dat, prediction)
dat <- dat %>% mutate(pred.expectation = mu.pred*(1-gamma.prediction))

# observation homogenization
Diffu.rst <- origin.rst
Diffu.rst[] <- diffu.rate.pred
Diffu.bar <- aggregate(Diffu.rst, fact = EPSILON, na.rm = TRUE,
                       fun = myharmonicmeanfun)

obs.homog <- NULL
for (i in 1:(TIME+1)) {
  obs <- dat$trees[which(dat$year==1984+i)]
  
  obs.rst <- origin.rst
  obs.rst[] <- obs
  
  C0 <- rescaling.rst
  C0[] <- raster::extract(obs.rst * Diffu.rst, SpatialPoints(rescaling.rst))
  
  Cmat <- as.matrix(C0[])
  Cmat <- Cmat * (1-stream.bar.vec)
  C.rst <- rescaling.rst
  C.rst[] <- as.vector(Cmat)
  obs.rst <- disaggregate(C.rst, EPSILON)/1
  obs.rst[] <- obs.rst[]/Diffu.rst[]
  obs.homog <- c(obs.homog, obs.rst[])
}
dat <- cbind(dat, obs.homog)

write.csv(dat, "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/prediction_real_data_300_by_300_homogenization_diffu_logis_growth_ZIB.csv")

pic_limit <- c(0,1)
dat85 <- dat %>% filter(year == 1985)
# png(filename = "observation in 1985.png")
# ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Observed Tree Coverage in 1985", fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal()
# dev.off()

# png(filename = "expectation in 1985.png")
# ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Expected Tree coverage in 1985") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "estimated mu in 1985.png")
# ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle(expression(paste("Estimated ", mu, " in 1985"))) +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# dat90 <- dat %>% filter(year == 1990)
# png(filename = "observation in 1990.png")
# ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Observed Tree coverage in 1990") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "expectation in 1990.png")
# ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Expected Tree coverage in 1990") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "estimated mu in 1990.png")
# ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle(expression(paste("Estimated ", mu, " in 1990"))) +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# dat95 <- dat %>% filter(year == 1995)
# png(filename = "observation in 1995.png")
# ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Observed Tree coverage in 1995") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "expectation in 1995.png")
# ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Expected Tree coverage in 1995") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "estimated mu in 1995.png")
# ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle(expression(paste("Estimated ", mu, " in 1995"))) +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

dat00 <- dat %>% filter(year == 2000)
# png(filename = "observation in 2000.png")
# ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Observed Tree coverage in 2000") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "expectation in 2000.png")
# ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Expected Tree coverage in 2000") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "estimated mu in 2000.png")
# ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle(expression(paste("Estimated ", mu, " in 2000"))) +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# dat05 <- dat %>% filter(year == 2005)
# png(filename = "observation in 2005.png")
# ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Observed Tree coverage in 2005") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "expectation in 2005.png")
# ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Expected Tree coverage in 2005") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "estimated mu in 2005.png")
# ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle(expression(paste("Estimated ", mu, " in 2005"))) +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

dat10 <- dat %>% filter(year == 2010)
# png(filename = "observation in 2010.png")
# ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Observed Tree coverage in 2010") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "expectation in 2010.png")
# ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Expected Tree coverage in 2010") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "estimated mu in 2010.png")
# ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle(expression(paste("Estimated ", mu, " in 2010"))) +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# dat15 <- dat %>% filter(year == 2015)
# png(filename = "observatation in 2015.png")
# ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Observed Tree coverage in 2015") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "expectation in 2015.png")
# ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Expected Tree coverage in 2015") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "estimated mu in 2015.png")
# ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle(expression(paste("Estimated ", mu, " in 2015"))) +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# dat16 <- dat %>% filter(year == 2016)
# png(filename = "observation in 2016.png")
# ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Observed Tree coverage in 2016") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "prediction in 2016.png")
# ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Predicted Tree coverage in 2016") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "estimated mu in 2016.png")
# ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle(expression(paste("Estimated ", mu, " in 2016"))) +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# dat18 <- dat %>% filter(year == 2018)
# png(filename = "observation in 2018.png")
# ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Observed Tree coverage in 2018") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "prediction in 2018.png")
# ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Predicted Tree coverage in 2018") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "estimated mu in 2018.png")
# ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle(expression(paste("Estimated ", mu, " in 2018"))) +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

dat20 <- dat %>% filter(year == 2020)
# png("observation in 2020.png")
# ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Observed Tree coverage in 2020") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "prediction in 2020.png")
# ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle("Predicted Tree coverage in 2020") +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = "estimated mu in 2020.png")
# ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   ggtitle(expression(paste("Estimated ", mu, " in 2020"))) +
#   labs(fill = "tree\ncoverage") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# histrogram
# histdat.train <- dat %>% filter(year >= 1986 & year <= 1985+TRAIN & edge == 0 & Stream == 0)
# histdat.test <- dat %>% filter(year > 1985+TRAIN & edge == 0 & Stream == 0)

# png(filename = "histogram comparison training dataset obs vs expectation.png")
# par(mfrow=c(1,2))
# hist(histdat.train[ ,"trees"], xlab = " ", main = "Distribution of Observed Tree Coverage \n(training)")
# hist(histdat.train[ ,"pred.expectation"], xlab = " ", main = "Distribution of Expected Tree Coverage \n(training)")
# dev.off()

# png(filename = "histogram training dataset obs.png")
# par(mfrow=c(1,1))
# hist(histdat.train[ ,"trees"], xlab = " ", main = "Distribution of Observed Tree Coverage \n(training)")
# dev.off()

# png(filename = "histogram training dataset expectation.png")
# par(mfrow=c(1,1))
# hist(histdat.train[ ,"pred.expectation"], xlab = " ", main = "Distribution of Expected Tree Coverage \n(training)")
# dev.off()

# par(mfrow=c(1,2))
# histdat.train <- dat %>% filter(year >= 1986 & year <= 1985+TRAIN & edge == 0)
# hist(histdat.train[ ,"trees"], ylim = c(0,80000), main = "Histogram of observations \n(training)")
# hist(histdat.train[ ,"expectation"], ylim = c(1,80000), main = "Histogram of expectation \n(training)")

# png(filename = "histogram comparison test dataset obs vs prediction.png")
# par(mfrow=c(1,2))
# hist(histdat.test[ ,"trees"], xlab = " ", main = "Distribution of Observed Tree Coverage \n(test)")
# hist(histdat.test[ ,"pred.expectation"], xlab = " ", main = "Distribution of Predicted Tree Coverage \n(test)")
# dev.off()

# png(filename = "histogram test dataset obs.png")
# par(mfrow=c(1,1))
# hist(histdat.test[ ,"trees"], xlab = " ", main = "Distribution of Observed Tree Coverage \n(test)")
# dev.off()

# png(filename = "histogram test dataset prediction.png")
# par(mfrow=c(1,1))
# hist(histdat.test[ ,"pred.expectation"], xlab = " ", main = "Distribution of Predicted Tree Coverage \n(test)")
# dev.off()

# Confidence interval
sink(file = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/credible interval.txt")
cat("beta diffusion intercept: ")
c(beta1.int_lb, beta1.int_ub)
cat("beta growth intercept: ")
c(beta3.int_lb, beta3.int_ub)
cat("beta gamma intercept: ")
c(betag.int_lb, betag.int_ub)
cat("beta diffusion dist: ")
c(beta1.dist_lb, beta1.dist_ub)
cat("beta growth dist: ")
c(beta3.dist_lb, beta3.dist_ub)
cat("beta gamma dist: ")
c(betag.dist_lb, betag.dist_ub)
cat("beta diffusion dm: ")
c(beta1.dm_lb, beta1.dm_ub)
cat("beta growth dm: ")
c(beta3.dm_lb, beta3.dm_ub)
cat("beta gamma dm: ")
c(betag.dm_lb, betag.dm_ub)
sink(file = NULL)

dat.without.edge.train <- dat %>% filter(edge == 0 & Stream == 0 & type == "train")
dat.without.edge.test <- dat %>% filter(edge == 0 & Stream == 0 & type == "test")
MSE.train <- mean((dat.without.edge.train$trees-dat.without.edge.train$prediction)^2)
MSE.test <- mean((dat.without.edge.test$trees-dat.without.edge.test$prediction)^2, na.rm = TRUE)

sink(file = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/MSE.txt")
cat("MSE training: ")
MSE.train
cat("MSE test: ")
MSE.test
sink(file = NULL)

MSE.train.exp <- mean((dat.without.edge.train$trees-dat.without.edge.train$pred.expectation)^2)
MSE.test.exp <- mean((dat.without.edge.test$trees-dat.without.edge.test$pred.expectation)^2, na.rm = TRUE)

sink(file = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/MSE with expectation.txt")
cat("MSE training: ")
MSE.train.exp
cat("MSE test: ")
MSE.test.exp
sink(file = NULL)

# MSE in terms of homogenization
dat.without.edge.homog.train <- dat %>% filter(type == "train" & edge.homog == 0 & stream.homog == 0)
dat.without.edge.homog.test <- dat %>% filter(type == "test" & edge.homog == 0 & stream.homog == 0)

MSE.train.PDE.homog <- mean((dat.without.edge.homog.train$trees-dat.without.edge.homog.train$prediction)^2, na.rm = TRUE)
MSE.test.PDE.homog <- mean((dat.without.edge.homog.test$trees-dat.without.edge.homog.test$prediction)^2, na.rm = TRUE)

sink(file = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/MSE_homog.txt")
cat("training: year = 1985-2015; test: time = 2016-2020.\n")
cat("MSE PDE training data: ")
MSE.train.PDE.homog
cat("MSE PDE test data: ")
MSE.test.PDE.homog
sink(file = NULL)

MSE.train.PDE.homog.exp <- mean((dat.without.edge.homog.train$trees-dat.without.edge.homog.train$pred.expectation)^2, na.rm = TRUE)
MSE.test.PDE.homog.exp <- mean((dat.without.edge.homog.test$trees-dat.without.edge.homog.test$pred.expectation)^2, na.rm = TRUE)

sink(file = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/MSE homog with expectation.txt")
cat("training: year = 1985-2015; test: time = 2016-2020.\n")
cat("MSE PDE training data: ")
MSE.train.PDE.homog.exp
cat("MSE PDE test data: ")
MSE.test.PDE.homog.exp
sink(file = NULL)

png(filename = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/estimated diffusion.png")
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = diffu.prediction)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  ggtitle("Diffusion Rate Coefficient") +
  labs(fill = "diffusion\nrate coefficient") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/estimated growth.png")
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = growth.prediction)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  ggtitle("Growth Rate Coefficient") +
  labs(fill = "growth\nrate coefficient") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/estimated gamma.png")
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = gamma.prediction)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  labs(title = "Estimated Gamma", fill = "gamma") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/probability0 in 2020.png")
ggplot(data = dat20, aes(s1*LoG, s2*LoG, 
                         fill = gamma.prediction+(1-gamma.prediction)*(1-mu.pred)^100)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  ggtitle("Probability of 0 Tree Coverage in 2020") +
  labs(fill = "probability") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/probability0 in 2000.png")
ggplot(data = dat00, aes(s1*LoG, s2*LoG, 
                         fill = gamma.prediction+(1-gamma.prediction)*(1-mu.pred)^100)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  ggtitle("Probability of 0 Tree Coverage in 2000") +
  labs(fill = "probability") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/probability0 in 2010.png")
ggplot(data = dat10, aes(s1*LoG, s2*LoG, 
                         fill = gamma.prediction+(1-gamma.prediction)*(1-mu.pred)^100)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  ggtitle("Probability of 0 Tree Coverage in 2010") +
  labs(fill = "probability") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "/home/yawenguan/jiaqichen/Region88/paper/real data/PDE/20250903 real data N300 initial a/20251123 prediction a/aggregated stream.png")
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = stream.homog)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  labs(title = "Aggregated Stream") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()