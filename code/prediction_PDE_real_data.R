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
# 7/16/26 parallelize prediction over MCMC samples (mclapply) based on v4.00d
# 7/17/26 prediction, MSE, predictive coverage with 3 levels (with edge and stream, without egde or stream, without edge.homog or stream.homog)

packages <- c("raster", "fields", "RColorBrewer", "rasterVis", "mvtnorm",
              "gridExtra", "fBasics", "coda", "tidyverse", "egg")
missing_packages <- packages[!packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

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
library(parallel)

# Set the path to the local project directory
working_dir <- "C:/Users/YourName/Documents/Spatio-Temporal_Dynamics_in_Invasive_Tree"

# Set working directory
setwd(working_dir)

# Get the working directory
wd <- getwd()

pde_dir <- file.path(wd, "real_data", "PDE")
prediction_dir <- file.path(pde_dir, "prediction")
plot_dir <- prediction_dir
result_dir <- file.path(wd, "real_data", "MSE_and_predictive_uncertainty")
dir.create(prediction_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(wd, "code", "functions.R"))
dat <- read.csv(file.path(wd, "real_data", "PDE", "dat_8520_small_area_300_by_300_train_and_test.csv"))
MCMC.result <- read.csv(file.path(pde_dir, "MCMC", "MCMC_result_20000_real_data_300_by_300_zibinom_homo_a.csv"))

# N_sp <- 10000
burnin <- 5002:nrow(MCMC.result)
# burnin <- 1

# number of MCMC samples to calculate prediction
N_sample <- 1000

N_iteration <- nrow(MCMC.result)
N_start <- 15002

set.seed(12345)
MCMC_sample <- sample(x = N_start:N_iteration, size = N_sample, replace = FALSE)

# number of cores for parallel prediction over MCMC samples
ncore <- max(1, detectCores() - 1)

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

sink(file = file.path(prediction_dir, "parameters_estimation.txt"))
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
sink(file = file.path(prediction_dir, "acceptance_rate.txt"))
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

diffu.rate.pred <- exp(-(beta1.int + 
                           beta1.dist * dist.to.str + 
                           beta1.dm * dem30m))

growth.rate.pred <- exp(-(beta3.int + 
                            beta3.dist * dist.to.str + 
                            beta3.dm * dem30m))

gamma.pred <- 1/(1+exp(-(betag.int + 
                           betag.dist * dist.to.str + 
                           betag.dm * dem30m)))

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

mu.pred <-  as.vector(process$result.mu)

dat <- cbind(dat, mu.pred)
diffu.prediction <- rep(diffu.rate.pred, TIME+1)
growth.prediction <- rep(growth.rate.pred, TIME+1)
gamma.prediction <- rep(gamma.pred, TIME+1)
dat <- cbind(dat, diffu.prediction, growth.prediction, gamma.prediction)

# MSE and Predictive Coverage
observation <- dat$trees
observation_train <- observation[which(dat$type=="train" )]
observation_test <- observation[which(dat$type=="test" )]
observation_train_wt_edge_stream <- observation[which(dat$type=="train" &
                                                        dat$Stream == 0 &
                                                        dat$edge == 0)]
observation_test_wt_edge_stream <- observation[which(dat$type=="test" &
                                                       dat$Stream == 0 &
                                                       dat$edge == 0)]
observation_train_wt_homo_edge_stream <- observation[which(dat$type=="train" &
                                                             dat$stream.homog == 0 &
                                                             dat$edge.homog == 0)]
observation_test_wt_homo_edge_stream <- observation[which(dat$type=="test" &
                                                            dat$stream.homog == 0 &
                                                            dat$edge.homog == 0)]

MSE_per_MCMC_sample_train <- rep(0, N_sample)
MSE_per_MCMC_sample_test <- rep(0, N_sample)
MSE_wt_edge_stream_per_MCMC_sample_train <- rep(0, N_sample)
MSE_wt_edge_stream_per_MCMC_sample_test <- rep(0, N_sample)
MSE_homo_per_MCMC_sample_train <- rep(0, N_sample)
MSE_homo_per_MCMC_sample_test <- rep(0, N_sample)

ZIB_random_draw_per_MCMC_sample <- NULL

# parallel over MCMC samples (Linux/Unix: fork via mclapply)
mcmc_pred_results <- mclapply(1:N_sample, function(k) {
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
  
  expectation <- mu.vec*(1-gamma.vec)
  expectation_train <- expectation[which(dat$type=="train")]
  expectation_test <- expectation[which(dat$type=="test")]
  expectation_train_wt_edge_stream <- expectation[which(dat$type=="train" &
                                                          dat$Stream == 0 &
                                                          dat$edge == 0)]
  expectation_test_wt_edge_stream <- expectation[which(dat$type=="test" &
                                                         dat$Stream == 0 &
                                                         dat$edge == 0)]
  expectation_train_wt_homo_edge_stream <- expectation[which(dat$type=="train" &
                                                               dat$stream.homog == 0 &
                                                               dat$edge.homog == 0)]
  expectation_test_wt_homo_edge_stream <- expectation[which(dat$type=="test" &
                                                              dat$stream.homog == 0 &
                                                              dat$edge.homog == 0)]
  
  MSE_train_k <- mean((observation_train-expectation_train)^2)
  MSE_test_k <- mean((observation_test-expectation_test)^2)
  MSE_wt_edge_stream_train_k <- mean((observation_train_wt_edge_stream-expectation_train_wt_edge_stream)^2)
  MSE_wt_edge_stream_test_k <- mean((observation_test_wt_edge_stream-expectation_test_wt_edge_stream)^2)
  MSE_homo_train_k <- mean((observation_train_wt_homo_edge_stream-expectation_train_wt_homo_edge_stream)^2)
  MSE_homo_test_k <- mean((observation_test_wt_homo_edge_stream-expectation_test_wt_homo_edge_stream)^2)
  
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
  
  list(MSE_train = MSE_train_k,
       MSE_test = MSE_test_k,
       MSE_wt_edge_stream_train = MSE_wt_edge_stream_train_k,
       MSE_wt_edge_stream_test = MSE_wt_edge_stream_test_k,
       MSE_homo_train = MSE_homo_train_k,
       MSE_homo_test = MSE_homo_test_k,
       pred_vec = pred_vec)
}, mc.cores = ncore)

MSE_per_MCMC_sample_train <- sapply(mcmc_pred_results, `[[`, "MSE_train")
MSE_per_MCMC_sample_test <- sapply(mcmc_pred_results, `[[`, "MSE_test")
MSE_wt_edge_stream_per_MCMC_sample_train <- sapply(mcmc_pred_results, `[[`, "MSE_wt_edge_stream_train")
MSE_wt_edge_stream_per_MCMC_sample_test <- sapply(mcmc_pred_results, `[[`, "MSE_wt_edge_stream_test")
MSE_homo_per_MCMC_sample_train <- sapply(mcmc_pred_results, `[[`, "MSE_homo_train")
MSE_homo_per_MCMC_sample_test <- sapply(mcmc_pred_results, `[[`, "MSE_homo_test")
ZIB_random_draw_per_MCMC_sample <- do.call(cbind, lapply(mcmc_pred_results, `[[`, "pred_vec"))
rm(mcmc_pred_results)

write.csv(ZIB_random_draw_per_MCMC_sample, file.path(prediction_dir, "ZIB_random_draws_1000.csv"))

prediction <- rowMeans(ZIB_random_draw_per_MCMC_sample)
dat <- cbind(dat, prediction)
pred.expectation <- mu.pred*(1-gamma.prediction)
dat <- dat %>% mutate(pred.expectation = pred.expectation)

pred.expectation_train <- pred.expectation[which(dat$type=="train")]
pred.expectation_test <- pred.expectation[which(dat$type=="test")]
pred.expectation_train_wt_edge_stream <- pred.expectation[which(dat$type=="train" &
                                                                  dat$Stream == 0 &
                                                                  dat$edge == 0)]
pred.expectation_test_wt_edge_stream <- pred.expectation[which(dat$type=="test" &
                                                                 dat$Stream == 0 &
                                                                 dat$edge == 0)]
pred.expectation_train_wt_homo_edge_stream <- pred.expectation[which(dat$type=="train" &
                                                                       dat$stream.homog == 0 &
                                                                       dat$edge.homog == 0)]
pred.expectation_test_wt_homo_edge_stream <- pred.expectation[which(dat$type=="test" &
                                                                      dat$stream.homog == 0 &
                                                                      dat$edge.homog == 0)]

MSE_pred.expectation_train <- mean((observation_train-pred.expectation_train)^2)
MSE_pred.expectation_test <- mean((observation_test-pred.expectation_test)^2)
MSE_pred.expectation_wt_edge_stream_train <- mean((observation_train_wt_edge_stream-pred.expectation_train_wt_edge_stream)^2)
MSE_pred.expectation_wt_edge_stream_test <- mean((observation_test_wt_edge_stream-pred.expectation_test_wt_edge_stream)^2)
MSE_homo_pred.expectation_train <- mean((observation_train_wt_homo_edge_stream-pred.expectation_train_wt_homo_edge_stream)^2)
MSE_homo_pred.expectation_test <- mean((observation_test_wt_homo_edge_stream-pred.expectation_test_wt_homo_edge_stream)^2)

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

write.csv(dat, file.path(prediction_dir, "prediction_real_data_300_by_300_homogenization_diffu_logis_growth_ZIB.csv"))

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
# png(filename = "observation in 2020.png")
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
sink(file = file.path(prediction_dir, "credible_interval.txt"))
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

MSE_train <- mean(MSE_per_MCMC_sample_train)
q_train <- quantile(MSE_per_MCMC_sample_train, probs = c(0.025, 0.975))
MSE_train_lower <- q_train[1]
MSE_train_upper <- q_train[2]

MSE_test <- mean(MSE_per_MCMC_sample_test)
q_test <- quantile(MSE_per_MCMC_sample_test, probs = c(0.025, 0.975))
MSE_test_lower <- q_test[1]
MSE_test_upper <- q_test[2]

MSE_wt_edge_stream_train <- mean(MSE_wt_edge_stream_per_MCMC_sample_train)
q_wt_edge_stream_train <- quantile(MSE_wt_edge_stream_per_MCMC_sample_train, probs = c(0.025, 0.975))
MSE_wt_edge_stream_train_lower <- q_wt_edge_stream_train[1]
MSE_wt_edge_stream_train_upper <- q_wt_edge_stream_train[2]

MSE_wt_edge_stream_test <- mean(MSE_wt_edge_stream_per_MCMC_sample_test)
q_wt_edge_stream_test <- quantile(MSE_wt_edge_stream_per_MCMC_sample_test, probs = c(0.025, 0.975))
MSE_wt_edge_stream_test_lower <- q_wt_edge_stream_test[1]
MSE_wt_edge_stream_test_upper <- q_wt_edge_stream_test[2]

MSE_homo_train <- mean(MSE_homo_per_MCMC_sample_train)
q_homo_train <- quantile(MSE_homo_per_MCMC_sample_train, probs = c(0.025, 0.975))
MSE_homo_train_lower <- q_homo_train[1]
MSE_homo_train_upper <- q_homo_train[2]

MSE_homo_test <- mean(MSE_homo_per_MCMC_sample_test)
q_homo_test <- quantile(MSE_homo_per_MCMC_sample_test, probs = c(0.025, 0.975))
MSE_homo_test_lower <- q_homo_test[1]
MSE_homo_test_upper <- q_homo_test[2]

sink(file = file.path(result_dir, "MSE_PDE_method.txt"))
cat("MSE with pred.expectation train: ")
MSE_pred.expectation_train
cat("\nMSE with pred.expectation test: ")
MSE_pred.expectation_test
cat("\nMSE with pred.expectation train (excluding edge and stream): ")
MSE_pred.expectation_wt_edge_stream_train
cat("\nMSE with pred.expectation test (excluding edge and stream): ")
MSE_pred.expectation_wt_edge_stream_test
cat("\nMSE with pred.expectation train (excluding edge.homog and stream.homog): ")
MSE_homo_pred.expectation_train
cat("\nMSE with pred.expectation test (excluding edge.homog and stream.homog): ")
MSE_homo_pred.expectation_test
cat("\nMSE training: ")
MSE_train
cat("\nlb: ")
MSE_train_lower
cat("ub: ")
MSE_train_upper
cat("\nMSE test: ")
MSE_test
cat("\nlb: ")
MSE_test_lower
cat("ub: ")
MSE_test_upper
cat("\nMSE train (excluding edge and stream): ")
MSE_wt_edge_stream_train
cat("\nlb: ")
MSE_wt_edge_stream_train_lower
cat("ub: ")
MSE_wt_edge_stream_train_upper
cat("\nMSE test (excluding edge and stream): ")
MSE_wt_edge_stream_test
cat("\nlb: ")
MSE_wt_edge_stream_test_lower
cat("ub: ")
MSE_wt_edge_stream_test_upper
cat("\nMSE train (excluding edge.homog and stream.homog): ")
MSE_homo_train
cat("\nlb: ")
MSE_homo_train_lower
cat("ub: ")
MSE_homo_train_upper
cat("\nMSE test (excluding edge.homog and stream.homog): ")
MSE_homo_test
cat("\nlb: ")
MSE_homo_test_lower
cat("ub: ")
MSE_homo_test_upper
sink(file = NULL)

q_pred <-  t(apply(ZIB_random_draw_per_MCMC_sample, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE))

sum_withinCI <- rep(0, nrow(dat))
lower <- q_pred[ ,1]
upper <- q_pred[ ,2]
lower_3 <- round(lower, 3)
upper_3 <- round(upper, 3)
lengthCI <- upper - lower

sum_withinCI <- as.integer(observation >= lower_3 & observation <= upper_3)

pred_coverage_train <- mean(sum_withinCI[which(dat$type=="train")])
pred_coverage_test <- mean(sum_withinCI[which(dat$type=="test")])

pred_coverage_wt_edge_stream_train <- mean(sum_withinCI[which(dat$type=="train" &
                                                 dat$Stream == 0 &
                                                 dat$edge == 0)])
pred_coverage_wt_edge_stream_test <- mean(sum_withinCI[which(dat$type=="test" &
                                                dat$Stream == 0 &
                                                dat$edge == 0)])

pred_coverage_homo_train <- mean(sum_withinCI[which(dat$type=="train" &
                                                      dat$stream.homog == 0 &
                                                      dat$edge.homog == 0)])
pred_coverage_homo_test <- mean(sum_withinCI[which(dat$type=="test" &
                                                     dat$stream.homog == 0 &
                                                     dat$edge.homog == 0)])

lengthCI_train <- mean(lengthCI[which(dat$type=="train")])
lengthCI_test <- mean(lengthCI[which(dat$type=="test")])

lengthCI_wt_edge_stream_train <- mean(lengthCI[which(dat$type=="train" &
                                        dat$Stream == 0 &
                                        dat$edge == 0)])
lengthCI_wt_edge_stream_test <- mean(lengthCI[which(dat$type=="test" &
                                       dat$Stream == 0 &
                                       dat$edge == 0)])

lengthCI_homo_train <- mean(lengthCI[which(dat$type=="train" &
                                             dat$stream.homog == 0 &
                                             dat$edge.homog == 0)])
lengthCI_homo_test <- mean(lengthCI[which(dat$type=="test" &
                                            dat$stream.homog == 0 &
                                            dat$edge.homog == 0)])


sink(file = file.path(result_dir, "prediction_coverage_PDE_method.txt"))
cat("predictive coverage training: ")
pred_coverage_train
cat("\naverage length: ")
lengthCI_train
cat("\npredictive coverage test: ")
pred_coverage_test
cat("\naverage length: ")
lengthCI_test
cat("\npredictive coverage training (excluding edge and stream): ")
pred_coverage_wt_edge_stream_train
cat("\naverage length: ")
lengthCI_wt_edge_stream_train
cat("\npredictive coverage test (excluding edge and stream): ")
pred_coverage_wt_edge_stream_test
cat("\naverage length: ")
lengthCI_wt_edge_stream_test
cat("\npredictive coverage train (excluding edge.homog and stream.homog): ")
pred_coverage_homo_train
cat("\naverage length: ")
lengthCI_homo_train
cat("\npredictive coverage test (excluding edge.homog and stream.homog): ")
pred_coverage_homo_test
cat("\naverage length: ")
lengthCI_homo_test
sink(file = NULL)

png(filename = file.path(plot_dir, "estimated_diffusion.png"))
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

png(filename = file.path(plot_dir, "estimated_growth.png"))
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

png(filename = file.path(plot_dir, "estimated_gamma.png"))
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

png(filename = file.path(plot_dir, "probability0_in_2020.png"))
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

png(filename = file.path(plot_dir, "probability0_in_2000.png"))
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

png(filename = file.path(plot_dir, "probability0_in_2010.png"))
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

png(filename = file.path(plot_dir, "aggregated_stream.png"))
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