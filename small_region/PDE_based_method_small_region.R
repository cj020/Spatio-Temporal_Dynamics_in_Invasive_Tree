# Parameter estimate of Bayesian Hierarchical Model
# Zero-inflated binomial data model
# Incorporate PDE with ecological diffusion and logistic growth
# Numerical PDE to approximate PDE solution
# Homogenization to reduce computational complexity 
# Multivariate normal parameter model
# MCMC (Metropolis-Hastings algorithm) used to sample from the posterior distribution of the parameter

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
# library(parallel)

setwd("/small_region")
source("functions/raster_and_homogenization_functions.R")

data_dir <- "/small_region/data"
output_dir <- "/small_region/output"

# Dataset loading
dat <- read.csv(file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized.csv"))

# Plot the small region
dat85 <- dat %>% filter(year == 1985)
ggplot(data = dat85, aes(s1, s2, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
  labs(title = "observation in 1985") +
  coord_equal() 

dat20 <- dat %>% filter(year == 2020) %>% mutate(elevation = original.DEM_30m)
ggplot(data = dat20, aes(s1, s2, fill = elevation)) +
  geom_tile() +
  ggtitle("Elevation in the Small Region") +
  theme_bw() +
  scale_fill_gradientn(colours = terrain.colors(7)) + 
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal()

ggplot(data = dat20, aes(s1, s2, fill = Stream)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
  ggtitle("Stream") +
  theme_bw()+
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  coord_equal()

# Parameters setup
# N (number of grids = N * N)
N <- 300
# Rescale ratio (epsilon) in homogenization method
EPSILON <- 15
# N (number of grids = N * N) after rescaling
M <- N/EPSILON
# Length of grid
LoG <- 30
# Time scale of 1 year
TSCALE <- 1
# Number of years
TIME <- 35
# Number of years as training data
TRAIN <- 30
# Time steps within 1 year
STEP <- 5
# Delta t: length of time of 1 step
DELTA_T <- TSCALE/STEP
# Number of covariates
K <- 2
# Size in binomial distribution
SIZE <- 100

# Set up edge label
dat <- dat %>% mutate(edge = ifelse(s1 == 1 |
                                      s1 == N |
                                      s2 == 1|
                                      s2 == N, 1, 0))

# Label training data and test data
dat <- dat %>% mutate(type=ifelse(year<=1985+TRAIN, "train", "test"))

# Write .csv file with training and test label
write.csv(dat, file.path(output_dir, "dat_8520_small_area_300_by_300_train_and_test.csv"))

# Select the columns, setup the training data
dat <- dat %>% select("index", "year", "s1", "s2", "Stream", "DEM_30m", "trees",
                      "distance.str", "edge", "type")
train <- dat %>% filter(year <= 1985+TRAIN)

# To indicate using the PDE with diffusion and growth, zero-inflated binomial distribution as data model
DIFFUSION <- 1
GROWTH <- 1
GAMMA <- 1

# Setup priors of parameters (unknown coefficients such as beta_D, beta_G, beta_gamma)
# Parameters for diffusion coefficient
mu1.int <- 0
Sigma1.int <- 100
mu1.dm <- 0
Sigma1.dm <- 100
mu1.dist <- 0
Sigma1.dist <- 100

# Parameters for growth coefficient
mu3.int <- 0
Sigma3.int <- 100
mu3.dm <- 0
Sigma3.dm <- 100
mu3.dist <- 0
Sigma3.dist <- 100

# Parameter for gamma parameter in zero-inflated binomial distribution
mug.int <- 0
Sigmag.int <- 100
mug.dm <- 0
Sigmag.dm <- 100
mug.dist <- 0
Sigmag.dist <- 100

# Label the edge grid after rescaling in homogenization as edge.homog
train <- train %>% mutate(edge.homog = ifelse(s1 <= EPSILON|
                                                s1 >= (N-EPSILON+1)|
                                                s2 <= EPSILON|
                                                s2 >= (N-EPSILON+1), 1, 0))

# Set up a zero matrix
zerovector <- rep(0, N*N)

# Set up raster-matrix for original scale and rescaling
onematrix.NbyN <- matrix(rep(1, N*N), nrow = N, ncol = N)
origin.rst <- raster(onematrix.NbyN, xmn = 0, xmx = 1, ymn = 0, ymx = 1, 
                     crs = NA)
rescaling.rst <- aggregate(origin.rst, fact = EPSILON, na.rm = TRUE, 
                           fun = myharmonicmeanfun)

# Obtain a matrix of neighbor information for each grid based on rescaling
neighbormatrix <- neighborhood(rescaling.rst)

# Covariate matrix, horizontal and vertical index
covariates <- train %>% filter(year == 1985) %>% 
  select("distance.str", "DEM_30m", "s1", "s2")
dist.to.str <- covariates[,"distance.str"]
dem30m <- covariates[,"DEM_30m"]
s10 <- covariates[,"s1"]
s20 <- covariates[,"s2"]

# Edge index of the first year with original scale and rescaling
edge85.homog <- as.vector(train[train$year==1985, "edge.homog"])
edge85 <- as.vector(train[train$year==1985, "edge"])

# Stream index of the first year with original scale and rescaling
stream85 <- as.vector(train[train$year==1985, "Stream"])
stream.rst <- origin.rst
stream.rst[] <- stream85
stream.bar <- aggregate(stream.rst, fact = EPSILON, na.rm = TRUE, 
                        fun = mean)
stream.bar.vec <- stream.bar[]
stream.bar.vec <- ifelse(stream.bar.vec > 0, 1, 0)

# Stream index of all the training data with rescaling
stream.bar.rst <- rescaling.rst
stream.bar.rst[] <- stream.bar.vec
stream.homog.rst <- disaggregate(stream.bar.rst, EPSILON)
stream85.homog <- stream.homog.rst[]
stream.homog <- rep(stream85.homog, TRAIN+1)

train <- cbind(train, stream.homog)

# Observation in the first year: treein85, as the initial value of the PDE
treein85 <- as.vector(train[train$year==1985, "trees"])


# Loglikelihood function for each step in MCMC
loglikelihood <- function(data, 
                          beta1.int, beta3.int, betag.int, 
                          beta1.dist, beta3.dist, betag.dist, 
                          beta1.dm, beta3.dm, betag.dm){
  
  dataset <- data
  loglikelihood.sum <- 0
  
  # Evalute diffusion coefficients for each grid
  diffu.rate <- exp(-(beta1.int + 
                        beta1.dist * dist.to.str + 
                        beta1.dm * dem30m))
  
  # Evalute growth coefficients for each grid
  growth.rate <- exp(-(beta3.int +
                         beta3.dist * dist.to.str +
                         beta3.dm * dem30m))
  
  # Evalute gamma value for each grid
  gamma.value <- 1/(1+exp(-(betag.int + 
                              betag.dist * dist.to.str + 
                              betag.dm * dem30m)))
  
  # Set the growth coefficient in the stream or in the edge as 0
  growth.rate[which(stream85 == 1 | edge85 == 1)] <- 0
  
  # Set the gamma value in the stream or in the edge, with respect to rescaling, as 1
  gamma.value[which(stream85.homog == 1 | edge85.homog == 1)] <- 1
  
  # Numerical solution of PDE
  process <- mynumericalPDEhomofun.Rcpp_mat(n = N, epsilon = EPSILON, 
                                            initial.vec = treein85, 
                                            stream.vec = stream.bar.vec,
                                            diffu.vec = diffu.rate, 
                                            growth.vec = growth.rate, 
                                            length_of_grid = LoG,
                                            time = TRAIN, step = STEP, delta_t = DELTA_T, 
                                            min.value = 0, max.value = 1)
  
  tree.mean <- as.vector(process$result.mu)
  dataset <- cbind(dataset, tree.mean)
  mean8620 <- as.vector(process$result.mu[,-1])
  
  # Check if there is tree cover greater than 0 or less than 0
  if(sum(mean8620 > 1) > 0 |
     sum(mean8620 < 0) > 0 |
     sum(is.na(mean8620)) > 0){
    loglikelihood.sum <- -Inf
  } else {
    
    # Calculate the likelihood for each grid in each year
    dataset <- dataset %>%
      mutate(loglik = ifelse(edge.homog==1|year==1985|stream.homog==1,
                             0,
                             myZIbinomlikelihood(x = trees * 100,
                                                 gamma = gamma.value,
                                                 p = tree.mean,
                                                 s = SIZE,
                                                 log = TRUE)))
    
    # Loglikelihood for priors
    loglikelihood.beta1.int <- dnorm(beta1.int, mean = mu1.int, sd = Sigma1.int, log = TRUE)
    loglikelihood.beta3.int <- dnorm(beta3.int, mean = mu3.int, sd = Sigma3.int, log = TRUE)
    loglikelihood.betag.int <- dnorm(betag.int, mean = mug.int, sd = Sigmag.int, log = TRUE)
    
    loglikelihood.beta1.dist <- dnorm(beta1.dist, mean = mu1.dist, sd = Sigma1.dist, log = TRUE)
    loglikelihood.beta3.dist <- dnorm(beta3.dist, mean = mu3.dist, sd = Sigma3.dist, log = TRUE)
    loglikelihood.betag.dist <- dnorm(betag.dist, mean = mug.dist, sd = Sigmag.dist, log = TRUE)
    
    loglikelihood.beta1.dm <- dnorm(beta1.dm, mean = mu1.dm, sd = Sigma1.dm, log = TRUE)
    loglikelihood.beta3.dm <- dnorm(beta3.dm, mean = mu3.dm, sd = Sigma3.dm, log = TRUE)
    loglikelihood.betag.dm <- dnorm(betag.dm, mean = mug.dm, sd = Sigmag.dm, log = TRUE)
    
    loglikelihood.sum <- sum(dataset$loglik) +
      loglikelihood.beta1.int +
      loglikelihood.beta3.int +
      loglikelihood.betag.int +
      loglikelihood.beta1.dist +
      loglikelihood.beta3.dist +
      loglikelihood.betag.dist +
      loglikelihood.beta1.dm +
      loglikelihood.beta3.dm +
      loglikelihood.betag.dm
  }
  
  return(loglikelihood.sum)
}


# Acceptance function in the MCMC algorithm
acceptance <- function(data, prev.loglik,
                       beta1.int.pv, beta3.int.pv, betag.int.pv, 
                       beta1.dist.pv, beta3.dist.pv, betag.dist.pv, 
                       beta1.dm.pv, beta3.dm.pv, betag.dm.pv, 
                       beta1.int.nw, beta3.int.nw, betag.int.nw, 
                       beta1.dist.nw, beta3.dist.nw, betag.dist.nw, 
                       beta1.dm.nw, beta3.dm.nw, betag.dm.nw){
  
  # Directly input the loglikelihood from the last sets of parameters as the previous loglikelihood
  loglikelihood.pv <- prev.loglik
  
  # Evaluate loglikelihood with new parameters
  loglikelihood.nw <- loglikelihood(data = data, 
                                    beta1.int = beta1.int.nw, beta3.int = beta3.int.nw,
                                    betag.int = betag.int.nw, 
                                    beta1.dist = beta1.dist.nw, beta3.dist = beta3.dist.nw, 
                                    betag.dist = betag.dist.nw, 
                                    beta1.dm = beta1.dm.nw, beta3.dm = beta3.dm.nw, 
                                    betag.dm = betag.dm.nw)
  
  # Proposal transition
  logTji.beta1.int <- dnorm(beta1.int.pv, mean = beta1.int.nw, sd = sigma.beta1.int, log = TRUE)
  logTji.beta3.int <- dnorm(beta3.int.pv, mean = beta3.int.nw, sd = sigma.beta3.int, log = TRUE)
  logTji.betag.int <- dnorm(betag.int.pv, mean = betag.int.nw, sd = sigma.betag.int, log = TRUE)
  
  logTji.beta1.dist <- dnorm(beta1.dist.pv, mean = beta1.dist.nw, sd = sigma.beta1.dist, log = TRUE)
  logTji.beta3.dist <- dnorm(beta3.dist.pv, mean = beta3.dist.nw, sd = sigma.beta3.dist, log = TRUE)
  logTji.betag.dist <- dnorm(betag.dist.pv, mean = betag.dist.nw, sd = sigma.betag.dist, log = TRUE)
  
  logTji.beta1.dm <- dnorm(beta1.dm.pv, mean = beta1.dm.nw, sd = sigma.beta1.dm, log = TRUE)
  logTji.beta3.dm <- dnorm(beta3.dm.pv, mean = beta3.dm.nw, sd = sigma.beta3.dm, log = TRUE)
  logTji.betag.dm <- dnorm(betag.dm.pv, mean = betag.dm.nw, sd = sigma.betag.dm, log = TRUE)
  
  logTij.beta1.int <- dnorm(beta1.int.nw, mean = beta1.int.pv, sd = sigma.beta1.int, log = TRUE)
  logTij.beta3.int <- dnorm(beta3.int.nw, mean = beta3.int.pv, sd = sigma.beta3.int, log = TRUE)
  logTij.betag.int <- dnorm(betag.int.nw, mean = betag.int.pv, sd = sigma.betag.int, log = TRUE)
  
  logTij.beta1.dist <- dnorm(beta1.dist.nw, mean = beta1.dist.pv, sd = sigma.beta1.dist, log = TRUE)
  logTij.beta3.dist <- dnorm(beta3.dist.nw, mean = beta3.dist.pv, sd = sigma.beta3.dist, log = TRUE)
  logTij.betag.dist <- dnorm(betag.dist.nw, mean = betag.dist.pv, sd = sigma.betag.dist, log = TRUE)
  
  logTij.beta1.dm <- dnorm(beta1.dm.nw, mean = beta1.dm.pv, sd = sigma.beta1.dm, log = TRUE)
  logTij.beta3.dm <- dnorm(beta3.dm.nw, mean = beta3.dm.pv, sd = sigma.beta3.dm, log = TRUE)
  logTij.betag.dm <- dnorm(betag.dm.nw, mean = betag.dm.pv, sd = sigma.betag.dm, log = TRUE)
  
  accpt <- loglikelihood.nw + 
    logTji.beta1.int + logTji.beta3.int + logTji.betag.int + 
    logTji.beta1.dist + logTji.beta3.dist + logTji.betag.dist + 
    logTji.beta1.dm + logTji.beta3.dm + logTji.betag.dm -
    loglikelihood.pv - 
    logTij.beta1.int - logTij.beta3.int - logTij.betag.int - 
    logTij.beta1.dist - logTij.beta3.dist - logTij.betag.dist - 
    logTij.beta1.dm - logTij.beta3.dm - logTij.betag.dm
  
  loglikelihood.vl <- loglikelihood.nw
  
  return(list(accept.vl = accpt,
              loglikelihood.vl = loglikelihood.vl))
}

# Metropolis-Hasting algorithm
# In each iteration, take a draw from the proposal for one parameter (update only one parameter each iteration)
# Calculate the acceptance probability
# Either accept or reject the proposal given the probability r
MetropolisHastingMCMC <- function(theta0, proposal_var_list, iter, burnin){
  
  # Initialize the chain
  chain <- matrix(NA, nrow = iter+1, ncol = length(theta0))
  chain[1, ] <- theta0
  acceptn <- matrix(0, nrow = iter, ncol = length(theta0))
  
  # Evaluate loglikelihood for initial value
  loglikelihood.prev <- loglikelihood(data = train, 
                                      beta1.int = chain[1,1], 
                                      beta3.int = chain[1,2],
                                      betag.int = chain[1,3], 
                                      beta1.dist = chain[1,4], 
                                      beta3.dist = chain[1,5], 
                                      betag.dist = chain[1,6], 
                                      beta1.dm = chain[1,7], 
                                      beta3.dm = chain[1,8], 
                                      betag.dm = chain[1,9])
  
  for (i in 1:iter) {
    
    cat(" Iteration:",i)
    
    # Create a vector for proposed parameters
    theta_star <- rep(0, length(theta0))
    
    # Update beta_diffu_int
    if(DIFFUSION == 1){
      
      # Generate a proposed parameter
      theta_star[1] <- rnorm(1, mean = chain[i,1], sd = proposal_var_list[[1]])
      
      # Evaluate acceptance function, plug in the new (proposed) parameter and old (unchanged) parameter
      accept.result <- acceptance(data = train, 
                                  prev.loglik = loglikelihood.prev,
                                  beta1.int.pv = chain[i,1], 
                                  beta3.int.pv = chain[i,2],
                                  betag.int.pv = chain[i,3],
                                  beta1.dist.pv = chain[i,4], 
                                  beta3.dist.pv = chain[i,5],
                                  betag.dist.pv = chain[i,6],
                                  beta1.dm.pv = chain[i,7], 
                                  beta3.dm.pv = chain[i,8],
                                  betag.dm.pv = chain[i,9],
                                  beta1.int.nw = theta_star[1], 
                                  beta3.int.nw = chain[i,2],
                                  betag.int.nw = chain[i,3],
                                  beta1.dist.nw = chain[i,4], 
                                  beta3.dist.nw = chain[i,5], 
                                  betag.dist.nw = chain[i,6],
                                  beta1.dm.nw = chain[i,7], 
                                  beta3.dm.nw = chain[i,8],
                                  betag.dm.nw = chain[i,9])
      
      # Evaluate acceptance probability r
      r1.int <- accept.result$accept.vl
      
      # Reject the proposal if fail to evaluate loglikehood (for example, loglikelihood = -infinity)
      r1.int <- ifelse(is.nan(r1.int)==1, -Inf, r1.int)
      
      u <- log(runif(1))
      
      # Determine if accept the proposal
      if(u < r1.int){
        chain[i+1,1] <- theta_star[1]
        acceptn[i,1] <- 1
        loglikelihood.prev <- accept.result$loglikelihood.vl
      } else {
        chain[i+1,1] <- chain[i,1]
      }
    } else if (DIFFUSION == 0){
      theta_star[1] <- beta_1.int
      r1.int <- 1
      chain[i+1,1] <- theta_star[1]
      acceptn[i,1] <- 1
    }
    
    cat(" r1.int:", r1.int)
    
    # Updat beta_growth_int
    if(GROWTH == 1){
      theta_star[2] <- rnorm(1, mean = chain[i,2], sd = proposal_var_list[[2]])
      
      accept.result <- acceptance(data = train, 
                                  prev.loglik = loglikelihood.prev,
                                  beta1.int.pv = chain[i+1,1], 
                                  beta3.int.pv = chain[i,2],
                                  betag.int.pv = chain[i,3],
                                  beta1.dist.pv = chain[i,4], 
                                  beta3.dist.pv = chain[i,5],
                                  betag.dist.pv = chain[i,6],
                                  beta1.dm.pv = chain[i,7], 
                                  beta3.dm.pv = chain[i,8],
                                  betag.dm.pv = chain[i,9],
                                  beta1.int.nw = chain[i+1,1], 
                                  beta3.int.nw = theta_star[2],
                                  betag.int.nw = chain[i,3],
                                  beta1.dist.nw = chain[i,4], 
                                  beta3.dist.nw = chain[i,5], 
                                  betag.dist.nw = chain[i,6],
                                  beta1.dm.nw = chain[i,7], 
                                  beta3.dm.nw = chain[i,8],
                                  betag.dm.nw = chain[i,9])
      
      r4.int <- accept.result$accept.vl
      
      r4.int <- ifelse(is.nan(r4.int)==1, -Inf, r4.int)
      
      u <- log(runif(1))
      
      if(u < r4.int){
        chain[i+1,2] <- theta_star[2]
        acceptn[i,2] <- 1
        loglikelihood.prev <- accept.result$loglikelihood.vl
      } else {
        chain[i+1,2] <- chain[i,2]
      }
    } else if (GROWTH == 0){
      theta_star[2] <- beta_3.int
      r4.int <- 1
      chain[i+1,2] <- theta_star[2]
      acceptn[i,2] <- 1
    }
    
    cat(" r4.int:", r4.int)
    
    # Update beta_gamma_int
    if(GAMMA == 1){
      theta_star[3] <- rnorm(1, mean = chain[i,3], sd = proposal_var_list[[3]])
      
      accept.result <- acceptance(data = train, 
                                  prev.loglik = loglikelihood.prev,
                                  beta1.int.pv = chain[i+1,1], 
                                  beta3.int.pv = chain[i+1,2],
                                  betag.int.pv = chain[i,3],
                                  beta1.dist.pv = chain[i,4], 
                                  beta3.dist.pv = chain[i,5],
                                  betag.dist.pv = chain[i,6],
                                  beta1.dm.pv = chain[i,7], 
                                  beta3.dm.pv = chain[i,8],
                                  betag.dm.pv = chain[i,9],
                                  beta1.int.nw = chain[i+1,1], 
                                  beta3.int.nw = chain[i+1,2],
                                  betag.int.nw = theta_star[3],
                                  beta1.dist.nw = chain[i,4], 
                                  beta3.dist.nw = chain[i,5], 
                                  betag.dist.nw = chain[i,6],
                                  beta1.dm.nw = chain[i,7], 
                                  beta3.dm.nw = chain[i,8],
                                  betag.dm.nw = chain[i,9])
      
      r6.int <- accept.result$accept.vl
      
      r6.int <- ifelse(is.nan(r6.int)==1, -Inf, r6.int)
      
      u <- log(runif(1))
      
      if(u < r6.int){
        chain[i+1,3] <- theta_star[3]
        acceptn[i,3] <- 1
        loglikelihood.prev <- accept.result$loglikelihood.vl
      } else {
        chain[i+1,3] <- chain[i,3]
      }
      
    } else if (GAMMA == 0){
      theta_star[3] <- beta_g.int
      r6.int <- 1
      chain[i+1,3] <- theta_star[3]
      acceptn[i,3] <- 1
    }
    
    cat(" r6.int:", r6.int)
    
    # Update beta_diffu_dist
    if(DIFFUSION == 1){
      theta_star[4] <- rnorm(1, mean = chain[i,4], sd = proposal_var_list[[4]])
      
      accept.result <- acceptance(data = train, 
                                  prev.loglik = loglikelihood.prev,
                                  beta1.int.pv = chain[i+1,1], 
                                  beta3.int.pv = chain[i+1,2],
                                  betag.int.pv = chain[i+1,3],
                                  beta1.dist.pv = chain[i,4], 
                                  beta3.dist.pv = chain[i,5],
                                  betag.dist.pv = chain[i,6],
                                  beta1.dm.pv = chain[i,7], 
                                  beta3.dm.pv = chain[i,8],
                                  betag.dm.pv = chain[i,9],
                                  beta1.int.nw = chain[i+1,1], 
                                  beta3.int.nw = chain[i+1,2],
                                  betag.int.nw = chain[i+1,3],
                                  beta1.dist.nw = theta_star[4], 
                                  beta3.dist.nw = chain[i,5], 
                                  betag.dist.nw = chain[i,6],
                                  beta1.dm.nw = chain[i,7], 
                                  beta3.dm.nw = chain[i,8],
                                  betag.dm.nw = chain[i,9])
      
      r1.dist <- accept.result$accept.vl
      
      r1.dist <- ifelse(is.nan(r1.dist)==1, -Inf, r1.dist)
      
      u <- log(runif(1))
      
      if(u < r1.dist){
        chain[i+1,4] <- theta_star[4]
        acceptn[i,4] <- 1
        loglikelihood.prev <- accept.result$loglikelihood.vl
      } else {
        chain[i+1,4] <- chain[i,4]
      }
    } else if (DIFFUSION == 0){
      theta_star[4] <- beta_1.dist
      r1.dist <- 1
      chain[i+1,4] <- theta_star[4]
      acceptn[i,4] <- 1
    }
    
    cat(" r1.dist:", r1.dist)
    
    # Update beta_growth_dist
    if(GROWTH == 1){
      theta_star[5] <- rnorm(1, mean = chain[i,5], sd = proposal_var_list[[5]])
      
      accept.result <- acceptance(data = train, 
                                  prev.loglik = loglikelihood.prev,
                                  beta1.int.pv = chain[i+1,1], 
                                  beta3.int.pv = chain[i+1,2],
                                  betag.int.pv = chain[i+1,3],
                                  beta1.dist.pv = chain[i+1,4], 
                                  beta3.dist.pv = chain[i,5],
                                  betag.dist.pv = chain[i,6],
                                  beta1.dm.pv = chain[i,7], 
                                  beta3.dm.pv = chain[i,8],
                                  betag.dm.pv = chain[i,9],
                                  beta1.int.nw = chain[i+1,1], 
                                  beta3.int.nw = chain[i+1,2],
                                  betag.int.nw = chain[i+1,3],
                                  beta1.dist.nw = chain[i+1,4], 
                                  beta3.dist.nw = theta_star[5], 
                                  betag.dist.nw = chain[i,6],
                                  beta1.dm.nw = chain[i,7], 
                                  beta3.dm.nw = chain[i,8],
                                  betag.dm.nw = chain[i,9])
      
      r4.dist <- accept.result$accept.vl
      
      r4.dist <- ifelse(is.nan(r4.dist)==1, -Inf, r4.dist)
      
      u <- log(runif(1))
      
      if(u < r4.dist){
        chain[i+1,5] <- theta_star[5]
        acceptn[i,5] <- 1
        loglikelihood.prev <- accept.result$loglikelihood.vl
      } else {
        chain[i+1,5] <- chain[i,5]
      }
    } else if (GROWTH == 0){
      theta_star[5] <- beta_3.dist
      r4.dist <- 1
      chain[i+1,5] <- theta_star[5]
      acceptn[i,5] <- 1
    }
    
    cat(" r4.dist:", r4.dist)
    
    # Update beta_gamma_dist
    if(GAMMA == 1){
      theta_star[6] <- rnorm(1, mean = chain[i,6], sd = proposal_var_list[[6]])
      
      accept.result <- acceptance(data = train, 
                                  prev.loglik = loglikelihood.prev,
                                  beta1.int.pv = chain[i+1,1], 
                                  beta3.int.pv = chain[i+1,2],
                                  betag.int.pv = chain[i+1,3],
                                  beta1.dist.pv = chain[i+1,4], 
                                  beta3.dist.pv = chain[i+1,5],
                                  betag.dist.pv = chain[i,6],
                                  beta1.dm.pv = chain[i,7], 
                                  beta3.dm.pv = chain[i,8],
                                  betag.dm.pv = chain[i,9],
                                  beta1.int.nw = chain[i+1,1], 
                                  beta3.int.nw = chain[i+1,2],
                                  betag.int.nw = chain[i+1,3],
                                  beta1.dist.nw = chain[i+1,4], 
                                  beta3.dist.nw = chain[i+1,5], 
                                  betag.dist.nw = theta_star[6],
                                  beta1.dm.nw = chain[i,7], 
                                  beta3.dm.nw = chain[i,8],
                                  betag.dm.nw = chain[i,9])
      
      r6.dist <- accept.result$accept.vl
      
      r6.dist <- ifelse(is.nan(r6.dist)==1, -Inf, r6.dist)
      
      u <- log(runif(1))
      
      if(u < r6.dist){
        chain[i+1,6] <- theta_star[6]
        acceptn[i,6] <- 1
        loglikelihood.prev <- accept.result$loglikelihood.vl
      } else {
        chain[i+1,6] <- chain[i,6]
      }
      
    } else if (GAMMA == 0){
      theta_star[6] <- beta_g.dist
      r6.dist <- 1
      chain[i+1,6] <- theta_star[6]
      acceptn[i,6] <- 1
    }
    
    cat(" r6.dist:", r6.dist)
    
    # Update beta_diffu_elev
    if(DIFFUSION == 1){
      theta_star[7] <- rnorm(1, mean = chain[i,7], sd = proposal_var_list[[7]])
      
      accept.result <- acceptance(data = train, 
                                  prev.loglik = loglikelihood.prev,
                                  beta1.int.pv = chain[i+1,1], 
                                  beta3.int.pv = chain[i+1,2],
                                  betag.int.pv = chain[i+1,3],
                                  beta1.dist.pv = chain[i+1,4], 
                                  beta3.dist.pv = chain[i+1,5],
                                  betag.dist.pv = chain[i+1,6],
                                  beta1.dm.pv = chain[i,7], 
                                  beta3.dm.pv = chain[i,8],
                                  betag.dm.pv = chain[i,9],
                                  beta1.int.nw = chain[i+1,1], 
                                  beta3.int.nw = chain[i+1,2],
                                  betag.int.nw = chain[i+1,3],
                                  beta1.dist.nw = chain[i+1,4], 
                                  beta3.dist.nw = chain[i+1,5], 
                                  betag.dist.nw = chain[i+1,6],
                                  beta1.dm.nw = theta_star[7], 
                                  beta3.dm.nw = chain[i,8],
                                  betag.dm.nw = chain[i,9])
      
      r1.dm <- accept.result$accept.vl
      
      r1.dm <- ifelse(is.nan(r1.dm)==1, -Inf, r1.dm)
      
      u <- log(runif(1))
      
      if(u < r1.dm){
        chain[i+1,7] <- theta_star[7]
        acceptn[i,7] <- 1
        loglikelihood.prev <- accept.result$loglikelihood.vl
      } else {
        chain[i+1,7] <- chain[i,7]
      }
    } else if (DIFFUSION == 0){
      theta_star[7] <- beta_1.dm
      r1.dm <- 1
      chain[i+1,7] <- theta_star[7]
      acceptn[i,7] <- 1
    }
    
    cat(" r1.dm:", r1.dm)
    
    # Update beta_growth_elev
    if(GROWTH == 1){
      theta_star[8] <- rnorm(1, mean = chain[i,8], sd = proposal_var_list[[8]])
      
      accept.result <- acceptance(data = train, 
                                  prev.loglik = loglikelihood.prev,
                                  beta1.int.pv = chain[i+1,1], 
                                  beta3.int.pv = chain[i+1,2],
                                  betag.int.pv = chain[i+1,3],
                                  beta1.dist.pv = chain[i+1,4], 
                                  beta3.dist.pv = chain[i+1,5],
                                  betag.dist.pv = chain[i+1,6],
                                  beta1.dm.pv = chain[i+1,7], 
                                  beta3.dm.pv = chain[i,8],
                                  betag.dm.pv = chain[i,9],
                                  beta1.int.nw = chain[i+1,1], 
                                  beta3.int.nw = chain[i+1,2],
                                  betag.int.nw = chain[i+1,3],
                                  beta1.dist.nw = chain[i+1,4], 
                                  beta3.dist.nw = chain[i+1,5], 
                                  betag.dist.nw = chain[i+1,6],
                                  beta1.dm.nw = chain[i+1,7], 
                                  beta3.dm.nw = theta_star[8],
                                  betag.dm.nw = chain[i,9])
      
      r4.dm <- accept.result$accept.vl
      
      r4.dm <- ifelse(is.nan(r4.dm)==1, -Inf, r4.dm)
      
      u <- log(runif(1))
      
      if(u < r4.dm){
        chain[i+1,8] <- theta_star[8]
        acceptn[i,8] <- 1
        loglikelihood.prev <- accept.result$loglikelihood.vl
      } else {
        chain[i+1,8] <- chain[i,8]
      }
    } else if (GROWTH == 0){
      theta_star[8] <- beta_3.dm
      r4.dm <- 1
      chain[i+1,8] <- theta_star[8]
      acceptn[i,8] <- 1
    }
    
    cat(" r4.dm:", r4.dm)
    
    # Update beta_gamma_elev
    if(GAMMA == 1){
      theta_star[9] <- rnorm(1, mean = chain[i,9], sd = proposal_var_list[[9]])
      
      accept.result <- acceptance(data = train, 
                                  prev.loglik = loglikelihood.prev,
                                  beta1.int.pv = chain[i+1,1], 
                                  beta3.int.pv = chain[i+1,2],
                                  betag.int.pv = chain[i+1,3],
                                  beta1.dist.pv = chain[i+1,4], 
                                  beta3.dist.pv = chain[i+1,5],
                                  betag.dist.pv = chain[i+1,6],
                                  beta1.dm.pv = chain[i+1,7], 
                                  beta3.dm.pv = chain[i+1,8],
                                  betag.dm.pv = chain[i,9],
                                  beta1.int.nw = chain[i+1,1], 
                                  beta3.int.nw = chain[i+1,2],
                                  betag.int.nw = chain[i+1,3],
                                  beta1.dist.nw = chain[i+1,4], 
                                  beta3.dist.nw = chain[i+1,5], 
                                  betag.dist.nw = chain[i+1,6],
                                  beta1.dm.nw = chain[i+1,7], 
                                  beta3.dm.nw = chain[i+1,8],
                                  betag.dm.nw = theta_star[9])
      
      r6.dm <- accept.result$accept.vl
      
      r6.dm <- ifelse(is.nan(r6.dm)==1, -Inf, r6.dm)
      
      u <- log(runif(1))
      
      if(u < r6.dm){
        chain[i+1,9] <- theta_star[9]
        acceptn[i,9] <- 1
        loglikelihood.prev <- accept.result$loglikelihood.vl
      } else {
        chain[i+1,9] <- chain[i,9]
      }
      
    } else if (GAMMA == 0){
      theta_star[9] <- beta_g.dm
      r6.dm <- 1
      chain[i+1,9] <- theta_star[9]
      acceptn[i,9] <- 1
    }
    
    cat(" r6.dm:", r6.dm)
    
  }
  
  # Acceptance rate for each paramter
  cat(" Acceptance rate b1.int:", mean(acceptn[burnin:iter,1]))
  cat(" Acceptance rate b3.int:", mean(acceptn[burnin:iter,2]))
  cat(" Acceptance rate bg.int:", mean(acceptn[burnin:iter,3]))
  cat(" Acceptance rate b1.dist:", mean(acceptn[burnin:iter,4]))
  cat(" Acceptance rate b3.dist:", mean(acceptn[burnin:iter,5]))
  cat(" Acceptance rate bg.dist:", mean(acceptn[burnin:iter,6]))
  cat(" Acceptance rate b1.dm:", mean(acceptn[burnin:iter,7]))
  cat(" Acceptance rate b3.dm:", mean(acceptn[burnin:iter,8]))
  cat(" Acceptance rate bg.dm:", mean(acceptn[burnin:iter,9]))
  
  chain <- cbind(chain, rbind(c(1,1,1,1,1,1,1,1,1), acceptn))
  
  return(chain[burnin:(iter+1),])
}


# Test parameters
# beta_diffu
beta_1.int <- rnorm(n = 1, mean = -5, sd = 0.1)
beta_1.dist <- rnorm(n = 1, mean = 0, sd = 0.1)
beta_1.dm <- rnorm(n = 1, mean = -1, sd = 0.1)

# beta_growth
beta_3.int <- rnorm(n = 1, mean = 3, sd = 0.1)
beta_3.dist <- rnorm(n = 1, mean = 1, sd = 0.1)
beta_3.dm <- rnorm(n = 1, mean = 1, sd = 0.1)

# beta_gamma
beta_g.int <- rnorm(n = 1, mean = -1, sd = 0.1)
beta_g.dist <- rnorm(n = 1, mean = 0, sd = 0.1)
beta_g.dm <- rnorm(n = 1, mean = 2, sd = 0.1)

# Transition parameters
sigma.beta1.int <- 0.025
sigma.beta3.int <- 0.005
sigma.betag.int <- 0.025
sigma.beta1.dist <- 0.025
sigma.beta3.dist <- 0.01
sigma.betag.dist <- 0.025
sigma.beta1.dm <- 0.01
sigma.beta3.dm <- 0.005
sigma.betag.dm <- 0.025

initial <- c(beta_1.int, beta_3.int, beta_g.int,
             beta_1.dist, beta_3.dist, beta_g.dist, 
             beta_1.dm, beta_3.dm, beta_g.dm)
varlist <- list(sigma.beta1.int, sigma.beta3.int, sigma.betag.int, 
                sigma.beta1.dist, sigma.beta3.dist, sigma.betag.dist, 
                sigma.beta1.dm, sigma.beta3.dm, sigma.betag.dm)

# Estimate the running time
system.time(MetropolisHastingMCMC(theta0 = initial,
                                  proposal_var_list = varlist,
                                  iter = 100,
                                  burnin = 1))

# Run 20000 iterations
par.table <- MetropolisHastingMCMC(theta0 = initial,
                                   proposal_var_list = varlist,
                                   iter = 20000,
                                   burnin = 1)

# MCMC trace plots
png(filename = file.path(output_dir, "MCMC_path_20000_real_data_300_by_300_zibinom_homo_diffu.png"))
par(mfrow=c(2,2))
plot(1:dim(par.table)[1], par.table[,1], type = 'l',
     main = expression(paste(beta, "diffu,int")),
     xlab = "iterations", ylab = "value")

plot(1:dim(par.table)[1], par.table[,4], type = 'l',
     main = expression(paste(beta, "diffu,dist")),
     xlab = "iterations", ylab = "value")

plot(1:dim(par.table)[1], par.table[,7], type = 'l',
     main = expression(paste(beta, "diffu,dm")),
     xlab = "iterations", ylab = "value")
dev.off()

png(filename = file.path(output_dir, "MCMC_path_5000_real_data_300_by_300_zibinom_homo_growth.png"))
par(mfrow=c(2,2))
plot(1:dim(par.table)[1], par.table[,2], type = 'l',
     main = expression(paste(beta, "growth,int")),
     xlab = "iterations", ylab = "value")

plot(1:dim(par.table)[1], par.table[,5], type = 'l',
     main = expression(paste(beta, "growth,dist")),
     xlab = "iterations", ylab = "value")

plot(1:dim(par.table)[1], par.table[,8], type = 'l',
     main = expression(paste(beta, "growth,dm")),
     xlab = "iterations", ylab = "value")
dev.off()

png(filename = file.path(output_dir, "MCMC_path_5000_real_data_300_by_300_zibinom_homo_gamma.png"))
par(mfrow=c(2,2))
plot(1:dim(par.table)[1], par.table[,3], type = 'l',
     main = expression(paste(beta[gamma], "int")),
     xlab = "iterations", ylab = "value")

plot(1:dim(par.table)[1], par.table[,6], type = 'l',
     main = expression(paste(beta[gamma], "dist")),
     xlab = "iterations", ylab = "value")

plot(1:dim(par.table)[1], par.table[,9], type = 'l',
     main = expression(paste(beta[gamma], "dm")),
     xlab = "iterations", ylab = "value")
dev.off()

# output the result
output <- par.table
write.csv(output, file.path(output_dir, "MCMC_result_20000_real_data_100_by_100_zibinom_homo.csv"))

