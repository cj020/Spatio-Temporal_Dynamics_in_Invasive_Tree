# Estimation and prediction with Bayesian Hierarchical Model
# Based on parameter posterior estimate with MCMC algorithm
# Burn in the first 5000 samples
# Zero-inflated binomial data model
# Incorporate PDE with ecological diffusion and logistic growth
# Numerical PDE to approximate PDE solution
# Homogenization to reduce computational complexity 
# Multivariate normal parameter model

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

setwd("/small_region")

source("functions/raster_and_homogenization_functions.R")

data_dir <- "/small_region/data"
output_dir <- "/small_region/output/PDE"
result_dir <- "/small_region/result/PDE"

dat <- read.csv(file.path(output_dir, "dat_8520_small_area_300_by_300_train_and_test.csv"))
result.table <- read.csv(file.path(output_dir, "MCMC_result_20000_real_data_300_by_300_zibinom_homo.csv"))

# Burnin the first 5000 samples
burnin <- 5002:nrow(result.table)

# Parameter estimate after burnin
beta1.int_est <- mean(result.table[burnin,2])
beta3.int_est <- mean(result.table[burnin,3])
betag.int_est <- mean(result.table[burnin,4])
beta1.dist_est <- mean(result.table[burnin,5])
beta3.dist_est <- mean(result.table[burnin,6])
betag.dist_est <- mean(result.table[burnin,7])
beta1.dm_est <- mean(result.table[burnin,8])
beta3.dm_est <- mean(result.table[burnin,9])
betag.dm_est <- mean(result.table[burnin,10])

# Save parameter estimate
sink(file = file.path(result_dir, "parameters estimation.txt"))
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

# Save acceptance rate
sink(file = file = file.path(result_dir, "acceptance rate.txt"))
cat("beta diffusion intercept: ")
mean(result.table[1:nrow(result.table),11])
cat("beta growth intercept: ")
mean(result.table[1:nrow(result.table),12])
cat("beta gamma intercept: ")
mean(result.table[1:nrow(result.table),13])
cat("beta diffusion dist: ")
mean(result.table[1:nrow(result.table),14])
cat("beta growth dist: ")
mean(result.table[1:nrow(result.table),15])
cat("beta gamma dist: ")
mean(result.table[1:nrow(result.table),16])
cat("beta diffusion dm: ")
mean(result.table[1:nrow(result.table),17])
cat("beta growth dm: ")
mean(result.table[1:nrow(result.table),18])
cat("beta gamma dm: ")
mean(result.table[1:nrow(result.table),19])
sink(file = NULL)

# Calculate 95% confidence interval
beta1.int_sort <- sort(result.table[burnin,2])
beta1.int_lb <- beta1.int_sort[floor(length(burnin)*0.025)]
beta1.int_ub <- beta1.int_sort[ceiling(length(burnin)*0.975)]

beta3.int_sort <- sort(result.table[burnin,3])
beta3.int_lb <- beta3.int_sort[floor(length(burnin)*0.025)]
beta3.int_ub <- beta3.int_sort[ceiling(length(burnin)*0.975)]

betag.int_sort <- sort(result.table[burnin,4])
betag.int_lb <- betag.int_sort[floor(length(burnin)*0.025)]
betag.int_ub <- betag.int_sort[ceiling(length(burnin)*0.975)]

beta1.dist_sort <- sort(result.table[burnin,5])
beta1.dist_lb <- beta1.dist_sort[floor(length(burnin)*0.025)]
beta1.dist_ub <- beta1.dist_sort[ceiling(length(burnin)*0.975)]

beta3.dist_sort <- sort(result.table[burnin,6])
beta3.dist_lb <- beta3.dist_sort[floor(length(burnin)*0.025)]
beta3.dist_ub <- beta3.dist_sort[ceiling(length(burnin)*0.975)]

betag.dist_sort <- sort(result.table[burnin,7])
betag.dist_lb <- betag.dist_sort[floor(length(burnin)*0.025)]
betag.dist_ub <- betag.dist_sort[ceiling(length(burnin)*0.975)]

beta1.dm_sort <- sort(result.table[burnin,8])
beta1.dm_lb <- beta1.dm_sort[floor(length(burnin)*0.025)]
beta1.dm_ub <- beta1.dm_sort[ceiling(length(burnin)*0.975)]

beta3.dm_sort <- sort(result.table[burnin,9])
beta3.dm_lb <- beta3.dm_sort[floor(length(burnin)*0.025)]
beta3.dm_ub <- beta3.dm_sort[ceiling(length(burnin)*0.975)]

betag.dm_sort <- sort(result.table[burnin,10])
betag.dm_lb <- betag.dm_sort[floor(length(burnin)*0.025)]
betag.dm_ub <- betag.dm_sort[ceiling(length(burnin)*0.975)]

#Parameter setup
# N (number of grids = N * N)
N <- 300
# Rescale epsilon
EPSILON <- 15
# Number of grids after rescaling
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
# Size in binomial
SIZE <- 100

# To indicate using the PDE with diffusion and growth, zero-inflated binomial distribution as data model
DIFFUSION <- 1
GROWTH <- 1
GAMMA <- 1

# Label the edge grid after rescaling in homogenization as edge.homog
dat <- dat %>% mutate(edge.homog = ifelse(s1 <= EPSILON|
                                      s1 >= (N-EPSILON+1)|
                                      s2 <= EPSILON|
                                      s2 >= (N-EPSILON+1), 1, 0))

# Set up raster-matrix for original scale and rescaling
onematrix.NbyN <- matrix(rep(1, N*N), nrow = N, ncol = N)
origin.rst <- raster(onematrix.NbyN, xmn = 0, xmx = 1, ymn = 0, ymx = 1, 
                     crs = NA)
rescaling.rst <- aggregate(origin.rst, fact = EPSILON, na.rm = TRUE, 
                           fun = myharmonicmeanfun)

# Obtain a matrix of neighbor information for each grid based on rescaling
neighbormatrix <- neighborhood(rescaling.rst)

# Set MCMC estimate as true parameters value
# Beta for diffusion coefficients
beta1.int <- beta1.int_est
beta1.dist <- beta1.dist_est
beta1.dm <- beta1.dm_est

# Beta for growth coefficients
beta3.int <- beta3.int_est
beta3.dist <- beta3.dist_est
beta3.dm <- beta3.dm_est

# Beta for gamma coefficients
betag.int <- betag.int_est
betag.dist <- betag.dist_est
betag.dm <- betag.dm_est

# Parameter calculation
covariates <- dat %>% filter(year == 1985) %>% 
  dplyr::select("distance.str","DEM_30m", "s1", "s2")
dist.to.str <- covariates[,"distance.str"]
dem30m <- covariates[,"DEM_30m"]
s10 <- covariates[,"s1"]
s20 <- covariates[,"s2"]

diffu.rate.pred <- exp(-(beta1.int + 
                           beta1.dist * dist.to.str + 
                           beta1.dm * dem30m))

growth.rate.pred <- exp(-(beta3.int + 
                            beta3.dist * dist.to.str + 
                            beta3.dm * dem30m))

gamma.pred <- 1/(1+exp(-(betag.int + 
                           betag.dist * dist.to.str + 
                           betag.dm * dem30m)))

# Edge index of the first year with original scale and rescaling
edge85 <- as.vector(dat[dat$year==1985, "edge"])
edge85.homog <- as.vector(dat[dat$year==1985, "edge.homog"])

# Stream index of the first year with original scale and rescaling
stream85 <- as.vector(dat[dat$year==1985, "Stream"])
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
stream.homog <- rep(stream85.homog, TIME+1)

dat <- cbind(dat, stream.homog)

# Observation in the first year: treein85, as the initial value of the PDE
treein85 <- as.vector(dat[dat$year==1985, "trees"])

# Set the growth coefficient in the stream or in the edge as 0
growth.rate.pred[which(stream85 ==1 | edge85 == 1)] <- 0

# Set the gamma value in the stream or in the edge, with respect to rescaling, as 1
gamma.pred[which(stream85.homog==1 | edge85.homog == 1)] <- 1

# Numerical solution of PDE
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

# Check if there is tree cover greater than 0 or less than 0
mean8620 <- dat[which(dat$edge.homog==0 & dat$year!=1985),"mu.pred"]

if(sum(mean8620 < 0) > 0){
  warning("There is value less than min value.")
}

if(sum(mean8620 > 1) > 0) {
  warning("There is value greater than max value.")
}

# Estimate/Predict the tree cover with zero-inflated binomial distribution
set.seed(1234)

prediction <- rep(0, N*N*(TIME+1))

for (i in 1:nrow(dat)) {
  
  # The grid with edge.homogs=1, stream.edge=1 or initial values in 1985 do not predict
  if(dat[i, "edge.homog"] == 1 | dat[i, "stream.homog"] == 1){
    prediction[i] <- 0
  } else{
    meansim <- dat[i, "mu.pred"]
    gamma <- dat[i, "gamma.prediction"]
    
    if(dat[i, "year"] == 1985){
      prediction[i] <- meansim
    } else{
      
      if(meansim >= 0 & meansim <= 1)
        prediction[i] <- myrZIbinom(n = 1, gamma = gamma, p = meansim, s = SIZE)/SIZE
      else{
        
        if(meansim > 1)
          prediction[i] <- 1
        if(meansim < 0)
          prediction[i] <- 0
      }
    }
  }
}

dat <- cbind(dat, prediction)
dat <- dat %>% mutate(pred.expectation = mu.pred*(1-gamma.prediction))

# Homogenizing the observation for comparison
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

# Save the result
write.csv(dat, file.path(result_dir, "prediction_real_data_300_by_300_homogenization_diffu_logis_growth_ZIB.csv"))

# Graph check
pic_limit <- c(0,1)
dat85 <- dat %>% filter(year == 1985)
png(filename = "observation in 1985.png")
ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Observed Tree Coverage in 1985", fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal()
dev.off()

png(filename = "expectation in 1985.png")
ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Expected Tree coverage in 1985") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "estimated mu in 1985.png")
ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle(expression(paste("Estimated ", mu, " in 1985"))) +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat90 <- dat %>% filter(year == 1990)
png(filename = "observation in 1990.png")
ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Observed Tree coverage in 1990") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "expectation in 1990.png")
ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Expected Tree coverage in 1990") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "estimated mu in 1990.png")
ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle(expression(paste("Estimated ", mu, " in 1990"))) +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat95 <- dat %>% filter(year == 1995)
png(filename = "observation in 1995.png")
ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Observed Tree coverage in 1995") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "expectation in 1995.png")
ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Expected Tree coverage in 1995") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "estimated mu in 1995.png")
ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle(expression(paste("Estimated ", mu, " in 1995"))) +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat00 <- dat %>% filter(year == 2000)
png(filename = "observation in 2000.png")
ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Observed Tree coverage in 2000") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "expectation in 2000.png")
ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Expected Tree coverage in 2000") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "estimated mu in 2000.png")
ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle(expression(paste("Estimated ", mu, " in 2000"))) +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat05 <- dat %>% filter(year == 2005)
png(filename = "observation in 2005.png")
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Observed Tree coverage in 2005") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "expectation in 2005.png")
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Expected Tree coverage in 2005") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "estimated mu in 2005.png")
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle(expression(paste("Estimated ", mu, " in 2005"))) +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat10 <- dat %>% filter(year == 2010)
png(filename = "observation in 2010.png")
ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Observed Tree coverage in 2010") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "expectation in 2010.png")
ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Expected Tree coverage in 2010") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "estimated mu in 2010.png")
ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle(expression(paste("Estimated ", mu, " in 2010"))) +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat15 <- dat %>% filter(year == 2015)
png(filename = "observatation in 2015.png")
ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Observed Tree coverage in 2015") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "expectation in 2015.png")
ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Expected Tree coverage in 2015") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "estimated mu in 2015.png")
ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle(expression(paste("Estimated ", mu, " in 2015"))) +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat16 <- dat %>% filter(year == 2016)
png(filename = "observation in 2016.png")
ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Observed Tree coverage in 2016") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "prediction in 2016.png")
ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Predicted Tree coverage in 2016") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "estimated mu in 2016.png")
ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle(expression(paste("Estimated ", mu, " in 2016"))) +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat18 <- dat %>% filter(year == 2018)
png(filename = "observation in 2018.png")
ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Observed Tree coverage in 2018") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "prediction in 2018.png")
ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Predicted Tree coverage in 2018") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "estimated mu in 2018.png")
ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle(expression(paste("Estimated ", mu, " in 2018"))) +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat20 <- dat %>% filter(year == 2020)
png("observation in 2020.png")
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Observed Tree coverage in 2020") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "prediction in 2020.png")
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Predicted Tree coverage in 2020") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "estimated mu in 2020.png")
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle(expression(paste("Estimated ", mu, " in 2020"))) +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = "estimated diffusion.png")
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

png(filename = "estimated growth.png")
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

png(filename = "estimated gamma.png")
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

png(filename = "probability0 in 2020.png")
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

png(filename = "probability0 in 2000.png")
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

png(filename = "probability0 in 2010.png")
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

png(filename = "aggregated stream.png")
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

# Histogram check
histdat.train <- dat %>% filter(year >= 1986 & year <= 1985+TRAIN & edge == 0 & Stream == 0)
histdat.test <- dat %>% filter(year > 1985+TRAIN & edge == 0 & Stream == 0)

png(filename = "histogram comparison training dataset obs vs expectation.png")
par(mfrow=c(1,2))
hist(histdat.train[ ,"trees"], xlab = " ", main = "Distribution of Observed Tree Coverage \n(training)")
hist(histdat.train[ ,"pred.expectation"], xlab = " ", main = "Distribution of Expected Tree Coverage \n(training)")
dev.off()

png(filename = "histogram training dataset obs.png")
par(mfrow=c(1,1))
hist(histdat.train[ ,"trees"], xlab = " ", main = "Distribution of Observed Tree Coverage \n(training)")
dev.off()

png(filename = "histogram training dataset expectation.png")
par(mfrow=c(1,1))
hist(histdat.train[ ,"pred.expectation"], xlab = " ", main = "Distribution of Expected Tree Coverage \n(training)")
dev.off()

png(filename = "histogram comparison test dataset obs vs prediction.png")
par(mfrow=c(1,2))
hist(histdat.test[ ,"trees"], xlab = " ", main = "Distribution of Observed Tree Coverage \n(test)")
hist(histdat.test[ ,"pred.expectation"], xlab = " ", main = "Distribution of Predicted Tree Coverage \n(test)")
dev.off()

png(filename = "histogram test dataset obs.png")
par(mfrow=c(1,1))
hist(histdat.test[ ,"trees"], xlab = " ", main = "Distribution of Observed Tree Coverage \n(test)")
dev.off()

png(filename = "histogram test dataset prediction.png")
par(mfrow=c(1,1))
hist(histdat.test[ ,"pred.expectation"], xlab = " ", main = "Distribution of Predicted Tree Coverage \n(test)")
dev.off()

# Confidence interval for parameters
sink(file = file.path(result_dir, "credible interval.txt"))
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

# MSE calculation
dat.without.edge.train <- dat %>% filter(edge == 0 & Stream == 0 & type == "train")
dat.without.edge.test <- dat %>% filter(edge == 0 & Stream == 0 & type == "test")
MSE.train <- mean((dat.without.edge.train$trees-dat.without.edge.train$pred.expectation)^2)
MSE.test <- mean((dat.without.edge.test$trees-dat.without.edge.test$pred.expectation)^2, na.rm = TRUE)

sink(file = file.path(result_dir, "MSE.txt"))
cat("MSE training: ")
MSE.train
cat("MSE test: ")
MSE.test
sink(file = NULL)

# MSE in terms of homogenization
dat.without.edge.homog.train <- dat %>% filter(type == "train" & edge.homog == 0 & stream.homog == 0)
dat.without.edge.homog.test <- dat %>% filter(type == "test" & edge.homog == 0 & stream.homog == 0)
MSE.train.PDE.homog <- mean((dat.without.edge.homog.train$trees-dat.without.edge.homog.train$pred.expectation)^2, na.rm = TRUE)
MSE.test.PDE.homog <- mean((dat.without.edge.homog.test$trees-dat.without.edge.homog.test$pred.expectation)^2, na.rm = TRUE)

sink(file = file.path(result_dir, "MSE_homog.txt"))
cat("training: year = 1985-2015; test: time = 2016-2020.\n")
cat("MSE PDE training data: ")
MSE.train.PDE.homog
cat("MSE PDE test data: ")
MSE.test.PDE.homog
sink(file = NULL)