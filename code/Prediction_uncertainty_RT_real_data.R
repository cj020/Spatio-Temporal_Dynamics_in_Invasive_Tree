# 5/14/25 gam
# 7/16/25 gam for simulated data 60 by 60
# 9/1/25 write result as csv file
# 11/16/25 prediction uncertainty

# Set the path to the local project directory
working_dir <- "C:/Users/YourName/Documents/Spatio-Temporal_Dynamics_in_Invasive_Tree"

# Set working directory
setwd(working_dir)

# Get the working directory
wd <- getwd()

data_dir <- file.path(wd, "real_data", "data", "N300")
uncertainty_dir <- file.path(wd, "real_data", "prediction_summary", "MSE_and_prediction_uncertainty")

if (!dir.exists(uncertainty_dir)) {
  dir.create(uncertainty_dir, recursive = TRUE)
}

required_packages <- c(
  "tidyverse",
  "mgcv",
  "ggplot2",
  "rpart"
)
missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(tidyverse)
library(mgcv)
library(ggplot2)
# library(caret)
library(rpart)
# library(rpart.plot)

# /home/yawenguan/jiaqichen/Region88/paper/real data/prediction uncertainty")
dat <- read.csv(file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized.csv"))


# bootstrap
set.seed(12345)

B <- 200
pred_mat <- matrix(NA, nrow = nrow(dat), ncol = B)

# parameters setting
# number of grids N by N
N <- 300
# number of years
TIME <- 35
# Length of grid
LoG <- 30
# number of years as training data
TRAIN <- 31
# number of covariates
K <- 2

time <- dat$year-1985
dat <- cbind(dat, time)
train <- dat %>% filter(year <= 1984+TRAIN)
dat <- dat %>% mutate(type = ifelse(year <= 1984+TRAIN, "train", "test"))

for (b in 1:B) {
  idx <- sample(nrow(train), replace = TRUE)
  
  rtree.fit <- rpart(trees ~ distance.str + DEM_30m + time + s1 + s2, 
                     data = train[idx, ],
                     method = "anova", #for regression tree
                     control = rpart.control(minsplit = 30,cp = 0.001))
  
  # prune the tree based on minimim xerror
  pruned.rtree.fit <- prune(rtree.fit, cp= rtree.fit$cptable[which.min(rtree.fit$cptable[,"xerror"]),"CP"])
  
  newdata <- data.frame(distance.str = dat$distance.str, DEM_30m = dat$DEM_30m,  
                        time = dat$time, s1 = dat$s1, s2 = dat$s2)
  # prediction.original.rt <- predict(rtree.fit, newdata = newdata)
  prediction.pruned.rt <- predict(pruned.rtree.fit, newdata = newdata)
  prediction_3 <- round(prediction.pruned.rt, 3)
  
  pred_mat <- cbind(pred_mat, prediction_3)
}

lower <- apply(pred_mat, 1, quantile, 0.025, na.rm = TRUE)
upper <- apply(pred_mat, 1, quantile, 0.975, na.rm = TRUE)

sum_withinCI <- rep(0, nrow(dat))

for (l in 1:nrow(dat)) {
  if(dat$trees[l] >= lower[l] & dat$trees[l] <= upper[l])
  sum_withinCI[l] <- 1
}

# write.csv(dat, "prediction_with_RT_result_real_data_300_by_300.csv")

sink(file = file.path(uncertainty_dir, "prediction_coverage_RT_method_B200.txt"))
cat("prediction coverage in training data: ")
mean(sum_withinCI[1:(N*N*(TRAIN))])
cat("prediction coverage in test data: ")
mean(sum_withinCI[(N*N*(TRAIN)+1):nrow(dat)])
sink(file = NULL)

