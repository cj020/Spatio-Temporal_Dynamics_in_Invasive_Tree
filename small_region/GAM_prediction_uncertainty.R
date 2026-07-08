# 5/14/25 gam
# 7/16/25 gam for simulated data 60 by 60
# 9/1/25 write result as csv file
# 11/11/25 prediction uncertainty
# 7/6/26 try "gam <- gam(cbind(tree100, size) ~ s(distance.str)+ s(DEM_30m) + s(time) + s(s1, s2), method = "REML", family = binomial(link = "logit"), data = train)"
# and "gam <- gam(cbind(tree100, size) ~ distance.str+ DEM_30m + time + s(s1, s2), method = "REML", family = binomial(link = "logit"), data = train)"
# 7/7/2026 model error fixed

library(tidyverse)
library(mgcv)
library(ggplot2)

data_dir <- "/home/appliedmath/jiaqichen/Region88/paper/real data/data/N300/"
output_dir <- "/home/appliedmath/jiaqichen/Region88/paper/real data/GAM/20260707 GAM for real data/model2"
result_dir <- "/home/appliedmath/jiaqichen/Region88/paper/real data/GAM/20260707 GAM for real data/model2"


# setwd("/home/yawenguan/jiaqichen/Region88/paper/real data/GAM/20251002 GAM for real data")
dat <- read.csv(file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized.csv"))

# parameters setting
# number of grids N by N
N <- 300
# number of years
TIME <- 35
# Length of grid
LoG <- 30
# number of years as training data
TRAIN <- 30
# number of covariates
K <- 2
# size in binomial
SIZE <- 100

tree100 <- 100 * dat$trees
size <- rep(SIZE, N*N*(TIME+1))
time <- dat$year-1985
dat <- cbind(dat, time, tree100, size)
train <- dat %>% filter(time <= TRAIN)
dat <- dat %>% mutate(type=ifelse(time <= 0+TRAIN, "train", "test"))

# model fitting
# Model 0
# gam <- gam(cbind(tree100, size-tree100) ~ s(distance.str) + DEM_30m + time + s(s1) + s(s2), method = "REML",
#            family = binomial(link = "logit"), data = train)

# Model 1
# gam <- gam(cbind(tree100, size-tree100) ~ s(distance.str) + s(DEM_30m) + s(time) + s(s1, s2), method = "REML",
#            family = binomial(link = "logit"), data = train)

# Model 2
gam <- gam(cbind(tree100, size-tree100) ~ distance.str + DEM_30m + time + s(s1, s2), method = "REML",
           family = binomial(link = "logit"), data = train)

newdata <- data.frame(distance.str = dat$distance.str, DEM_30m = dat$DEM_30m, 
                      time = dat$time, s1 = dat$s1, s2=dat$s2, size = dat$size)
pred <- predict(gam, newdata = newdata, type = "link", se.fit = TRUE)

# 95% CI
crit <- qnorm(0.975)

prediction <- plogis(pred$fit)

lower_CI <- plogis(pred$fit - crit * pred$se.fit)
upper_CI <- plogis(pred$fit + crit * pred$se.fit)
lower_CI_3 <- round(lower_CI, 3)
upper_CI_3 <- round(upper_CI, 3)

dat <- cbind(dat, prediction, lower_CI_3, upper_CI_3)

# write.csv(dat, "prediction_with_GAM_result_simulated_data_60_by_60.csv")

sum_withinCI <- rep(0, nrow(dat))
for (j in 1:nrow(dat)) {
  if(dat[j,"trees"] <= dat[j,"upper_CI_3"] & dat[j,"trees"] >= dat[j,"lower_CI_3"])
    sum_withinCI[j] <- 1
}

sink(file = file.path(result_dir, "prediction coverage GAM method.txt"))
cat("prediction coverage in training data: ")
mean(sum_withinCI[1:(N*N*(TRAIN+1))])
cat("prediction coverage in test data: ")
mean(sum_withinCI[(N*N*(TRAIN+1)+1):nrow(dat)])
sink(file = NULL)
