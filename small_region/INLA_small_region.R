# INLA in small region with spatial dependent, dist+elevation+time+s1+s2

library(spdep)
library(INLA)
library(tidyverse)
library(ggplot2)

setwd("/small_region")

data_dir <- "/small_region/data"
output_dir <- "/small_region/output/INLA"
result_dir <- "/small_region/result/INLA"

dat <- read.csv(file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized.csv"))

# parameters setting
# N (number of grids = N * N)
N <- 300
# Number of years
TIME <- 35
# Length of grid
LoG <- 30
# Number of years as training data
TRAIN <- 31
# Number of covariates
K <- 2
# Size in binomial
SIZE <- 100

# Set up edge label
dat <- dat %>% mutate(edge = ifelse(s1 == 1 |
                                      s1 == N |
                                      s2 == 1|
                                      s2 == N, 1, 0))

# Plot check - tree cover in 1985
dat85 <- dat %>% filter(year == 1985)
ggplot(data = dat85, aes(s1, s2, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
  labs(title = "observation at 1985") +
  coord_equal() 

# Setup spatial index
spatial_idx <- rep((1:(N*N)), TIME+1)

tree100 <- as.integer(100 * dat$trees)
response100 <- tree100
response100[((N*N*TRAIN)+1):(N*N*(TIME+1))] <- NA
size <- rep(SIZE, N*N*(TIME+1))
time <- dat$year-1985
dat <- cbind(dat, spatial_idx, tree100, response100, size, time)
train <- dat %>% filter(time <= 0+TRAIN)
dat <- dat %>% mutate(type=ifelse(time <= 0+TRAIN, "train", "test"))

nb <- spdep::cell2nb(N, N, type = "rook")
class(nb)
spdep::nb2INLA(file = "grid_adj.graph", nb = nb)

dat$space_id <- (dat$s1 - 1) * N + dat$s2       # 1..N^2
dat$t_id <- as.integer(factor(dat$time))        # 1..(TIME+1)

# Model fitting - INLA with binomial (sp) 
formula_sp <- response100 ~ distance.str + DEM_30m + time + s1 + s2 +
  f(space_id, model = "besag", graph="grid_adj.graph")
INLA_fit_sp <- inla(formula_sp, data = dat,
                    family = "binomial",
                    Ntrials = size,
                    control.fixed = list(
                      mean = 0,
                      prec = 1 / 10^2,          # precision = 1/variance
                      mean.intercept = 0,
                      prec.intercept = 1 / 10^2),
                    control.predictor = list(compute = TRUE, link = 1),
                    control.compute = list(dic = TRUE, return.marginals.predictor = TRUE))

sink(file = file.path(result_dir, "summary_sp.txt"))
INLA_fit_sp$summary.fixed
sink(NULL)

prediction_sp <- INLA_fit_sp$summary.fitted.values$mean
lower_CI_sp <- INLA_fit_sp$summary.fitted.values$`0.025quant`
upper_CI_sp <- INLA_fit_sp$summary.fitted.values$`0.975quant`
lower_CI_3_sp <- round(lower_CI_sp, 3)
upper_CI_3_sp <- round(upper_CI_sp, 3)

dat <- cbind(dat, prediction_sp, lower_CI_3_sp, upper_CI_3_sp)

write.csv(dat, file.path(result_dir, "prediction_with_INLA_result_real_data_300_by_300.csv"))

# Graph check
pic_limit <- c(0,1)
dat90 <- dat %>% filter(year == 1990)
png(filename = file.path(output_dir, "observation in 1990.png"))
ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Observed Tree Coverage at year 1990") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "estimation in 1990 sp.png"))
ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = prediction_sp)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage at year 1990 (sp)") +
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
png(filename = file.path(output_dir, "observation in 2000.png"))
ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Observed Tree Coverage at year 2000") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "estimation in 2000 sp.png"))
ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = prediction_sp)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage at year 2000 (sp)") +
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
png(filename = file.path(output_dir, "observation in 2010.png"))
ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Observed Tree Coverage at year 2010") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "estimation in 2010 sp.png"))
ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = prediction_sp)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage at year 2010 (sp)") +
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
png(filename = file.path(output_dir, "observation in 2020.png"))
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Observed Tree Coverage at year 2020") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "estimation in 2020 sp.png"))
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = prediction_sp)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage at year 2020 (sp)") +
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
png(filename = file.path(output_dir, "histogram comparison training dataset obs vs est.png"),
    height = 1000, width = 1600, res = 100)
par(mfrow=c(1,2),
    cex.lab = 2,
    cex.axis = 1.8, 
    mar = c(5, 5, 2, 2),
    oma = c(1, 1, 1, 1))
histdat.train <- dat %>% filter(type == "train" & edge == 0 & Stream == 0)
hist(histdat.train[ ,"trees"], main = "Histogram of Observations \n(training)")
hist(histdat.train[ ,"prediction_sp"], main = "Histogram of Estimation (sp)\n(training)")
dev.off()

png(filename = file.path(output_dir, "histogram comparison test dataset obs vs est.png"),
    height = 1000, width = 1600, res = 100)
par(mfrow=c(1,2),
    cex.lab = 2,
    cex.axis = 1.8, 
    mar = c(5, 5, 2, 2),
    oma = c(1, 1, 1, 1))
histdat.test <- dat %>% filter(type == "test" & edge == 0 & Stream == 0)
hist(histdat.test[ ,"trees"], main = "Histogram of Observations \n(test)")
hist(histdat.test[ ,"prediction_sp"], main = "Histogram of Prediction (sp)\n(test)")
dev.off()

dat.without.edge.train <- dat %>% filter(type == "train" & edge == 0 & Stream == 0)
dat.without.edge.test <- dat %>% filter(type == "test" & edge == 0 & Stream == 0)
MSE.train.sp <- mean((dat.without.edge.train$trees-dat.without.edge.train$prediction_sp)^2, na.rm = TRUE)
MSE.test.sp <- mean((dat.without.edge.test$trees-dat.without.edge.test$prediction_sp)^2, na.rm = TRUE)

# MSE
sink(file = file.path(result_dir, "MSE.txt"))
cat("training: time 0-8; test: time 9-10.\n")
cat("MSE of training dataset (sp): ")
MSE.train.sp
cat("MSE of test dataset (sp): ")
MSE.test.sp
sink(file = NULL)

sum_withinCI_sp <- rep(0, nrow(dat))
for (j in 1:nrow(dat)) {
  if(dat[j,"trees"] <= dat[j,"upper_CI_3_sp"] & dat[j,"trees"] >= dat[j,"lower_CI_3_sp"])
    sum_withinCI_sp[j] <- 1
}

# Prediction coverage
sink(file = file.path(result_dir, "prediction coverage INLA method.txt"))
cat("prediction coverage in training data (sp): ")
mean(sum_withinCI_sp[1:(N*N*(TRAIN+1))])
cat("prediction coverage in test data (sp): ")
mean(sum_withinCI_sp[(N*N*(TRAIN+1)+1):nrow(dat)])
sink(file = NULL)

sink(file = file.path(result_dir, "/prediction coverage INLA method.txt"))
cat("prediction coverage in training data (sp): ")
mean(sum_withinCI_sp[1:(N*N*(TRAIN+1))])
cat("prediction coverage in test data (sp): ")
mean(sum_withinCI_sp[(N*N*(TRAIN+1)+1):nrow(dat)])
sink(file = NULL)