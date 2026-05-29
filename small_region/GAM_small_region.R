# GAM in small region

library(tidyverse)
library(mgcv)
library(ggplot2)

setwd("/small_region")

data_dir <- "/small_region/data"
output_dir <- "/small_region/output/GAM"
result_dir <- "/small_region/result/GAM"

dat <- read.csv(file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized.csv"))

# Parameters setting
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


# Graph check - tree cover in 1985
dat85 <- dat %>% filter(year == 1985)
ggplot(data = dat85, aes(s1, s2, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
  labs(title = "observation at 1985") +
  coord_equal() 

# Convert observation to integers between 0 to 100
tree100 <- 100 * dat$trees
size <- rep(SIZE, N*N*(TIME+1))
time <- dat$year-1985
dat <- cbind(dat, time, tree100, size)
train <- dat %>% filter(year <= 1984+TRAIN)
dat <- dat %>% mutate(type = ifelse(year <= 1984+TRAIN, "train", "test"))

# Model fitting
gam <- gam(cbind(tree100, size) ~ s(distance.str) + DEM_30m + time + s(s1) + s(s2), method = "REML", 
           family = binomial(link = "logit"), data = train)
sink(file = file.path(result_dir, "summary.txt"))
summary(gam)
sink(file = NULL)

# Estimate and prediction of tree cover
newdata <- data.frame(distance.str = dat$distance.str, DEM_30m = dat$DEM_30m, 
                      time = dat$time, s1 = dat$s1, s2=dat$s2, size = dat$size)
prediction <- predict(gam, newdata = newdata, type = "response")

dat <- cbind(dat, prediction)

# Save data
write.csv(dat, file.path(result_dir, "prediction_with_GAM_real_data_300_by_300.csv"))

# Graph Check
pic_limit <- c(0,1)
dat85 <- dat %>% filter(year == 1985)
png(filename = file.path(output_dir, "observation in 1985.png"))
ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Observed Tree Coverage in 1985") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "estimation in 1985.png"))
ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = prediction)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimateed Tree Coverage in 1985") +
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
png(filename = file.path(output_dir, "observation in 1995.png"))
ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Observed Tree Coverage in 1995") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "estimation in 1995.png"))
ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = prediction)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage in 1995") +
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
png(filename = file.path(output_dir, "observation in 2005.png"))
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Observed Tree Coverage in 2005") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal()  
dev.off()

png(filename = file.path(output_dir, "estimation in 2005.png"))
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = prediction)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage in 2005") +
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
png(filename = file.path(output_dir, "observation in 2015.png"))
ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Observed Tree Coverage in 2015") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "estimation in 2015.png"))
ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = prediction)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage in 2015") +
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
png(filename = file.path(output_dir, "observation in 2018.png"))
ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Observed Tree Coverage 2018") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "prediction in 2018.png"))
ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = prediction)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage in 2018") +
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
  labs(title = "Observed Tree Coverage in 2020") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "prediction in 2020.png"))
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = prediction)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage in 2020") +
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
png(filename = file.path(output_dir, "histogram comparison training dataset obs vs est.png"))
par(mfrow=c(1,2))
histdat.train <- dat %>% filter(type == "train" & edge == 0 & Stream == 0)
hist(histdat.train[ ,"trees"], main = "Histogram of Observations \n(training)")
hist(histdat.train[ ,"prediction"], main = "Histogram of Estimation \n(training)")
dev.off()

png(filename = file.path(output_dir, "histogram comparison test dataset obs vs est.png"))
par(mfrow=c(1,2))
histdat.test <- dat %>% filter(type == "test" & edge == 0 & Stream == 0)
hist(histdat.test[ ,"trees"], main = "Histogram of Observations \n(test)")
hist(histdat.test[ ,"prediction"], main = "Histogram of Prediction \n(test)")
dev.off()

png(filename = file.path(output_dir, "histogram training dataset obs.png"))
par(mfrow=c(1,1))
hist(histdat.train[ ,"trees"], main = "Histogram of Observations \n(training)")
dev.off()

png(filename = file.path(output_dir, "histogram training dataset est.png"))
par(mfrow=c(1,1))
hist(histdat.train[ ,"prediction"], main = "Histogram of Estimation \n(training)")
dev.off()

png(filename = file.path(output_dir, "histogram test dataset obs.png"))
par(mfrow=c(1,1))
hist(histdat.test[ ,"trees"], main = "Histogram of Observations \n(test)")
dev.off()

png(filename = file.path(output_dir, "histogram test dataset pred.png"))
par(mfrow=c(1,1))
hist(histdat.test[ ,"prediction"], main = "Histogram of Prediction \n(test)")
dev.off()

dat.without.edge.train <- dat %>% filter(type == "train" & edge == 0 & Stream == 0)
dat.without.edge.test <- dat %>% filter(type == "test" & edge == 0 & Stream == 0)
MSE.train <- mean((dat.without.edge.train$trees-dat.without.edge.train$prediction)^2, na.rm = TRUE)
MSE.test <- mean((dat.without.edge.test$trees-dat.without.edge.test$prediction)^2, na.rm = TRUE)

sink(file = file.path(result_dir, "MSE.txt"))
cat("training: year 1985-2015; test: year 2016-2020.\n")
cat("MSE of training dataset: ")
MSE.train
cat("MSE of test dataset: ")
MSE.test
sink(file = NULL)
