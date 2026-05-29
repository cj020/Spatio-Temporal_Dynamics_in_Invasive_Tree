# RT GAM in small region

library(tidyverse)
library(mgcv)
library(ggplot2)
library(caret)
library(rpart)
library(rpart.plot)

setwd("/small_region")

data_dir <- "/small_region/data"
output_dir <- "/small_region/output/RT"
result_dir <- "/small_region/result/RT"

dat <- read.csv(file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized.csv"))

# Parameters setting
# N (number of grids = N * N)
N <- 300
# number of years
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

time <- dat$year-1985
dat <- cbind(dat, time)
train <- dat %>% filter(year <= 1984+TRAIN)
dat <- dat %>% mutate(type = ifelse(year <= 1984+TRAIN, "train", "test"))

# Model fitting
rtree.fit <- rpart(trees ~ distance.str + DEM_30m + time + s1 + s2, 
                   data = train,
                   method = "anova", #for regression tree
                   control = rpart.control(minsplit = 30,cp = 0.001))
printcp(rtree.fit)
rsq.rpart(rtree.fit)

png(filename = file.path(result_dir, "Tree Size.png"))
plotcp(rtree.fit)
dev.off()

summary(rtree.fit)

# Plot the tree
png(filename = file.path(result_dir, "Original Regression Tree.png"))
plot(rtree.fit, uniform=TRUE, 
     main="Regression Tree for Tree Coverage")
text(rtree.fit, use.n=TRUE, all=TRUE, cex=.3)
dev.off()

png(filename = file.path(result_dir, "Orignial Regression Tree a.png"))
prp(rtree.fit, main = "Original Regression Tree for Tree Coverage")
dev.off()

# Prune the tree based on minimim xerror
pruned.rtree.fit <- prune(rtree.fit, cp= rtree.fit$cptable[which.min(rtree.fit$cptable[,"xerror"]),"CP"])

# Plot the pruned tree using prp() in the rpart.plot package
png(filename = file.path(result_dir, "Pruned Regression Tree.png"))
prp(pruned.rtree.fit, main = "Pruned Regression Tree for Tree Coverage")
dev.off()

# Estimate and predict tree cover
newdata <- data.frame(distance.str = dat$distance.str, DEM_30m = dat$DEM_30m, 
                      time = dat$time, s1 = dat$s1, s2 = dat$s2)
prediction.original.rt <- predict(rtree.fit, newdata = newdata)
prediction.pruned.rt <- predict(pruned.rtree.fit, newdata = newdata)

dat <- cbind(dat, prediction.original.rt, prediction.pruned.rt)

write.csv(dat, file.path(result_dir, "prediction_with_RT_real_data_300_by_300.csv"))

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

png(filename = file.path(output_dir, "estimation with original regression tree in 1985.png"))
ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = prediction.original.rt)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage with Original R.T. in 1985") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "estimation with pruned regression tree in 1985.png"))
ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = prediction.pruned.rt)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage with Pruned R.T. in 1985") +
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

png(filename = file.path(output_dir, "estimation with original regression tree in 1995.png"))
ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = prediction.original.rt)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage with Original R.T. in 1995") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "estimation with pruned regression tree in 1995.png"))
ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = prediction.pruned.rt)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage with Pruned R.T. in 1995") +
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

png(filename = file.path(output_dir, "estimation with original regression tree in 2005.png"))
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = prediction.original.rt)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage with Original R.T. in 2005") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "estimation with pruned regression tree in 2005.png"))
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = prediction.pruned.rt)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage with Pruned R.T. in 2005") +
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

png(filename = file.path(output_dir, "estimation with original regression tree in 2015.png"))
ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = prediction.original.rt)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage with Original R.T. in 2015") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "estimation with pruned regression tree in 2015.png"))
ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = prediction.pruned.rt)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage with Pruned R.T. in 2015") +
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
  labs(title = "Observed Tree Coverage in 2018") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "estimation with original regression tree in 2018.png"))
ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = prediction.original.rt)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage with Original R.T.in 2018") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "estimation with pruned regression tree in 2018.png"))
ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = prediction.pruned.rt)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Estimated Tree Coverage with Pruned R.T. in 2018") +
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

png(filename = file.path(output_dir, "prediction with original regression tree in 2020.png"))
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = prediction.original.rt)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Predicted Tree Coverage with Original R.T. in 2020") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename = file.path(output_dir, "prediction with pruned regression tree in 2020.png"))
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = prediction.pruned.rt)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Predicted Tree Coverage with Pruned R.T. in 2020") +
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
png(filename = file.path(output_dir, "histogram comparison training dataset obs vs original rt est.png"))
par(mfrow=c(1,2))
histdat.train <- dat %>% filter(type == "train" & edge == 0 & Stream == 0)
hist(histdat.train[ ,"trees"], main = "Histogram of observations \n(training)")
hist(histdat.train[ ,"prediction.original.rt"], main = "Histogram of estimation with orignial rt \n(training)")
dev.off()

png(filename = file.path(output_dir, "histogram training dataset obs.png"))
par(mfrow=c(1,1))
hist(histdat.train[ ,"trees"], main = "Histogram of observations \n(training)")
dev.off()

png(filename = file.path(output_dir, "histogram training dataset original rt est.png"))
par(mfrow=c(1,1))
hist(histdat.train[ ,"prediction.original.rt"], main = "Histogram of estimation with orignial rt \n(training)")
dev.off()

png(filename = file.path(output_dir, "histogram comparison training dataset obs vs pruned rt est.png"))
par(mfrow=c(1,2))
hist(histdat.train[ ,"trees"], main = "Histogram of observations \n(training)")
hist(histdat.train[ ,"prediction.pruned.rt"], main = "Histogram of estimation with pruned rt \n(training)")
dev.off()

png(filename = file.path(output_dir, "histogram training dataset pruned rt est.png"))
par(mfrow=c(1,1))
hist(histdat.train[ ,"prediction.pruned.rt"], main = "Histogram of estimation with pruned rt \n(training)")
dev.off()

png(filename = file.path(output_dir, "histogram comparison test dataset obs vs original rt est.png"))
par(mfrow=c(1,2))
histdat.test <- dat %>% filter(type == "test" & edge == 0 & Stream == 0)
hist(histdat.test[ ,"trees"], main = "Histogram of observations \n(test)")
hist(histdat.test[ ,"prediction.original.rt"], main = "Histogram of prediction with original rt \n(test)")
dev.off()

png(filename = file.path(output_dir, "histogram test dataset obs.png"))
par(mfrow=c(1,1))
hist(histdat.test[ ,"trees"], main = "Histogram of observations \n(test)")
dev.off()

png(filename = file.path(output_dir, "histogram test dataset original rt est.png"))
par(mfrow=c(1,1))
hist(histdat.test[ ,"prediction.original.rt"], main = "Histogram of prediction with original rt \n(test)")
dev.off()

png(filename = file.path(output_dir, "histogram comparison test dataset obs vs pruned rt est.png"))
par(mfrow=c(1,2))
hist(histdat.test[ ,"trees"], main = "Histogram of observations \n(test)")
hist(histdat.test[ ,"prediction.pruned.rt"], main = "Histogram of prediction with pruned rt \n(test)")
dev.off()

png(filename = file.path(output_dir, "histogram test dataset pruned rt est.png"))
par(mfrow=c(1,1))
hist(histdat.test[ ,"prediction.pruned.rt"], main = "Histogram of prediction with pruned rt \n(test)")
dev.off()

dat.without.edge.train <- dat %>% filter(type == "train" & edge == 0 & Stream == 0)
dat.without.edge.test <- dat %>% filter(type == "test" & edge == 0 & Stream == 0)
MSE.train.original.rt <- mean((dat.without.edge.train$trees-dat.without.edge.train$prediction.original.rt)^2, na.rm = TRUE)
MSE.train.pruned.rt <- mean((dat.without.edge.train$trees-dat.without.edge.train$prediction.pruned.rt)^2, na.rm = TRUE)
MSE.test.original.rt <- mean((dat.without.edge.test$trees-dat.without.edge.test$prediction.original.rt)^2, na.rm = TRUE)
MSE.test.pruned.rt <- mean((dat.without.edge.test$trees-dat.without.edge.test$prediction.pruned.rt)^2, na.rm = TRUE)

sink(file = file.path(result_dir, "MSE.txt"))
cat("training: year = 1985-2015; test: time = 2016-2020.\n")
cat("MSE training data original RT: ")
MSE.train.original.rt
cat("MSE training data pruned RT: ")
MSE.train.pruned.rt
cat("MSE test data original RT: ")
MSE.test.original.rt
cat("MSE test data pruned RT: ")
MSE.test.pruned.rt
sink(file = NULL)