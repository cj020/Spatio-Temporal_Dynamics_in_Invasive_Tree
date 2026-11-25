# 5/14/25 gam
# 7/16/25 gam for simulated data 60 by 60
# 7/30/25 real covariates
# 9/2/25 write result as a csv file
# 10/2/25 real data

# Set the path to the local project directory
working_dir <- "C:/Users/YourName/Documents/Spatio-Temporal_Dynamics_in_Invasive_Tree"

# Set working directory
setwd(working_dir)

# Get the working directory
wd <- getwd()

data_dir <- file.path(wd, "real_data", "data", "N300")
output_dir <- file.path(wd, "real_data", "RT")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

required_packages <- c(
  "tidyverse",
  "mgcv",
  "ggplot2",
  "caret",
  "rpart",
  "rpart.plot"
)
missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(tidyverse)
library(mgcv)
library(ggplot2)
library(caret)
library(rpart)
library(rpart.plot)

dat <- read.csv(file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized.csv"))

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
# size in binomial
SIZE <- 100

dat <- dat %>% mutate(edge = ifelse(s1 == 1 |
                                      s1 == N |
                                      s2 == 1|
                                      s2 == N, 1, 0))


# plot
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

# model fitting
rtree.fit <- rpart(trees ~ distance.str + DEM_30m + time + s1 + s2, 
                   data = train,
                   method = "anova", #for regression tree
                   control = rpart.control(minsplit = 30,cp = 0.001))
printcp(rtree.fit)
rsq.rpart(rtree.fit)

png(filename = file.path(output_dir, "Tree_Size.png"))
plotcp(rtree.fit)
dev.off()

summary(rtree.fit)
# plot the tree
png(filename = file.path(output_dir, "Original_Regression_Tree.png"))
plot(rtree.fit, uniform=TRUE, 
     main="Regression Tree for Tree Coverage")
text(rtree.fit, use.n=TRUE, all=TRUE, cex=.3)
dev.off()

png(filename = file.path(output_dir, "Orignial_Regression_Tree_a.png"))
prp(rtree.fit, main = "Original Regression Tree for Tree Coverage")
dev.off()

# prune the tree based on minimim xerror
pruned.rtree.fit <- prune(rtree.fit, cp= rtree.fit$cptable[which.min(rtree.fit$cptable[,"xerror"]),"CP"])
# plot the pruned tree using prp() in the rpart.plot package
png(filename = file.path(output_dir, "Pruned_Regression_Tree.png"))
prp(pruned.rtree.fit, main = "Pruned Regression Tree for Tree Coverage")
dev.off()

newdata <- data.frame(distance.str = dat$distance.str, DEM_30m = dat$DEM_30m, 
                      time = dat$time, s1 = dat$s1, s2 = dat$s2)
prediction.original.rt <- predict(rtree.fit, newdata = newdata)
prediction.pruned.rt <- predict(pruned.rtree.fit, newdata = newdata)

dat <- cbind(dat, prediction.original.rt, prediction.pruned.rt)

write.csv(dat, file.path(output_dir, "prediction_with_RT_real_data_300_by_300.csv"))

pic_limit <- c(0,1)
dat85 <- dat %>% filter(year == 1985)
png(filename=file.path(output_dir, "observation_in_1985.png"))
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

png(filename=file.path(output_dir, "estimation_with_original_regression_tree_in_1985.png"))
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

png(filename=file.path(output_dir, "estimation_with_pruned_regression_tree_in_1985.png"))
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
png(filename=file.path(output_dir, "observation_in_1995.png"))
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

png(filename=file.path(output_dir, "estimation_with_original_regression_tree_in_1995.png"))
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

png(filename=file.path(output_dir, "estimation_with_pruned_regression_tree_in_1995.png"))
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
png(filename=file.path(output_dir, "observation_in_2005.png"))
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

png(filename=file.path(output_dir, "estimation_with_original_regression_tree_in_2005.png"))
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

png(filename=file.path(output_dir, "estimation_with_pruned_regression_tree_in_2005.png"))
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
png(filename=file.path(output_dir, "observation_in_2015.png"))
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

png(filename=file.path(output_dir, "estimation_with_original_regression_tree_in_2015.png"))
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

png(filename=file.path(output_dir, "estimation_with_pruned_regression_tree_in_2015.png"))
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
png(filename=file.path(output_dir, "observation_in_2018.png"))
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

png(filename=file.path(output_dir, "estimation_with_original_regression_tree_in_2018.png"))
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

png(filename=file.path(output_dir, "estimation_with_pruned_regression_tree_in_2018.png"))
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
png(filename=file.path(output_dir, "observation_in_2020.png"))
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

png(filename=file.path(output_dir, "prediction_with_original_regression_tree_in_2020.png"))
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

png(filename=file.path(output_dir, "prediction_with_pruned_regression_tree_in_2020.png"))
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

# png(filename = file.path(output_dir, "tree_02.png"))
# ggarrange(obs0p, pred0p,
#           obs2p, pred2p, 
#           nrow = 2, ncol = 2, byrow = T)
# dev.off()
# 
# png(filename = file.path(output_dir, "tree_46.png"))
# ggarrange(obs4p, pred4p,
#           obs6p, pred6p, 
#           nrow = 2, ncol = 2, byrow = T)
# dev.off()
# 
# png(filename = file.path(output_dir, "tree_810.png"))
# ggarrange(obs8p, pred8p,
#           obs10p, pred10p, 
#           nrow = 2, ncol = 2, byrow = T)
# dev.off()

# histrogram
png(filename=file.path(output_dir, "histogram_comparison_training_dataset_obs_vs_original_rt_est.png"))
par(mfrow=c(1,2))
histdat.train <- dat %>% filter(type == "train" & edge == 0 & Stream == 0)
hist(histdat.train[ ,"trees"], main = "Histogram of observations \n(training)")
hist(histdat.train[ ,"prediction.original.rt"], main = "Histogram of estimation with orignial rt \n(training)")
dev.off()

png(filename=file.path(output_dir, "histogram_training_dataset_obs.png"))
par(mfrow=c(1,1))
hist(histdat.train[ ,"trees"], main = "Histogram of observations \n(training)")
dev.off()

png(filename=file.path(output_dir, "histogram_training_dataset_original_rt_est.png"))
par(mfrow=c(1,1))
hist(histdat.train[ ,"prediction.original.rt"], main = "Histogram of estimation with orignial rt \n(training)")
dev.off()

png(filename=file.path(output_dir, "histogram_comparison_training_dataset_obs_vs_pruned_rt_est.png"))
par(mfrow=c(1,2))
hist(histdat.train[ ,"trees"], main = "Histogram of observations \n(training)")
hist(histdat.train[ ,"prediction.pruned.rt"], main = "Histogram of estimation with pruned rt \n(training)")
dev.off()

png(filename=file.path(output_dir, "histogram_training_dataset_pruned_rt_est.png"))
par(mfrow=c(1,1))
hist(histdat.train[ ,"prediction.pruned.rt"], main = "Histogram of estimation with pruned rt \n(training)")
dev.off()

png(filename=file.path(output_dir, "histogram_comparison_test_dataset_obs_vs_original_rt_est.png"))
par(mfrow=c(1,2))
histdat.test <- dat %>% filter(type == "test" & edge == 0 & Stream == 0)
hist(histdat.test[ ,"trees"], main = "Histogram of observations \n(test)")
hist(histdat.test[ ,"prediction.original.rt"], main = "Histogram of prediction with original rt \n(test)")
dev.off()

png(filename=file.path(output_dir, "histogram_test_dataset_obs.png"))
par(mfrow=c(1,1))
hist(histdat.test[ ,"trees"], main = "Histogram of observations \n(test)")
dev.off()

png(filename=file.path(output_dir, "histogram_test_dataset_original_rt_est.png"))
par(mfrow=c(1,1))
hist(histdat.test[ ,"prediction.original.rt"], main = "Histogram of prediction with original rt \n(test)")
dev.off()

png(filename=file.path(output_dir, "histogram_comparison_test_dataset_obs_vs_pruned_rt_est.png"))
par(mfrow=c(1,2))
hist(histdat.test[ ,"trees"], main = "Histogram of observations \n(test)")
hist(histdat.test[ ,"prediction.pruned.rt"], main = "Histogram of prediction with pruned rt \n(test)")
dev.off()

png(filename=file.path(output_dir, "histogram_test_dataset_pruned_rt_est.png"))
par(mfrow=c(1,1))
hist(histdat.test[ ,"prediction.pruned.rt"], main = "Histogram of prediction with pruned rt \n(test)")
dev.off()

# par(mfrow=c(1,2))
# histdat.test <- dat %>% filter(time > 0+TRAIN & edge.homog == 0 & stream.homog == 0)
# hist(histdat.test[ ,"trees"], main = "Histogram of observations \n(test)")
# hist(histdat.test[ ,"prediction"], main = "Histogram of forecast \n(test)")

dat.without.edge.train <- dat %>% filter(type == "train" & edge == 0 & Stream == 0)
dat.without.edge.test <- dat %>% filter(type == "test" & edge == 0 & Stream == 0)
MSE.train.original.rt <- mean((dat.without.edge.train$trees-dat.without.edge.train$prediction.original.rt)^2, na.rm = TRUE)
MSE.train.pruned.rt <- mean((dat.without.edge.train$trees-dat.without.edge.train$prediction.pruned.rt)^2, na.rm = TRUE)
MSE.test.original.rt <- mean((dat.without.edge.test$trees-dat.without.edge.test$prediction.original.rt)^2, na.rm = TRUE)
MSE.test.pruned.rt <- mean((dat.without.edge.test$trees-dat.without.edge.test$prediction.pruned.rt)^2, na.rm = TRUE)

sink(file = file.path(output_dir, "MSE_RT.txt"))
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
