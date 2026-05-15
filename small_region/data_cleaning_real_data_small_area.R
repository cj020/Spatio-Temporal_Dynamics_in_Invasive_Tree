# Data cleaning: smooth, calculate distance to stream of each grid, standardize the covariates

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

# N (number of grids = N * N)
N <- 300
# Length of grid
LoG <- 30
# number of years 1986-2020
TIME <- 35
# size in binomial
SIZE <- 100

setwd("/small_region/data/")
dat.original <- read.csv("/small_region/data/dat_8520_small_area_300_by_300_original.csv")

dat <- dat.original[, c("year", "s1", "s2", "Longitude", "Latitude", "Stream", "DEM_30m", "trees")]
index <- rep(1:(N*N), TIME+1)

dat <- cbind(index, dat)

# Smoothing by averaging consecutive 3 years
trees.smth <- c()
prestyr <- dat[which(dat$year==1985), "trees"]
nextyr <- dat[which(dat$year==1986), "trees"]
smooth.trees.prest <- (prestyr + nextyr)/2
trees.smth <- c(trees.smth, round(smooth.trees.prest, digits = 0))

for (i in 1:(TIME-1)) {
  prest <- 1985+i
  prestyr <- dat[which(dat$year==prest), "trees"]
  prevyr <- dat[which(dat$year==(prest-1)), "trees"]
  nextyr <- dat[which(dat$year==(prest+1)), "trees"]
  
  smooth.trees.prest <- (prevyr + prestyr + nextyr)/3
  trees.smth <- c(trees.smth, round(smooth.trees.prest, digits = 0))
}

prestyr <- dat[which(dat$year==2020), "trees"]
prevyr <- dat[which(dat$year==2019), "trees"]
smooth.trees.prest <- (prevyr + prestyr)/2
trees.smth <- c(trees.smth, round(smooth.trees.prest, digits = 0))

dat$original.trees <- dat$trees
dat$trees <- trees.smth/SIZE

write.csv(dat, "/small_region/data/dat_8520_small_area_300_by_300_smoothed_integered.csv")

# Add a covariate as distance to the nearest stream
dat <- read.csv("/small_region/data/dat_8520_small_area_300_by_300_smoothed_integered.csv")
data <- dat %>% filter(year == 1985)
coord <- cbind(data["s1"], data["s2"])
stream.cell <- data %>% filter(Stream == 1)
coord.stream <- cbind(stream.cell["s1"], stream.cell["s2"])
knnx <- get.knnx(coord.stream, coord, k = 1)
dist.str0 <- LoG*knnx$nn.dist
direc.h0 <- rep(0, dim(data)[1])
direc.v0 <- rep(0, dim(data)[1])
for (i in 1:dim(data)[1]) {
  str <- knnx$nn.index[i]
  direc.h0[i] <- LoG*(coord[i,"s1"]-coord.stream[str,"s1"])
  direc.v0[i] <- LoG*(coord[i,"s2"]-coord.stream[str,"s2"])
}

distance.str <- rep(dist.str0, TIME+1)
direc.to.str.h <- rep(direc.h0, TIME+1)
direc.to.str.v <- rep(direc.v0, TIME+1)

dat <- cbind(dat, distance.str, direc.to.str.h, direc.to.str.v)

# Standardize dist, Dem_30, Stream
dat$original.dist <- dat$distance.str
dat$original.DEM_30m <- dat$DEM_30m

dat_85 <- dat %>% filter(year == 1985)
dist_min <- min(dat_85$distance.str)
dist_max <- max(dat_85$distance.str)
dem30m_min <- min(dat_85$DEM_30m)
dem30m_max <- max(dat_85$DEM_30m)

dist85.std <- (dat_85$distance.str-dist_min)/(dist_max-dist_min)
dem30m.std <- (dat_85$DEM_30m-dem30m_min)/(dem30m_max-dem30m_min)

dat$distance.str <- rep(dist85.std, TIME+1)
dat$DEM_30m <- rep(dem30m.std, TIME+1)
dat$Stream <- rep(dat_85$Stream, TIME+1)

write.csv(dat, "/small_region/data/dat_8520_small_area_300_by_300_smoothed_integered_standardized.csv")

# Graph
dat85 <- dat %>% filter(year == 1985)
png(filename = "/small_region/data/plots/observation in 1985 smoothed 300 by 300.png")
ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 1985") +
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
png(filename = "/small_region/data/plots/observation in 1990 smoothed 300 by 300.png")
ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 1990") +
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
png(filename = "/small_region/data/plots/observation in 1995 smoothed 300 by 300.png")
ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 1995") +
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
png(filename = "/small_region/data/plots/observation in 2000 smoothed 300 by 300.png")
ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 2000") +
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
png(filename = "/small_region/data/plots/observation in 2005 smoothed 300 by 300.png")
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 2005") +
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
png(filename = "/small_region/data/plots/observation in 2010 smoothed 300 by 300.png")
ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 2010") +
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
png(filename = "/small_region/data/plots/observation in 2015 smoothed 300 by 300.png")
ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 2015") +
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
png(filename = "/small_region/data/plots/observation in 2020 smoothed 300 by 300.png")
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 2020") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename="/small_region/data/plots/stream.png")
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = Stream)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  ggtitle("Stream") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  coord_equal() 
dev.off()

png(filename="/small_region/data/plots/distance to stream.png")
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = original.dist)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  labs(fill = "distance\nto\nstream") +
  ggtitle("Distance to Stream") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename="/small_region/data/plots/elevation.png")
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = original.DEM_30m)) +
  geom_tile() +
  ggtitle("Elevation") +
  scale_fill_gradientn(colours = terrain.colors(7)) +
  labs(fill = "Elevation\n(DEM_30m)") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal()
dev.off()

# aggregate to 100 by 100
N.100 <- 100
LoG.100 <- 90

Lmatrix <- matrix(0, nrow = N.100, ncol = N)
for (k in 1:N.100) {
  Lmatrix[k, 3*k-2] <- 1/3
  Lmatrix[k, 3*k-1] <- 1/3
  Lmatrix[k, 3*k] <- 1/3
}
Rmatrix <- t(Lmatrix)

# Stream
stream.matrix.85 <- matrix(dat85$Stream, nrow = N, ncol = N, byrow = F)
stream.matrix100.85 <- Lmatrix %*% stream.matrix.85 %*% Rmatrix
stream.vec100.85 <- as.vector(stream.matrix100.85)
stream.vec100.85 <- ceiling(stream.vec100.85)

# Dem30m
dem30m.matrix.85 <- matrix(dat85$original.DEM_30m, nrow = N, ncol = N, byrow = F)
dem30m.matrix100.85 <- Lmatrix %*% dem30m.matrix.85 %*% Rmatrix
dem30m.vec100.85 <- as.vector(dem30m.matrix100.85)

# Trees
trees <- c()
for (j in 1:(TIME+1)) {
  tree.vec300 <- dat[which(dat$year==1984+j), "original.trees"]
  tree.matrix <- matrix(tree.vec300, nrow = N, ncol = N, byrow = F)
  tree.matrix100 <- Lmatrix %*% tree.matrix %*% Rmatrix
  trees <- c(trees, as.vector(tree.matrix100))
}

s1.100 <- rep(rep(c(N.100:1), each = N.100), TIME+1)
s2.100 <- rep(rep(c(1:N.100), N.100), TIME+1)
dem30m <- rep(dem30m.vec100.85, TIME+1)
stream <- rep(stream.vec100.85, TIME+1)
year <- rep(c(1985:2020), each = N.100*N.100)
dat8520.N100 <- data.frame(trees = trees, DEM_30m = dem30m,
                      Stream = stream, year = year, s1 = s1.100, s2 = s2.100)

dat.N100 <- dat8520.N100

# Convert coordinates to index, trees to percentage
dat.N100 <- dat.N100 %>% mutate(trees = round(trees, digits = 0)/100)

index.N100 <- rep(c(1:(N.100 * N.100)), TIME+1)
dat.N100 <- cbind(index.N100, dat.N100)

write.csv(dat.N100, "/small_region/data/dat_8520_small_area_100_by_100_smoothed_integered.csv")

pic_limit <- c(0,1)

data2015 <- dat.N100 %>% filter(year == 2015) %>% mutate(elevation = DEM_30m)
png(filename = "/small_region/data/plots/observation in 2015 with average aggreagation 100 by 100.png")
ggplot(data = data2015, aes(s1*LoG.100, s2*LoG.100, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Observation in 2015 with Average Aggregation 100 by 100") +
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

png(filename = "/small_region/data/plots/stream with average aggregation 100 by 100.png")
ggplot(data = data2015, aes(s1*LoG.100, s2*LoG.100, fill = Stream)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  ggtitle("Stream with Average Aggregation 100 by 100") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  coord_equal()
dev.off()

png(filename = "/small_region/data/plots/elevation with average aggregation 100 by 100.png")
ggplot(data = data2015, aes(s1*LoG.100, s2*LoG.100, fill = elevation)) +
  geom_tile() +
  ggtitle("Elevation with Average Aggregation 100 by 100") +
  labs(fill = "elevation\n(DEM_30m)") +
  theme_bw() +
  scale_fill_gradientn(colours = terrain.colors(7)) + 
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal()
dev.off()

# dat <- read.csv("dat_8520_small_area_100_by_100_smoothed_integered_standardized.csv")

# dat <- read.csv("dat_8520_small_area_100_by_100_smoothed.csv")
