# Data cleaning: smooth, calculate distance to stream of each grid, standardize the covariates, N=100

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
N <- 100
# Length of grid
LoG <- 90
# number of years 1986-2020
TIME <- 35
# size in binomial
SIZE <- 100

setwd("/small_region/data/")
data_small_region_dir <- "/small_region/data"
plot_small_region_dir <- "/small_region/data/plots"

dat.original <- read.csv(file.path(data_small_region_dir, "dat_8520_small_area_100_by_100_smoothed_integered.csv"))

dat <- dat.original[, c("year", "s1", "s2", "Stream", "DEM_30m", "trees")]
index <- rep(1:(N*N), TIME+1)

dat <- cbind(index, dat)

# Add a covariate as distance to the nearest stream
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

# Standardize dist, dem_30, Stream
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

write.csv(dat, file.path(data_small_region_dir, "dat_8520_small_area_100_by_100_smoothed_integered_standardized.csv"))

# graph
dat85 <- dat %>% filter(year == 1985)
png(filename = file.path(plot_small_region_dir, "observation in 1985 smoothed 100 by 100.png"))
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
png(filename = file.path(plot_small_region_dir, "observation in 1990 smoothed 100 by 100.png"))
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
png(filename = file.path(plot_small_region_dir, "observation in 1995 smoothed 100 by 100.png"))
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
png(filename = file.path(plot_small_region_dir, "observation in 2000 smoothed 100 by 100.png"))
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
png(filename = file.path(plot_small_region_dir, "observation in 2005 smoothed 100 by 100.png"))
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
png(filename = file.path(plot_small_region_dir, "observation in 2010 smoothed 100 by 100.png"))
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
png(filename = file.path(plot_small_region_dir, "observation in 2015 smoothed 100 by 100.png"))
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
png(filename = file.path(plot_small_region_dir, "observation in 2020 smoothed 100 by 100.png"))
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

png(filename = file.path(plot_small_region_dir, "stream 100 by 100.png"))
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

png(filename = file.path(plot_small_region_dir, "distance to stream 100 by 100.png"))
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

png(filename = file.path(plot_small_region_dir, "elevation 100 by 100.png"))
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

