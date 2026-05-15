# draw the plots of the small region selected

library(tidyverse)
library(mvtnorm)
library(FNN)
library(egg)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)

setwd("/small_region/raw_data/plots")

# Load dataset 
dat1985 <- read.csv("/raw_data/1985.csv")
dat1986 <- read.csv("/raw_data/1986.csv")
dat1987 <- read.csv("/raw_data/1987.csv")
dat1988 <- read.csv("/raw_data/1988.csv")
dat1989 <- read.csv("/raw_data/1989.csv")
dat1990 <- read.csv("/raw_data/1990.csv")
dat1991 <- read.csv("/raw_data/1991.csv")
dat1992 <- read.csv("/raw_data/1992.csv")
dat1993 <- read.csv("/raw_data/1993.csv")
dat1994 <- read.csv("/raw_data/1994.csv")
dat1995 <- read.csv("/raw_data/1995.csv")
dat1996 <- read.csv("/raw_data/1996.csv")
dat1997 <- read.csv("/raw_data/1997.csv")
dat1998 <- read.csv("/raw_data/1998.csv")
dat1999 <- read.csv("/raw_data/1999.csv")
dat2000 <- read.csv("/raw_data/2000.csv")
dat2001 <- read.csv("/raw_data/2001.csv")
dat2002 <- read.csv("/raw_data/2002.csv")
dat2003 <- read.csv("/raw_data/2003.csv")
dat2004 <- read.csv("/raw_data/2004.csv")
dat2005 <- read.csv("/raw_data/2005.csv")
dat2006 <- read.csv("/raw_data/2006.csv")
dat2007 <- read.csv("/raw_data/2007.csv")
dat2008 <- read.csv("/raw_data/2008.csv")
dat2009 <- read.csv("/raw_data/2009.csv")
dat2010 <- read.csv("/raw_data/2010.csv")
dat2011 <- read.csv("/raw_data/2011.csv")
dat2012 <- read.csv("/raw_data/2012.csv")
dat2013 <- read.csv("/raw_data/2013.csv")
dat2014 <- read.csv("/raw_data/2014.csv")
dat2015 <- read.csv("/raw_data/2015.csv")
dat2016 <- read.csv("/raw_data/2016.csv")
dat2017 <- read.csv("/raw_data/2017.csv")
dat2018 <- read.csv("/raw_data/2018.csv")
dat2019 <- read.csv("/raw_data/2019.csv")
dat2020 <- read.csv("/raw_data/2020.csv")

# Sort data in 2019 and 2020 by X
dat2019 <- dat2019 %>% arrange(X)
dat2020 <- dat2020 %>% arrange(X)

# Range of longitude and latitude
range(dat1990$Longitude)
range(dat1990$Latitude)

# Check the dimension for 2010
longtd <- dat2010[1000,]$Longitude
nrow(dat2010 %>% filter(Longitude == longtd))
latitd <- dat2010[1000,]$Latitude
nrow(dat2010 %>% filter(Latitude == latitd))

# N (number of grids = N * N)
N <- 300
# Length of grid
LoG <- 30
# number of years 1986-2020
TIME <- 35

# value limit in the plot
pic_limit <- c(0,1)

# Some plots 
# Tree cover in 1985
dat85.all <- dat1985 %>% filter(Longitude >= -100.68 & Longitude <= -100.20 &  
                                  Latitude >= 42.81 & Latitude <= 42.99)
s1.85.all <- rep(c(N1:1), each = N2)
s2.85.all <- rep(c(1:N2), N1)
dat85.all <- dat85.all %>% mutate(s1 = s1.85.all, s2=s2.85.all, trees = trees/100)
png(filename = "/raw_data/plots/observation in 1985 full.png")
ggplot(data = dat85.all, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = c(0,1)) +
  ggtitle(" ") +
  labs(fill = "tree\ncoverage") +
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

# Tree cover in 1995
dat95.all <- dat1995 %>% filter(Longitude >= -100.68 & Longitude <= -100.20 &  
                                  Latitude >= 42.81 & Latitude <= 42.99)
s1.95.all <- rep(c(N1:1), each = N2)
s2.95.all <- rep(c(1:N2), N1)
dat95.all <- dat95.all %>% mutate(s1 = s1.95.all, s2=s2.95.all, trees = trees/100)
png(filename = "/raw_data/plots/observation in 1995 full.png")
ggplot(data = dat95.all, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = c(0,1)) +
  ggtitle(" ") +
  labs(fill = "tree\ncoverage") +
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

# Tree cover in 2005
dat05.all <- dat2005 %>% filter(Longitude >= -100.68 & Longitude <= -100.20 &  
                                  Latitude >= 42.81 & Latitude <= 42.99)
s1.05.all <- rep(c(N1:1), each = N2)
s2.05.all <- rep(c(1:N2), N1)
dat05.all <- dat05.all %>% mutate(s1 = s1.05.all, s2=s2.05.all, trees = trees/100)
png(filename = "/raw_data/plots/observation in 2005 full.png")
ggplot(data = dat05.all, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = c(0,1)) +
  ggtitle(" ") +
  labs(fill = "tree\ncoverage") +
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

# Tree cover in 2015
dat15.all <- dat2015 %>% filter(Longitude >= -100.68 & Longitude <= -100.20 &  
                                  Latitude >= 42.81 & Latitude <= 42.99)
s1.15.all <- rep(c(N1:1), each = N2)
s2.15.all <- rep(c(1:N2), N1)
dat15.all <- dat15.all %>% mutate(s1 = s1.15.all, s2=s2.15.all, trees = trees/100)
png(filename = "/raw_data/plots/observation in 2015 full.png")
ggplot(data = dat15.all, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = c(0,1)) +
  ggtitle(" ") +
  labs(fill = "tree\ncoverage") +
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

# Tree cover in 2020
dat20.all <- dat2020 %>% filter(Longitude >= -100.68 & Longitude <= -100.20 &  
                                  Latitude >= 42.81 & Latitude <= 42.99)
s1.20.all <- rep(c(N1:1), each = N2)
s2.20.all <- rep(c(1:N2), N1)
dat20.all <- dat20.all %>% mutate(s1 = s1.20.all, s2=s2.20.all, trees = trees/100)
png(filename = "/raw_data/plots/observation in 2020 full.png")
ggplot(data = dat20.all, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = c(0,1)) +
  ggtitle(" ") +
  labs(fill = "tree\ncoverage") +
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

# Stream
stream.all <- dat20.all %>% mutate(s1 = s1.20.all, s2=s2.20.all, stream = Stream)
png(filename = "/raw_data/plots/stream full.png")
ggplot(data = stream.all, aes(s1*LoG, s2*LoG, fill = stream)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
  ggtitle(" ") +
  theme_bw()+
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  coord_equal()
dev.off()

# Elevation
elevation.all <- dat20.all %>% mutate(s1 = s1.20.all, s2=s2.20.all, elevation = DEM_30m)
png(filename = "/raw_data/plots/elevation full.png")
ggplot(data = elevation.all, aes(s1*LoG, s2*LoG, fill = elevation)) +
  geom_tile() +
  ggtitle(" ") +
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

# Common setting of the plots
common_fill_obs <- scale_fill_gradient(
  low="white", high="green4", limits = pic_limit,
  name = "tree\ncoverage",
  breaks = seq(0, 1, by = 0.2),
  labels = seq(0, 1, by = 0.2),
  guide = guide_colorbar(
    barheight = unit(30, "cm"),
    barwidth  = unit(1.5, "cm"),
    ticks = TRUE, 
    ticks.colour = "black",
    title.position = "top",
    title.hjust = 0.5,
    title.vjust = 3
  )
)

plot.obs85 <- ggplot(data = dat85.all, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_obs +  # use the common scale here
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 28),
    legend.text  = element_text(size = 24),
    legend.margin = margin(20, 40, 20, 40)) +
  coord_equal()

plot.obs15 <- ggplot(data = dat15.all, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_obs +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal()


title.plot.obs <- theme(plot.title = element_text(size = 28, 
                                                 color = "black", 
                                                 face = "bold", 
                                                 hjust = 0.5, vjust = 1))

plot.obs85.tl <- plot.obs85 + labs(title = "Tree Coverage in the Sandhills Region (1985)") + title.plot.obs
plot.obs15.tl <- plot.obs15 + labs(title = "Tree Coverage in the Sandhills Region (2015)") + title.plot.obs

combined_plot_obs_8515_tl <- plot.obs85.tl + plot.obs15.tl + 
  plot_layout(nrow = 2, ncol = 1, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

# Comparison tree cover between 1985 and 2015
png(filename = "/raw_data/plots/obs_8515_inside_title.png", height = 5000, width = 6000, res = 300)
combined_plot_obs_8515_tl
dev.off()

plot.obs85.str <- plot.obs85.tl + geom_point(data = subset(dat85.all, Stream == 1),
                                               aes(x = s1*LoG, y = s2*LoG),
                                               color = "black",
                                               size = 0.5,
                                               inherit.aes = FALSE)

plot.obs15.str <- plot.obs15.tl + geom_point(data = subset(dat15.all, Stream == 1),
                                             aes(x = s1*LoG, y = s2*LoG),
                                             color = "black",
                                             size = 0.5,
                                             inherit.aes = FALSE)

combined_plot_obs_8515_str_tl <- plot.obs85.str + plot.obs15.str + 
  plot_layout(nrow = 2, ncol = 1, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

# Comparison tree cover between 1985 and 2015
png(filename = "/raw_data/plots/obs_8515_str_inside_title.png", height = 5000, width = 6000, res = 300)
combined_plot_obs_8515_str_tl
dev.off()

obs.diff.8515 <- dat15.all$trees-dat85.all$trees
dat.diff.8515 <- cbind(dat85.all, obs.diff.8515)

# title for single plot
title.plot <- theme(plot.title = element_text(size = 28, 
                                              color = "black", 
                                              face = "bold", 
                                              hjust = 0.5, vjust = 0.5))

# Common setting
plot.obsdiff.8515 <- ggplot(data = dat.diff.8515, aes(s1*LoG, s2*LoG, fill = obs.diff.8515)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = rev(brewer.pal(11, "RdBu")),
    limits = c(-1, 1),
    name = "tree\ncoverage",
    breaks = seq(-1, 1, by = 0.5),
    labels = seq(-1, 1, by = 0.5),
    guide = guide_colorbar(
      barheight = unit(15, "cm"),
      barwidth  = unit(1, "cm"),
      ticks = TRUE, 
      ticks.colour = "black",
      title.position = "top",
      title.hjust = 0.5,
      title.vjust = 3
    )
  ) +
  theme_bw() +
  labs(title = "Change in Tree Coverage in the Sandhills Region (1985 - 2015)") +
  title.plot +
  theme(
    plot.background = element_blank(),
    axis.text = element_text(size = 24),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 28),
    legend.text  = element_text(size = 24),
    legend.margin = margin(20, 40, 20, 40)) +
  coord_equal()

png(filename = "/raw_data/plots/obs_diff_8515_no_stream_inside_title.png", height = 2500, width = 6000, res = 300)
plot.obsdiff.8515
dev.off()

plot.obsdiff <- plot.obsdiff.8515 + geom_point(data = subset(dat.diff.8515, Stream == 1),
                                               aes(x = s1*LoG, y = s2*LoG),
                                               color = "black",
                                               size = 0.5,
                                               inherit.aes = FALSE)

png(filename = "/raw_data/plots/obs_diff_8515_inside_title.png", height = 2500, width = 6000, res = 300)
plot.obsdiff
dev.off()

plot.elevation <- ggplot(data = elevation.all, aes(s1*LoG, s2*LoG, fill = DEM_30m)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = terrain.colors(7),
    name = "m",
    # breaks = seq(0, 1, 0.1),
    # labels = seq(0, 1, 0.1),
    guide = guide_colorbar(
      barheight = unit(13.5, "cm"),
      barwidth  = unit(1, "cm"),
      ticks = TRUE,
      ticks.colour = "black",
      title.position = "top",
      title.hjust = 0.5,
      title.vjust = 2
    )
  ) +
  theme_bw() +
  theme(plot.background = element_blank(),
        axis.text = element_text(size = 24),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 28),
        legend.text  = element_text(size = 24)) +
  coord_equal()

# Plot of elevation
png(filename="/raw_data/plots/elevation_full_no_title.png", height = 2200, width = 6000, res = 300)
plot.elevation
dev.off()

plot.elevation.tl <- plot.elevation + labs(title = "Elevation in the Sandhills Region") + title.plot

png(filename="/raw_data/plots/elevation_full.png", height = 2500, width = 6000, res = 300)
plot.elevation.tl
dev.off()

plot.elevation.str <- plot.elevation.tl + geom_point(data = subset(dat85.all, Stream == 1),
                                                     aes(x = s1*LoG, y = s2*LoG),
                                                     color = "black",
                                                     size = 0.5,
                                                     inherit.aes = FALSE)

# Plot of elevation with stream
png(filename="/raw_data/plots/elevation_str_full.png", height = 2500, width = 6000, res = 300)
plot.elevation.str
dev.off()

plot.stream.blue <- ggplot(data = dat85.all, aes(s1*LoG, s2*LoG, fill = Stream)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
  theme_bw() +
  theme(plot.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  coord_equal() 

# Plot of stream
png(filename="/raw_data/plots/stream_full_blue_no_title.png", height = 2200, width = 5500, res = 300)
plot.stream.blue
dev.off()

plot.stream.blue.tl <- plot.stream.blue + labs(title = "Stream in the Sandhills Region") + title.plot

png(filename="/raw_data/plots/stream_full_blue.png", height = 2500, width = 5500, res = 300)
plot.stream.blue.tl
dev.off()

plot.stream.black <- ggplot(data = dat85.all, aes(s1*LoG, s2*LoG, fill = Stream)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="black", limit = c(0,1)) +
  theme_bw() +
  theme(plot.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  coord_equal() 

png(filename="/raw_data/plots/stream_full_black_no_title.png", height = 2200, width = 5500, res = 300)
plot.stream.black
dev.off()

plot.stream.black.tl <- plot.stream.black + labs(title = "Stream in the Sandhills Region") + title.plot

png(filename="/raw_data/plots/stream_full_black.png", height = 2500, width = 5500, res = 300)
plot.stream.black.tl
dev.off()

# Setup the small region
mid.longitude <- dat15.all[1000*300+1, "Longitude"]
(dat15.all %>% filter(Longitude == mid.longitude))$Latitude[251]
(dat15.all %>% filter(Longitude == mid.longitude))$Latitude[500]
mid.latitude <- dat15.all[1000*300+1, "Latitude"]
(dat15.all %>% filter(Latitude == mid.latitude))$Longitude[201]
(dat15.all %>% filter(Latitude == mid.latitude))$Longitude[500]

dat15 <- dat2015 %>% filter(Longitude >= -100.3346 & Longitude <= -100.254 &  
                              Latitude >= 42.86405 & Latitude <= 42.94464)
nrow(table(dat15$Longitude))
nrow(table(dat15$Latitude))
s1.15 <- rep(c(N:1), each = N)
s2.15 <- rep(c(1:N), N)
dat15 <- dat15 %>% mutate(s1 = s1.15, s2=s2.15, trees = trees/100)
ggplot(data = dat15, aes(s1, s2, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
  labs(title = "observation in 2015") +
  coord_equal()

# Setting up the range of the small area
range_x <- c(-100.3346, -100.254)
range_y <- c(42.86405, 42.94464)

sink(file = "/small_region/data/longitude and latitude.txt")
cat("dimension N: ")
N
cat("longitude range: ")
range_x
cat("latitude range: ")
range_y
sink(file = NULL)

# Extract the small region for each year
subset1985 <- dat1985 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1985 <- subset1985 %>% mutate(year = 1985)

subset1986 <- dat1986 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1986 <- subset1986 %>% mutate(year = 1986)

subset1987 <- dat1987 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1987 <- subset1987 %>% mutate(year = 1987)

subset1988 <- dat1988 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1988 <- subset1988 %>% mutate(year = 1988)

subset1989 <- dat1989 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1989 <- subset1989 %>% mutate(year = 1989)

subset1990 <- dat1990 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1990 <- subset1990 %>% mutate(year = 1990)

subset1991 <- dat1991 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1991 <- subset1991 %>% mutate(year = 1991)

subset1992 <- dat1992 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1992 <- subset1992 %>% mutate(year = 1992)

subset1993 <- dat1993 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1993 <- subset1993 %>% mutate(year = 1993)

subset1994 <- dat1994 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1994 <- subset1994 %>% mutate(year = 1994)

subset1995 <- dat1995 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1995 <- subset1995 %>% mutate(year = 1995)

subset1996 <- dat1996 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1996 <- subset1996 %>% mutate(year = 1996)

subset1997 <- dat1997 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1997 <- subset1997 %>% mutate(year = 1997)

subset1998 <- dat1998 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1998 <- subset1998 %>% mutate(year = 1998)

subset1999 <- dat1999 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1999 <- subset1999 %>% mutate(year = 1999)

subset2000 <- dat2000 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2000 <- subset2000 %>% mutate(year = 2000)

subset2001 <- dat2001 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2001 <- subset2001 %>% mutate(year = 2001)

subset2002 <- dat2002 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2002 <- subset2002 %>% mutate(year = 2002)

subset2003 <- dat2003 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2003 <- subset2003 %>% mutate(year = 2003)

subset2004 <- dat2004 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2004 <- subset2004 %>% mutate(year = 2004)

subset2005 <- dat2005 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2005 <- subset2005 %>% mutate(year = 2005)

subset2006 <- dat2006 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2006 <- subset2006 %>% mutate(year = 2006)

subset2007 <- dat2007 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2007 <- subset2007 %>% mutate(year = 2007)

subset2008 <- dat2008 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2008 <- subset2008 %>% mutate(year = 2008)

subset2009 <- dat2009 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2009 <- subset2009 %>% mutate(year = 2009)

subset2010 <- dat2010 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2010 <- subset2010 %>% mutate(year = 2010)

subset2011 <- dat2011 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2011 <- subset2011 %>% mutate(year = 2011)

subset2012 <- dat2012 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2012 <- subset2012 %>% mutate(year = 2012)

subset2013 <- dat2013 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2013 <- subset2013 %>% mutate(year = 2013)

subset2014 <- dat2014 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2014 <- subset2014 %>% mutate(year = 2014)

subset2015 <- dat2015 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2015 <- subset2015 %>% mutate(year = 2015)

subset2016 <- dat2016 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2016 <- subset2016 %>% mutate(year = 2016)

subset2017 <- dat2017 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2017 <- subset2017 %>% mutate(year = 2017)

subset2018 <- dat2018 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2018 <- subset2018 %>% mutate(year = 2018)

subset2019 <- dat2019 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2019 <- subset2019 %>% 
  dplyr::select(!c(X.1, POINT_X, POINT_Y, Climate_coordinate)) %>% 
  rename(vpdma = X2019_vpdma, tmax = X2019_tmax, tdmea = X2019_tdmea, 
         tmean = X2019_tmean, tmin = X2019_tmin, vpdmi = X2019_vpdmi, 
         ppt = X2019_ppt)
subset2019 <- subset2019 %>% mutate(year = 2019)

subset2020 <- dat2020 %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2020 <- subset2020 %>% 
  dplyr::select(!c(X.1, POINT_X, POINT_Y, Climate_coordinate)) %>% 
  rename(vpdma = X2020_vpdma, tmax = X2020_tmax, tdmea = X2020_tdmea, 
         tmean = X2020_tmean, tmin = X2020_tmin, vpdmi = X2020_vpdmi, 
         ppt = X2020_ppt)
subset2020 <- subset2020 %>% mutate(year = 2020)

df <- NULL
df <- rbind(subset1985, 
            subset1986, subset1987, subset1988, subset1989, subset1990,
            subset1991, subset1992, subset1993, subset1994, subset1995, 
            subset1996, subset1997, subset1998, subset1999, subset2000,
            subset2001, subset2002, subset2003, subset2004, subset2005,
            subset2006, subset2007, subset2008, subset2009, subset2010,
            subset2011, subset2012, subset2013, subset2014, subset2015,
            subset2016, subset2017, subset2018, subset2019, subset2020)

s1 <- rep(rep(c(N:1), each = N), TIME+1)
s2 <- rep(rep(c(1:N), N), TIME+1)
df <- cbind(df, s1, s2)

# Plot of small area
df2015 <- df %>% filter(year == 2015) %>% mutate(tree_coverage = trees/100, elevation = DEM_30m)
png(filename = "/small_region/data/plots/observation in 2015 original 300 by 300.png")
ggplot(data = df2015, aes(s1*LoG, s2*LoG, fill = tree_coverage)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = c(0,1)) +
  labs(title = "observation in 2015 Original 300 by 300") +
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

# Plot of stream in small area
png(filename = "/small_region/data/plots/stream in 2015 original 300 by 300.png")
ggplot(data = df2015, aes(s1*LoG, s2*LoG, fill = Stream)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  labs(title = "Stream 300 by 300") +
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

# Plot of elevation in smal area
png(filename = "/small_region/data/plots/elevation in 2015 original 300 by 300.png")
ggplot(data = df2015, aes(s1*LoG, s2*LoG, fill = elevation)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = c(0,1)) +
  labs(title = "Elevation 300 by 300") +
  labs(fill = "Elevation\n(Dem_30m)") +
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

# write.csv(df, "/small_region/data/dat_8520_small_area_300_by_300_original.csv")

subset1985a <- dat85.all %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset1985a <- subset1985a %>% mutate(year = 1985)

subset2015a <- dat15.all %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2015a <- subset2015a %>% mutate(year = 2015)

subset2020a <- dat20.all %>% 
  filter(Longitude >= range_x[1] & Longitude <= range_x[2]) %>% 
  filter(Latitude >= range_y[1] & Latitude <= range_y[2])
subset2020a <- subset2020a %>% mutate(year = 2020)

box_x <- c(min(subset1985a$s1), max(subset1985a$s1))
box_y <- c(min(subset1985a$s2), max(subset1985a$s2))

plot.obs85.str.box <- plot.obs85.str + 
  geom_point(data = subset(dat85.all, 
                           (s1 >= box_x[1] & s1 <= box_x[2] &
                              s2 >= box_y[1] & s2 <= box_y[2]) &
                             (s1 == box_x[1] | s1 == box_x[2] |
                                s2 == box_y[1] | s2 == box_y[2])),
             aes(x = s1*LoG, y = s2*LoG),
             color = "red",
             size = 0.5,
             inherit.aes = FALSE)

plot.obs15.str.box <- plot.obs15.str + 
  geom_point(data = subset(dat15.all, 
                           (s1 >= box_x[1] & s1 <= box_x[2] &
                              s2 >= box_y[1] & s2 <= box_y[2]) &
                             (s1 == box_x[1] | s1 == box_x[2] |
                                s2 == box_y[1] | s2 == box_y[2])),
             aes(x = s1*LoG, y = s2*LoG),
             color = "red",
             size = 0.5,
             inherit.aes = FALSE)

combined_plot_obs_8515_str_box_tl <- plot.obs85.str.box + plot.obs15.str.box + 
  plot_layout(nrow = 2, ncol = 1, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

png(filename = "/small_region/data/plots/obs_8515_str_box_inside_title.png", height = 5000, width = 6000, res = 300)
combined_plot_obs_8515_str_box_tl
dev.off()

plot.obs85.str.boxgrey <- plot.obs85.str + 
  geom_point(data = subset(dat85.all, 
                           (s1 >= box_x[1] & s1 <= box_x[2] &
                              s2 >= box_y[1] & s2 <= box_y[2]) &
                             (s1 == box_x[1] | s1 == box_x[2] |
                                s2 == box_y[1] | s2 == box_y[2])),
             aes(x = s1*LoG, y = s2*LoG),
             color = "gray50",
             size = 0.5,
             inherit.aes = FALSE)

plot.obs15.str.boxgrey <- plot.obs15.str + 
  geom_point(data = subset(dat15.all, 
                           (s1 >= box_x[1] & s1 <= box_x[2] &
                              s2 >= box_y[1] & s2 <= box_y[2]) &
                             (s1 == box_x[1] | s1 == box_x[2] |
                                s2 == box_y[1] | s2 == box_y[2])),
             aes(x = s1*LoG, y = s2*LoG),
             color = "gray50",
             size = 0.5,
             inherit.aes = FALSE)

combined_plot_obs_8515_str_boxgrey_tl <- plot.obs85.str.boxgrey + plot.obs15.str.boxgrey + 
  plot_layout(nrow = 2, ncol = 1, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

png(filename = "/small_region/data/plots/obs_8515_str_boxgrey_inside_title.png", height = 5000, width = 6000, res = 300)
combined_plot_obs_8515_str_boxgrey_tl
dev.off()

plot.obsdiff.boxgrey <- plot.obsdiff + 
  geom_point(data = subset(dat.diff.8515, 
                           (s1 >= box_x[1] & s1 <= box_x[2] &
                              s2 >= box_y[1] & s2 <= box_y[2]) &
                             (s1 == box_x[1] | s1 == box_x[2] |
                                s2 == box_y[1] | s2 == box_y[2])),
             aes(x = s1*LoG, y = s2*LoG),
             color = "gray50",
             size = 0.5,
             inherit.aes = FALSE)

png(filename = "/small_region/data/plots/obs_diff_8515_boxgrey_inside_title.png", height = 2500, width = 6000, res = 300)
plot.obsdiff.boxgrey
dev.off()

plot.elevation.str.box <- plot.elevation.str +
  geom_point(data = subset(elevation.all, 
                           (s1 >= box_x[1] & s1 <= box_x[2] &
                              s2 >= box_y[1] & s2 <= box_y[2]) &
                             (s1 == box_x[1] | s1 == box_x[2] |
                                s2 == box_y[1] | s2 == box_y[2])),
             aes(x = s1*LoG, y = s2*LoG),
             color = "red",
             size = 0.5,
             inherit.aes = FALSE)

png(filename="/small_region/data/plots/elevation_str_box_full.png", height = 2500, width = 6000, res = 300)
plot.elevation.str.box
dev.off()

plot.elevation.str.boxgrey <- plot.elevation.str +
  geom_point(data = subset(elevation.all, 
                           (s1 >= box_x[1] & s1 <= box_x[2] &
                              s2 >= box_y[1] & s2 <= box_y[2]) &
                             (s1 == box_x[1] | s1 == box_x[2] |
                                s2 == box_y[1] | s2 == box_y[2])),
             aes(x = s1*LoG, y = s2*LoG),
             color = "gray50",
             size = 0.5,
             inherit.aes = FALSE)

png(filename="/small_region/data/plots/elevation_str_boxgrey_full.png", height = 2500, width = 6000, res = 300)
plot.elevation.str.boxgrey
dev.off()

# Aggregate the small area to size 100 * 100
# Setup the left matrix and right matrix
Lmatrix <- matrix(0, nrow = 100, ncol = N)
for (k in 1:100) {
  Lmatrix[k, 3*k-2] <- 1/3
  Lmatrix[k, 3*k-1] <- 1/3
  Lmatrix[k, 3*k] <- 1/3
}
Rmatrix <- t(Lmatrix)

# Aggregate the stream
stream.matrix.85 <- matrix(subset1985$Stream, nrow = N, ncol = N, byrow = F)
stream.matrix100.85 <- Lmatrix %*% stream.matrix.85 %*% Rmatrix
stream.vec100.85 <- as.vector(stream.matrix100.85)
stream.vec100.85 <- ceiling(stream.vec100.85)

# Aggregate Dem30m
dem30m.matrix.85 <- matrix(subset1985$DEM_30m, nrow = N, ncol = N, byrow = F)
dem30m.matrix100.85 <- Lmatrix %*% dem30m.matrix.85 %*% Rmatrix
dem30m.vec100.85 <- as.vector(dem30m.matrix100.85)

# Aggregate the slope
slope30m.matrix.85 <- matrix(subset1985$Slope_30m, nrow = N, ncol = N, byrow = F)
slope30m.matrix100.85 <- Lmatrix %*% slope30m.matrix.85 %*% Rmatrix
slope30m.vec100.85 <- as.vector(slope30m.matrix100.85)

# Aggregate tree cover in each year
tree.matrix.85 <- matrix(subset1985$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.85 <- Lmatrix %*% tree.matrix.85 %*% Rmatrix
tree.vec100.85 <- as.vector(tree.matrix100.85)

tree.matrix.86 <- matrix(subset1986$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.86 <- Lmatrix %*% tree.matrix.86 %*% Rmatrix
tree.vec100.86 <- as.vector(tree.matrix100.86)

tree.matrix.87 <- matrix(subset1987$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.87 <- Lmatrix %*% tree.matrix.87 %*% Rmatrix
tree.vec100.87 <- as.vector(tree.matrix100.87)

tree.matrix.88 <- matrix(subset1988$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.88 <- Lmatrix %*% tree.matrix.88 %*% Rmatrix
tree.vec100.88 <- as.vector(tree.matrix100.88)

tree.matrix.89 <- matrix(subset1989$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.89 <- Lmatrix %*% tree.matrix.89 %*% Rmatrix
tree.vec100.89 <- as.vector(tree.matrix100.89)

tree.matrix.90 <- matrix(subset1990$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.90 <- Lmatrix %*% tree.matrix.90 %*% Rmatrix
tree.vec100.90 <- as.vector(tree.matrix100.90)

tree.matrix.91 <- matrix(subset1991$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.91 <- Lmatrix %*% tree.matrix.91 %*% Rmatrix
tree.vec100.91 <- as.vector(tree.matrix100.91)

tree.matrix.92 <- matrix(subset1992$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.92 <- Lmatrix %*% tree.matrix.92 %*% Rmatrix
tree.vec100.92 <- as.vector(tree.matrix100.92)

tree.matrix.93 <- matrix(subset1993$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.93 <- Lmatrix %*% tree.matrix.93 %*% Rmatrix
tree.vec100.93 <- as.vector(tree.matrix100.93)

tree.matrix.94 <- matrix(subset1994$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.94 <- Lmatrix %*% tree.matrix.94 %*% Rmatrix
tree.vec100.94 <- as.vector(tree.matrix100.94)

tree.matrix.95 <- matrix(subset1995$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.95 <- Lmatrix %*% tree.matrix.95 %*% Rmatrix
tree.vec100.95 <- as.vector(tree.matrix100.95)

tree.matrix.96 <- matrix(subset1996$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.96 <- Lmatrix %*% tree.matrix.96 %*% Rmatrix
tree.vec100.96 <- as.vector(tree.matrix100.96)

tree.matrix.97 <- matrix(subset1997$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.97 <- Lmatrix %*% tree.matrix.97 %*% Rmatrix
tree.vec100.97 <- as.vector(tree.matrix100.97)

tree.matrix.98 <- matrix(subset1998$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.98 <- Lmatrix %*% tree.matrix.98 %*% Rmatrix
tree.vec100.98 <- as.vector(tree.matrix100.98)

tree.matrix.99 <- matrix(subset1999$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.99 <- Lmatrix %*% tree.matrix.99 %*% Rmatrix
tree.vec100.99 <- as.vector(tree.matrix100.99)

tree.matrix.00 <- matrix(subset2000$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.00 <- Lmatrix %*% tree.matrix.00 %*% Rmatrix
tree.vec100.00 <- as.vector(tree.matrix100.00)

tree.matrix.01 <- matrix(subset2001$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.01 <- Lmatrix %*% tree.matrix.01 %*% Rmatrix
tree.vec100.01 <- as.vector(tree.matrix100.01)

tree.matrix.02 <- matrix(subset2002$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.02 <- Lmatrix %*% tree.matrix.02 %*% Rmatrix
tree.vec100.02 <- as.vector(tree.matrix100.02)

tree.matrix.03 <- matrix(subset2003$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.03 <- Lmatrix %*% tree.matrix.03 %*% Rmatrix
tree.vec100.03 <- as.vector(tree.matrix100.03)

tree.matrix.04 <- matrix(subset2004$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.04 <- Lmatrix %*% tree.matrix.04 %*% Rmatrix
tree.vec100.04 <- as.vector(tree.matrix100.04)

tree.matrix.05 <- matrix(subset2005$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.05 <- Lmatrix %*% tree.matrix.05 %*% Rmatrix
tree.vec100.05 <- as.vector(tree.matrix100.05)

tree.matrix.06 <- matrix(subset2006$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.06 <- Lmatrix %*% tree.matrix.06 %*% Rmatrix
tree.vec100.06 <- as.vector(tree.matrix100.06)

tree.matrix.07 <- matrix(subset2007$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.07 <- Lmatrix %*% tree.matrix.07 %*% Rmatrix
tree.vec100.07 <- as.vector(tree.matrix100.07)

tree.matrix.08 <- matrix(subset2008$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.08 <- Lmatrix %*% tree.matrix.08 %*% Rmatrix
tree.vec100.08 <- as.vector(tree.matrix100.08)

tree.matrix.09 <- matrix(subset2009$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.09 <- Lmatrix %*% tree.matrix.09 %*% Rmatrix
tree.vec100.09 <- as.vector(tree.matrix100.09)

tree.matrix.10 <- matrix(subset2010$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.10 <- Lmatrix %*% tree.matrix.10 %*% Rmatrix
tree.vec100.10 <- as.vector(tree.matrix100.10)

tree.matrix.11 <- matrix(subset2011$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.11 <- Lmatrix %*% tree.matrix.11 %*% Rmatrix
tree.vec100.11 <- as.vector(tree.matrix100.11)

tree.matrix.12 <- matrix(subset2012$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.12 <- Lmatrix %*% tree.matrix.12 %*% Rmatrix
tree.vec100.12 <- as.vector(tree.matrix100.12)

tree.matrix.13 <- matrix(subset2013$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.13 <- Lmatrix %*% tree.matrix.13 %*% Rmatrix
tree.vec100.13 <- as.vector(tree.matrix100.13)

tree.matrix.14 <- matrix(subset2014$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.14 <- Lmatrix %*% tree.matrix.14 %*% Rmatrix
tree.vec100.14 <- as.vector(tree.matrix100.14)

tree.matrix.15 <- matrix(subset2015$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.15 <- Lmatrix %*% tree.matrix.15 %*% Rmatrix
tree.vec100.15 <- as.vector(tree.matrix100.15)

tree.matrix.16 <- matrix(subset2016$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.16 <- Lmatrix %*% tree.matrix.16 %*% Rmatrix
tree.vec100.16 <- as.vector(tree.matrix100.16)

tree.matrix.17 <- matrix(subset2017$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.17 <- Lmatrix %*% tree.matrix.17 %*% Rmatrix
tree.vec100.17 <- as.vector(tree.matrix100.17)

tree.matrix.18 <- matrix(subset2018$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.18 <- Lmatrix %*% tree.matrix.18 %*% Rmatrix
tree.vec100.18 <- as.vector(tree.matrix100.18)

tree.matrix.19 <- matrix(subset2019$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.19 <- Lmatrix %*% tree.matrix.19 %*% Rmatrix
tree.vec100.19 <- as.vector(tree.matrix100.19)

tree.matrix.20 <- matrix(subset2020$trees, nrow = N, ncol = N, byrow = F)
tree.matrix100.20 <- Lmatrix %*% tree.matrix.20 %*% Rmatrix
tree.vec100.20 <- as.vector(tree.matrix100.20)

# Build a dataset
N.100 <- 100
s1.100 <- rep(rep(c(N.100:1), each = N.100), TIME+1)
s2.100 <- rep(rep(c(1:N.100), N.100), TIME+1)
trees <- c(tree.vec100.85, tree.vec100.86, tree.vec100.87, tree.vec100.88, tree.vec100.89,
           tree.vec100.90, tree.vec100.91, tree.vec100.92, tree.vec100.93, tree.vec100.94,
           tree.vec100.95, tree.vec100.96, tree.vec100.97, tree.vec100.98, tree.vec100.99,
           tree.vec100.00, tree.vec100.01, tree.vec100.02, tree.vec100.03, tree.vec100.04,
           tree.vec100.05, tree.vec100.06, tree.vec100.07, tree.vec100.08, tree.vec100.09,
           tree.vec100.10, tree.vec100.11, tree.vec100.12, tree.vec100.13, tree.vec100.14,
           tree.vec100.15, tree.vec100.16, tree.vec100.17, tree.vec100.18, tree.vec100.19,
           tree.vec100.20)
dem30m <- rep(dem30m.vec100.85, TIME+1)
slope30m <- rep(slope30m.vec100.85, TIME+1)
stream <- rep(stream.vec100.85, TIME+1)
year <- rep(c(1985:2020), each = N.100*N.100)
dat8520 <- data.frame(trees = trees, DEM_30m = dem30m, Slope_30m = slope30m, 
                      Stream = stream, year = year, s1 = s1.100, s2 = s2.100)

dat <- dat8520

# Convert coordinates to index, trees to percentage
dat <- dat %>% mutate(trees = round(trees, digits = 0)/100)

index <- rep(c(1:10000), TIME+1)
dat <- cbind(index, dat)

# write.csv(dat, "/small_region/data/dat_8520_small_area_100_by_100_average_aggregate_from_300_by_300.csv")

LoG.100 <- LoG * 3

# Some plots
data2015 <- dat %>% filter(year == 2015) %>% mutate(elevation = DEM_30m)
# png(filename = "/small_region/data/plots/observation in 2015 with average aggregation 100 by 100.png")
ggplot(data = data2015, aes(s1*LoG.100, s2*LoG.100, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  labs(title = "Observation in 2015 with Average Aggregation 100 by 100") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal()
# dev.off()

# png(filename = "/small_region/data/plots/stream with average aggregation 100 by 100.png")
ggplot(data = data2015, aes(s1*LoG.100, s2*LoG.100, fill = Stream)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  labs(title = "Stream with Average Aggregation 100 by 100") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  coord_equal()
# dev.off()

# png(filename = "/small_region/data/plots/elevation with average aggregation 100 by 100.png")
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
# dev.off()

# png(filename = "/small_region/data/plots/slope with average aggregation 100 by 100.png")
ggplot(data = data2015, aes(s1*LoG.100, s2*LoG.100, fill = Slope_30m)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  labs(title = "Slope with Average Aggregation 100 by 100") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal()
# dev.off()