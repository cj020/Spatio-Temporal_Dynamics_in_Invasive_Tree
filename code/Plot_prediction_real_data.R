# 9/11/25 plots
# 10/3/25 add GAM, RT
# 1/11/26 add INLA 

# install.packages("ggplot2", dependencies = TRUE, repos = "https://cloud.r-project.org")
# install.packages("patchwork", dependencies = TRUE, repos = "https://cloud.r-project.org")


# Set the path to the local project directory
working_dir <- "C:/Users/YourName/Documents/Spatio-Temporal_Dynamics_in_Invasive_Tree"

# Set working directory
setwd(working_dir)

# Get the working directory
wd <- getwd()

plot_dir <- file.path(wd, "real_data", "prediction_summary", "plots")
uncertainty_dir <- file.path(wd, "real_data", "prediction_summary", "MSE_and_prediction_uncertainty")

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
if (!dir.exists(uncertainty_dir)) {
  dir.create(uncertainty_dir, recursive = TRUE)
}

required_packages <- c(
  "fields",
  "RColorBrewer",
  "mvtnorm",
  "fBasics",
  "coda",
  "tidyverse",
  "egg",
  "ggplot2",
  "patchwork",
  "cowplot"
)
missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(fields)
library(RColorBrewer)
library(mvtnorm)
# library(gridExtra)
library(fBasics)
library(coda)
library(tidyverse)
library(egg)
library(ggplot2)
library(patchwork)
library(cowplot)
library(grid)

source(file.path(wd, "code", "Functions.R"))
PDE_dat <- read.csv(file.path(wd, "real_data", "PDE", "prediction", "prediction_real_data_300_by_300_homogenization_diffu_logis_growth_ZIB.csv"))
GAM_dat <- read.csv(file.path(wd, "real_data", "GAM", "prediction_with_GAM_real_data_300_by_300.csv"))
RT_dat <- read.csv(file.path(wd, "real_data", "RT", "prediction_with_RT_real_data_300_by_300.csv"))
INLA_dat <- read.csv(file.path(wd, "real_data", "INLA", "prediction_with_INLA_result_real_data_300_by_300.csv"))
dat <- PDE_dat
GAM.pred <- GAM_dat[,"prediction"]
INLA.pred <- INLA_dat[,"prediction_sp"]
RT.original.pred <- RT_dat[,"prediction.original.rt"]
RT.pruned.pred <- RT_dat[,"prediction.pruned.rt"]
dat <- cbind(dat, GAM.pred, RT.original.pred, RT.pruned.pred, INLA.pred)

# number of grids N by N
N <- 300
# rescale epsilon
EPSILON <- 15
# number of grids after rescaling
M <- N/EPSILON
# Length of grid
LoG <- 30
# time scale of 1 year
TSCALE <- 1
# number of years
TIME <- 35
# number of years as training data
TRAIN <- 30
# Time steps within 1 year
STEP <- 5
# delta t: length of time of 1 step
DELTA_T <- TSCALE/STEP
# number of covariates
K <- 2
# size in binomial
SIZE <- 100

DIFFUSION <- 1
GROWTH <- 1
GAMMA <- 1

pic_limit <- c(0,1)

dat85 <- dat %>% filter(year == 1985)
dat90 <- dat %>% filter(year == 1990)
dat95 <- dat %>% filter(year == 1995)
dat00 <- dat %>% filter(year == 2000)
dat05 <- dat %>% filter(year == 2005)
dat08 <- dat %>% filter(year == 2008)
dat10 <- dat %>% filter(year == 2010)
dat15 <- dat %>% filter(year == 2015)
dat16 <- dat %>% filter(year == 2016)
dat18 <- dat %>% filter(year == 2018)
dat20 <- dat %>% filter(year == 2020)

# title for single plot
title.plot <- theme(plot.title = element_text(size = 20, 
                                              color = "black", 
                                              face = "bold", 
                                              hjust = 0.5, vjust = 0.5))

# initial values
plot.observation0 <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit, 
                      name = "tree\ncoverage",
                      breaks = seq(0, 1, by = 0.2),
                      labels = seq(0, 1, by = 0.2),
                      guide = guide_colorbar(
                        barheight = unit(13.5, "cm"),    
                        barwidth  = unit(1, "cm"),
                        ticks = TRUE, 
                        ticks.colour = "black",
                        title.position = "top",
                        title.hjust = 0.5,
                        title.vjust = 2,
                        label = TRUE
                        )) +
  theme_bw() +
  theme(plot.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14)) +
  coord_equal()

# png(filename = file.path(plot_dir, "initial_values_no_title.png"), height = 2000, width = 2400, res = 300)
# plot.observation0
# dev.off()

plot.observation0.tl <- plot.observation0 + labs(title = "Tree Coverage (1985)") + title.plot
png(filename = file.path(plot_dir, "initial_values.png"), height = 2200, width = 2400, res = 300)
plot.observation0.tl
dev.off()

# elevation
plot.elevation <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = original.DEM_30m)) +
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
        axis.text = element_text(size = 14),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14)) +
  coord_equal()

# png(filename = file.path(plot_dir, "elevation_no_title.png"), height = 2000, width = 2400, res = 300)
# plot.elevation
# dev.off()

plot.elevation.tl <- plot.elevation + labs(title = "Elevation") + title.plot

png(filename = file.path(plot_dir, "elevation.png"), height = 2200, width = 2400, res = 300)
plot.elevation.tl
dev.off()

# stream
plot.stream <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = Stream)) +
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

# png(filename = file.path(plot_dir, "stream_no_title.png"), height = 2000, width = 2000, res = 300)
# plot.stream
# dev.off()

plot.stream.tl <- plot.stream + labs(title = "Stream") + title.plot

png(filename = file.path(plot_dir, "stream.png"), height = 2200, width = 2000, res = 300)
plot.stream.tl
dev.off()

# dist
plot.dist <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = original.dist)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", 
                      name = "m", 
                      # breaks = seq(0, 7000, by = 1000),
                      # labels = seq(0, 7000, by = 1000),
                      guide = guide_colorbar(
                        barheight = unit(13.5, "cm"),    
                        barwidth  = unit(1, "cm"),
                        ticks = TRUE, 
                        ticks.colour = "black",
                        title.position = "top",
                        title.hjust = 0.5,
                        title.vjust = 2,
                        label = TRUE
                        )) +
  theme_bw() +
  theme(plot.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14)) +
  coord_equal()

# png(filename = file.path(plot_dir, "distance_to_stream_no_title.png"), height = 2000, width = 2400, res = 300)
# plot.dist
# dev.off()

plot.dist.tl <- plot.dist + labs(title = "Distance to Stream") + title.plot

png(filename = file.path(plot_dir, "distance_to_stream.png"), height = 2200, width = 2400, res = 300)
plot.dist.tl
dev.off()

# training data
common_fill_training <- scale_fill_gradient(
  low="white", high="green4", limits=pic_limit,
  name = "tree\ncoverage",
  breaks = seq(0, 1, by = 0.2),
  labels = seq(0, 1, by = 0.2),
  guide = guide_colorbar(
    barheight = unit(60, "cm"),
    barwidth  = unit(2, "cm"),
    ticks = TRUE, 
    ticks.colour = "black",
    title.position = "top",
    title.hjust = 0.5,
    title.vjust = 3
    )
)

plot.obs85 <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_training +  # use the common scale here
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

plot.pde85 <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_training +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  )  + 
  coord_equal()

plot.gam85 <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.INLA85 <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.rt85 <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = RT.pruned.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.obs90 <- ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.pde90 <- ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.gam90 <- ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.INLA90 <- ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.rt90 <- ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = RT.pruned.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.obs95 <- ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.pde95 <- ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.gam95 <- ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.INLA95 <- ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.rt95 <- ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = RT.pruned.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.obs00 <- ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.pde00 <- ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.gam00 <- ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.INLA00 <- ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.rt00 <- ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = RT.pruned.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.obs05 <- ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.pde05 <- ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.gam05 <- ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.INLA05 <- ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.rt05 <- ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = RT.pruned.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.obs10 <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.pde10 <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.gam10 <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.INLA10 <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.rt10 <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = RT.pruned.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.obs15 <- ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.pde15 <- ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.gam15 <- ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.INLA15 <- ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.rt15 <- ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = RT.pruned.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

# combined_plot_obs_PDE_GAM_RT_training <- plot.obs0 + plot.pde0 + plot.gam0 + plot.rt0 +
#   plot.obs2 + plot.pde2 + plot.gam2 + plot.rt2 +
#   plot.obs4 + plot.pde4 + plot.gam4 + plot.rt4 +
#   plot.obs6 + plot.pde6 + plot.gam6 + plot.rt6 +
#   plot.obs8 + plot.pde8 + plot.gam8 + plot.rt8 + 
#   plot_layout(nrow = 5, ncol = 4, guides = "collect") & 
#   theme(axis.text = element_text(size = 24),
#         plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

combined_plot_obs_PDE_training <- plot.obs85 + plot.pde85 +  
  plot.obs90 + plot.pde90 + 
  plot.obs95 + plot.pde95 + 
  plot.obs00 + plot.pde00 + 
  plot.obs05 + plot.pde05 + 
  plot.obs10 + plot.pde10 + 
  plot.obs15 + plot.pde15 + 
  plot_layout(nrow = 7, ncol = 2, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

# png(filename = file.path(plot_dir, "obs_PDE_training_no_annotation.png"), height = 14700, width = 5000, res = 300)
# combined_plot_obs_PDE_training
# dev.off()

# column labels as tiny ggplots
# col_labels_training <- c(
#   " ",
#   "Observation", 
#   "Estimates with\nPDE-Based Method", 
#   "Estimates with\nGAM", 
#   "Estimates with\nRT"
# )

col_labels_training <- c(
  " ",
  "Observation", 
  "Estimates with\nPDE-Based Method")


col_label_tr <- lapply(col_labels_training, function(txt) {
  ggplot() +
    geom_text(aes(x = 0, y = 0, label = txt), size = 16, fontface = "bold") +
    theme_void() + coord_fixed()
})

# Row labels as tiny ggplots
# row_labels_training <- c("T = 0", "T = 2", "T = 4", "T = 6", "T = 8")

row_labels_training <- c("1985", "1990", "1995", "2000", "2005", "2010", "2015")
row_label_tr <- lapply(row_labels_training, function(txt) {
  ggplot() +
    geom_text(aes(x = 0, y = 0, label = txt), size = 16, fontface = "bold", hjust = 1) +
    theme_void() + coord_fixed()
})

axis_theme_training <- theme(axis.text = element_text(size = 24),
                             plot.margin = margin(t = 10, r = 10, b = 10, l = 10))


plot.obs85.lb <- plot.obs85 + axis_theme_training
plot.pde85.lb <- plot.pde85 + axis_theme_training

plot.obs90.lb <- plot.obs90 + axis_theme_training
plot.pde90.lb <- plot.pde90 + axis_theme_training

plot.obs95.lb <- plot.obs95 + axis_theme_training
plot.pde95.lb <- plot.pde95 + axis_theme_training

plot.obs00.lb <- plot.obs00 + axis_theme_training
plot.pde00.lb <- plot.pde00 + axis_theme_training

plot.obs05.lb <- plot.obs05 + axis_theme_training
plot.pde05.lb <- plot.pde05 + axis_theme_training

plot.obs10.lb <- plot.obs10 + axis_theme_training
plot.pde10.lb <- plot.pde10 + axis_theme_training

plot.obs15.lb <- plot.obs15 + axis_theme_training
plot.pde15.lb <- plot.pde15 + axis_theme_training


combined_training_with_labels <- col_label_tr[[1]] + col_label_tr[[2]] + 
  col_label_tr[[3]] + 
  row_label_tr[[1]] + plot.obs85.lb + plot.pde85.lb + 
  row_label_tr[[2]] + plot.obs90.lb + plot.pde90.lb + 
  row_label_tr[[3]] + plot.obs95.lb + plot.pde95.lb + 
  row_label_tr[[4]] + plot.obs00.lb + plot.pde00.lb + 
  row_label_tr[[5]] + plot.obs05.lb + plot.pde05.lb + 
  row_label_tr[[6]] + plot.obs10.lb + plot.pde10.lb + 
  row_label_tr[[7]] + plot.obs15.lb + plot.pde15.lb + 
  plot_layout(nrow = 8, ncol = 3, guides = "collect") 

# png(filename = file.path(plot_dir, "obs_PDE_training_with_labels.png"), height = 16500, width = 7500, res = 300)
# combined_training_with_labels
# dev.off()

title.plot.tr <- theme(plot.title = element_text(size = 28, 
                                                 color = "black", 
                                                 face = "bold", 
                                                 hjust = 0.5, vjust = 1))

plot.obs85.tl <- plot.obs85 + labs(title = "Observation\n(1985)") + title.plot.tr
plot.pde85.tl <- plot.pde85 + labs(title = "Estimates with\nPDE-Based Method (1985)") + title.plot.tr
plot.gam85.tl <- plot.gam85 + labs(title = "Estimates with\nGAM (1985)") + title.plot.tr
plot.INLA85.tl <- plot.INLA85 + labs(title = "Estimates with\nINLA (1985)") + title.plot.tr
plot.rt85.tl <- plot.rt85 + labs(title = "Estimates with\nRT (1985)") + title.plot.tr

plot.obs90.tl <- plot.obs90 + labs(title = "Observation\n(1990)") + title.plot.tr
plot.pde90.tl <- plot.pde90 + labs(title = "Estimates with\nPDE-Based Method (1990)") + title.plot.tr
plot.gam90.tl <- plot.gam90 + labs(title = "Estimates with\nGAM (1990)") + title.plot.tr
plot.INLA90.tl <- plot.INLA90 + labs(title = "Estimates with\nINLA (1990)") + title.plot.tr
plot.rt90.tl <- plot.rt90 + labs(title = "Estimates with\nRT (1990)") + title.plot.tr

plot.obs95.tl <- plot.obs95 + labs(title = "Observation\n(1995)") + title.plot.tr
plot.pde95.tl <- plot.pde95 + labs(title = "Estimates with\nPDE-Based Method (1995)") + title.plot.tr
plot.gam95.tl <- plot.gam95 + labs(title = "Estimates with\nGAM (1995)") + title.plot.tr
plot.INLA95.tl <- plot.INLA95 + labs(title = "Estimates with\nINLA (1995)") + title.plot.tr
plot.rt95.tl <- plot.rt95 + labs(title = "Estimates with\nRT (1995)") + title.plot.tr

plot.obs00.tl <- plot.obs00 + labs(title = "Observation\n(2000)") + title.plot.tr
plot.pde00.tl <- plot.pde00 + labs(title = "Estimates with\nPDE-Based Method (2000)") + title.plot.tr
plot.gam00.tl <- plot.gam00 + labs(title = "Estimates with\nGAM (2000)") + title.plot.tr
plot.INLA00.tl <- plot.INLA00 + labs(title = "Estimates with\nINLA (2000)") + title.plot.tr
plot.rt00.tl <- plot.rt00 + labs(title = "Estimates with\nRT (2000)") + title.plot.tr

plot.obs05.tl <- plot.obs05 + labs(title = "Observation\n(2005)") + title.plot.tr
plot.pde05.tl <- plot.pde05 + labs(title = "Estimates with\nPDE-Based Method (2005)") + title.plot.tr
plot.gam05.tl <- plot.gam05 + labs(title = "Estimates with\nGAM (2005)") + title.plot.tr
plot.INLA05.tl <- plot.INLA05 + labs(title = "Estimates with\nINLA (2005)") + title.plot.tr
plot.rt05.tl <- plot.rt05 + labs(title = "Estimates with\nRT (2005)") + title.plot.tr

plot.obs10.tl <- plot.obs10 + labs(title = "Observation\n(2010)") + title.plot.tr
plot.pde10.tl <- plot.pde10 + labs(title = "Estimates with\nPDE-Based Method (2010)") + title.plot.tr
plot.gam10.tl <- plot.gam10 + labs(title = "Estimates with\nGAM (2010)") + title.plot.tr
plot.INLA10.tl <- plot.INLA10 + labs(title = "Estimates with\nINLA (2010)") + title.plot.tr
plot.rt10.tl <- plot.rt10 + labs(title = "Estimates with\nRT (2010)") + title.plot.tr

plot.obs15.tl <- plot.obs15 + labs(title = "Observation\n(2015)") + title.plot.tr
plot.pde15.tl <- plot.pde15 + labs(title = "Estimates with\nPDE-Based Method (2015)") + title.plot.tr
plot.gam15.tl <- plot.gam15 + labs(title = "Estimates with\nGAM (2015)") + title.plot.tr
plot.INLA15.tl <- plot.INLA15 + labs(title = "Estimates with\nINLA (2015)") + title.plot.tr
plot.rt15.tl <- plot.rt15 + labs(title = "Estimates with\nRT (2015)") + title.plot.tr

combined_plot_obs_PDE_training_tl <- plot.obs85.tl + plot.pde85.tl + 
  plot.obs90.tl + plot.pde90.tl + 
  plot.obs95.tl + plot.pde95.tl + 
  plot.obs00.tl + plot.pde00.tl + 
  plot.obs05.tl + plot.pde05.tl + 
  plot.obs10.tl + plot.pde10.tl + 
  plot.obs15.tl + plot.pde15.tl + 
  plot_layout(nrow = 7, ncol = 2, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

png(filename = file.path(plot_dir, "obs_PDE_training_inside_title.png"), height = 18900, width = 6000, res = 300)
combined_plot_obs_PDE_training_tl
dev.off()

combined_plot_obs_PDE_GAM_INLA_RT_training_tl <- plot.obs85.tl + plot.pde85.tl + plot.gam85.tl + plot.INLA85.tl + plot.rt85.tl +
  plot.obs90.tl + plot.pde90.tl + plot.gam90.tl + plot.INLA90.tl + plot.rt90.tl +
  plot.obs95.tl + plot.pde95.tl + plot.gam95.tl + plot.INLA95.tl + plot.rt95.tl +
  plot.obs00.tl + plot.pde00.tl + plot.gam00.tl + plot.INLA00.tl + plot.rt00.tl +
  plot.obs05.tl + plot.pde05.tl + plot.gam05.tl + plot.INLA05.tl + plot.rt05.tl +
  plot.obs10.tl + plot.pde10.tl + plot.gam10.tl + plot.INLA10.tl + plot.rt10.tl +
  plot.obs15.tl + plot.pde15.tl + plot.gam15.tl + plot.INLA15.tl + plot.rt15.tl +
  plot_layout(nrow = 7, ncol = 5, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

png(filename = file.path(plot_dir, "obs_PDE_GAM_INLA_RT_training_inside_title.png"), height = 18900, width = 13500, res = 300)
combined_plot_obs_PDE_GAM_INLA_RT_training_tl
dev.off()

combined_plot_obs_PDE_training_tl_4yrs <- plot.obs85.tl + plot.pde85.tl + 
  plot.obs95.tl + plot.pde95.tl + 
  plot.obs05.tl + plot.pde05.tl + 
  plot.obs15.tl + plot.pde15.tl + 
  plot_layout(nrow = 4, ncol = 2, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

png(filename = file.path(plot_dir, "obs_PDE_training_inside_title_4yrs.png"), height = 10800, width = 6000, res = 300)
combined_plot_obs_PDE_training_tl_4yrs
dev.off()

combined_plot_obs_PDE_GAM_INLA_RT_training_tl_4yrs <- plot.obs85.tl + plot.pde85.tl + plot.gam85.tl + plot.INLA85.tl + plot.rt85.tl +
  plot.obs95.tl + plot.pde95.tl + plot.gam95.tl + plot.INLA95.tl + plot.rt95.tl +
  plot.obs05.tl + plot.pde05.tl + plot.gam05.tl + plot.INLA05.tl + plot.rt05.tl +
  plot.obs15.tl + plot.pde15.tl + plot.gam15.tl + plot.INLA15.tl + plot.rt15.tl +
  plot_layout(nrow = 4, ncol = 5, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

png(filename = file.path(plot_dir, "obs_PDE_GAM_INLA_RT_training_inside_title_4yrs.png"), height = 10800, width = 13500, res = 300)
combined_plot_obs_PDE_GAM_INLA_RT_training_tl_4yrs
dev.off()

# test data
common_fill_test <- scale_fill_gradient(
  low="white", high="green4", limits=pic_limit,
  name = "tree\ncoverage",
  breaks = seq(0, 1, by = 0.2),
  labels = seq(0, 1, by = 0.2),
  guide = guide_colorbar(
    barheight = unit(40, "cm"),
    barwidth  = unit(2, "cm"),
    ticks = TRUE, 
    ticks.colour = "black",
    title.position = "top",
    title.hjust = 0.5,
    title.vjust = 3
    )
)

plot.obs16 <- ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_test +  # use the common scale here
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

plot.obs16a <- ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_test +  # use the common scale here
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 28),
    legend.text  = element_text(size = 24),
    legend.margin = margin(400, 40, 0, 40)) +
  coord_equal()

plot.pde16 <- ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_test +  # use the common scale
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

plot.gam16 <- ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.INLA16 <- ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.rt16 <- ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = RT.pruned.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.obs18 <- ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_test +  # use the common scale
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

plot.pde18 <- ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_test +  # use the common scale
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

plot.gam18 <- ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.INLA18 <- ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.rt18 <- ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = RT.pruned.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.obs20 <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_test +  # use the common scale
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

plot.pde20 <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_test +  # use the common scale
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

plot.gam20 <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.INLA20 <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

plot.rt20 <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = RT.pruned.pred)) +
  geom_tile() +
  common_fill_training +  # use the common scale
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

combined_plot_obs_PDE_test <- plot.obs16 + plot.pde16 + 
  plot.obs18 + plot.pde18 + plot.obs20 + plot.pde20 + 
  plot_layout(nrow = 3, ncol = 2, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

# png(filename = file.path(plot_dir, "obs_PDE_test_no_annotation.png"), height = 6300, width = 6000, res = 300)
# combined_plot_obs_PDE_test
# dev.off()

# column labels as tiny ggplots
col_labels_test <- c(
  " ",
  "Observation", 
  "Prediction with\nPDE-Based Method"
)

col_label_ts <- lapply(col_labels_test, function(txt) {
  ggplot() +
    geom_text(aes(x = 0, y = 0, label = txt), size = 16, fontface = "bold") +
    theme_void() + coord_fixed()
})

# Row labels as tiny ggplots
row_labels_test <- c("2016", "2018", "2020")
row_label_ts <- lapply(row_labels_test, function(txt) {
  ggplot() +
    geom_text(aes(x = 0, y = 0, label = txt), size = 16, fontface = "bold", hjust = 1) +
    theme_void() + coord_fixed()
})

axis_theme_test <- theme(axis.text = element_text(size = 24),
                         plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

plot.obs16a.lb <- plot.obs16a + axis_theme_test 
plot.pde16.lb <- plot.pde16 + axis_theme_test 

plot.obs18.lb <- plot.obs18 + axis_theme_test 
plot.pde18.lb <- plot.pde18 + axis_theme_test 

plot.obs20.lb <- plot.obs20 + axis_theme_test 
plot.pde20.lb <- plot.pde20 + axis_theme_test 

combined_test_with_labels <- col_label_ts[[1]] + col_label_ts[[2]] +
  col_label_ts[[3]] + 
  row_label_ts[[1]] + plot.obs16a.lb + plot.pde16.lb + 
  row_label_ts[[2]] + plot.obs18.lb + plot.pde18.lb + 
  row_label_ts[[3]] + plot.obs20.lb + plot.pde20.lb + 
  plot_layout(nrow = 4, ncol = 3, guides = "collect")

# png(filename = file.path(plot_dir, "obs_PDE_test_with_labels.png"), height = 8300, width = 7000, res = 300)
# combined_test_with_labels
# dev.off()

title.plot.ts <- theme(plot.title = element_text(size = 28, 
                                                 color = "black", 
                                                 face = "bold", 
                                                 hjust = 0.5, vjust = 1))

plot.obs16.tl <- plot.obs16 + labs(title = "Observation\n(2016)") + title.plot.ts
plot.pde16.tl <- plot.pde16 + labs(title = "Prediction with\nPDE-Based Method (2016)") + title.plot.ts
plot.gam16.tl <- plot.gam16 + labs(title = "Prediction with\nGAM (2016)") + title.plot.ts
plot.INLA16.tl <- plot.INLA16 + labs(title = "Prediction with\nINLA (2016)") + title.plot.ts
plot.rt16.tl <- plot.rt16 + labs(title = "Prediction with\nRT (2016)") + title.plot.ts

plot.obs18.tl <- plot.obs18 + labs(title = "Observation\n(2018)") + title.plot.ts
plot.pde18.tl <- plot.pde18 + labs(title = "Prediction with\nPDE-Based Method (2018)") + title.plot.ts
plot.gam18.tl <- plot.gam18 + labs(title = "Prediction with\nGAM (2018)") + title.plot.ts
plot.INLA18.tl <- plot.INLA18 + labs(title = "Prediction with\nINLA (2018)") + title.plot.ts
plot.rt18.tl <- plot.rt18 + labs(title = "Prediction with\nRT (2018)") + title.plot.ts

plot.obs20.tl <- plot.obs20 + labs(title = "Observation\n(2020)") + title.plot.ts
plot.pde20.tl <- plot.pde20 + labs(title = "Prediction with\nPDE-Based Method (2020)") + title.plot.ts
plot.gam20.tl <- plot.gam20 + labs(title = "Prediction with\nGAM (2020)") + title.plot.ts
plot.INLA20.tl <- plot.INLA20 + labs(title = "Prediction with\nINLA (2020)") + title.plot.ts
plot.rt20.tl <- plot.rt20 + labs(title = "Prediction with\nRT (2020)") + title.plot.ts

combined_plot_obs_PDE_test_tl <- plot.obs16.tl + plot.pde16.tl + 
  plot.obs18.tl + plot.pde18.tl + 
  plot.obs20.tl + plot.pde20.tl + 
  plot_layout(nrow = 3, ncol = 2, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

png(filename = file.path(plot_dir, "obs_PDE_test_inside_title.png"), height = 9000, width = 5000, res = 300)
combined_plot_obs_PDE_test_tl
dev.off()

combined_plot_obs_PDE_GAM_RT_test_tl <- plot.obs16.tl + plot.pde16.tl + plot.gam16.tl + plot.INLA16.tl + plot.rt16.tl + 
  plot.obs18.tl + plot.pde18.tl + plot.gam18.tl + plot.INLA18.tl + plot.rt18.tl + 
  plot.obs20.tl + plot.pde20.tl + plot.gam20.tl + plot.INLA20.tl + plot.rt20.tl +
  plot_layout(nrow = 3, ncol = 5, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

png(filename = file.path(plot_dir, "obs_PDE_GAM_RT_test_inside_title.png"), height = 9000, width = 13500, res = 300)
combined_plot_obs_PDE_GAM_RT_test_tl
dev.off()

# for slides
common_fill_sl <- scale_fill_gradient(
  low="white", high="green4", limits=pic_limit,
  name = "tree\ncoverage",
  breaks = seq(0, 1, by = 0.2),
  labels = seq(0, 1, by = 0.2),
  guide = guide_colorbar(
    barheight = unit(30, "cm"),
    barwidth  = unit(2, "cm"),
    ticks = TRUE, 
    ticks.colour = "black",
    title.position = "top",
    title.hjust = 0.5,
    title.vjust = 3
  )
)

title.plot.sl <- theme(plot.title = element_text(size = 28, 
                                                 color = "black", 
                                                 face = "bold", 
                                                 hjust = 0.5, vjust = 1))

plot.obs10.sl <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_sl +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 28),
    legend.text  = element_text(size = 24),
    legend.margin = margin(20, 40, 20, 40)
  ) +
  coord_equal() +
  labs(title = "Observation\n(2010)") + 
  title.plot.sl

plot.pde10.sl <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_sl +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Estimates with\nPDE-Based Method (2010)") +
  title.plot.sl

plot.gam10.sl <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_sl +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Estimates with\nGAM (2010)") + 
  title.plot.sl

plot.INLA10.sl <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_sl +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Estimates with\nINLA (2010)") + 
  title.plot.sl

plot.rt10.sl <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = RT.original.pred)) +
  geom_tile() +
  common_fill_sl +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() + 
  labs(title = "Estimates with\nRT (2010)") + 
  title.plot.sl

plot.obs20.sl <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_sl +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) + 
  coord_equal() +
  labs(title = "Observation\n(2020)") +
  title.plot.sl

plot.pde20.sl <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_sl +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) + 
  coord_equal() + 
  labs(title = "Prediction with\nPDE-Based Method (2020)") +
  title.plot.sl

plot.gam20.sl <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_test +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Prediction with\nGAM (2020)") + 
  title.plot.sl

plot.INLA20.sl <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_test +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Prediction with\nINLA (2020)") + 
  title.plot.sl

plot.rt20.sl <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = RT.original.pred)) +
  geom_tile() +
  common_fill_sl +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Prediction with\nRT (2020)") + 
  title.plot.sl

combined_plot_obs_PDE_GAM_INLA_RT_test_sl <- plot.obs10.sl + plot.pde10.sl + plot.gam10.sl + plot.INLA10.sl + plot.rt10.sl +
  plot.obs20.sl + plot.pde20.sl + plot.gam20.sl + plot.INLA20.sl + plot.rt20.sl +
  plot_layout(nrow = 2, ncol = 5, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

png(filename = file.path(plot_dir, "obs_PDE_GAM_INLA_RT_inside_slides.png"), height = 6200, width = 13500, res = 300)
combined_plot_obs_PDE_GAM_INLA_RT_test_sl
dev.off()

# 2+2
common_fill_dist <- scale_fill_gradient(
  low="white", high="green4", limits=pic_limit,
  name = "tree cover\n(proportion)",
  breaks = seq(0, 1, by = 0.2),
  labels = seq(0, 1, by = 0.2),
  guide = guide_colorbar(
    barheight = unit(60, "cm"),
    barwidth  = unit(2, "cm"),
    ticks = TRUE, 
    ticks.colour = "black",
    title.position = "top",
    title.hjust = 0.5,
    title.vjust = 3
  )
)

title.plot.dist <- theme(plot.title = element_text(size = 32, 
                                                   color = "black", 
                                                   face = "bold", 
                                                   hjust = 0.5, vjust = 1))

plot.obs00.dist <- ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 30),
    legend.text  = element_text(size = 28),
    legend.margin = margin(20, 40, 20, 40)
  ) +
  coord_equal() +
  labs(title = "Observed\n(2000)") + 
  title.plot.dist

plot.pde00.dist <- ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Fitted, PDE-based model\n(2000)") +
  title.plot.dist

plot.gam00.dist <- ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Fitted, GAM\n(2000)") + 
  title.plot.dist

plot.INLA00.dist <- ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Fitted, INLA\n(2000)") + 
  title.plot.dist

plot.rt00.dist <- ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = RT.original.pred)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() + 
  labs(title = "Fitted, RT\n(2000)") + 
  title.plot.dist

plot.obs10.dist <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) + 
  coord_equal() +
  labs(title = "Observed\n(2010)") +
  title.plot.dist

plot.pde10.dist <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) + 
  coord_equal() + 
  labs(title = "Fitted, PDE-based model\n(2010)") +
  title.plot.dist

plot.gam10.dist <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Fitted, GAM\n(2010)") + 
  title.plot.dist

plot.INLA10.dist <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Fitted, INLA\n(2010)") + 
  title.plot.dist

plot.rt10.dist <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = RT.original.pred)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Fitted, RT\n(2010)") + 
  title.plot.dist

plot.obs16.dist <- ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) + 
  coord_equal() +
  labs(title = "Observed\n(2016)") +
  title.plot.dist

plot.pde16.dist <- ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) + 
  coord_equal() + 
  labs(title = "Predicted, PDE-based model\n(2016)") +
  title.plot.dist

plot.gam16.dist <- ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Predicted, GAM\n(2016)") + 
  title.plot.dist

plot.INLA16.dist <- ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Predicted, INLA\n(2016)") + 
  title.plot.dist

plot.rt16.dist <- ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = RT.original.pred)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Predicted, RT\n(2016)") + 
  title.plot.dist

plot.obs20.dist <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) + 
  coord_equal() +
  labs(title = "Observed\n(2020)") +
  title.plot.dist

plot.pde20.dist <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = pred.expectation)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) + 
  coord_equal() + 
  labs(title = "Predicted, PDE-based model\n(2020)") +
  title.plot.dist

plot.gam20.dist <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = GAM.pred)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Predicted, GAM\n(2020)") + 
  title.plot.dist

plot.INLA20.dist <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = INLA.pred)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Predicted, INLA\n(2020)") + 
  title.plot.dist

plot.rt20.dist <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = RT.original.pred)) +
  geom_tile() +
  common_fill_dist +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  ) +
  coord_equal() +
  labs(title = "Predicted, RT\n(2020)") + 
  title.plot.dist

combined_plot_obs_PDE_GAM_INLA_RT_4yrs <- plot.obs00.dist + plot.pde00.dist + plot.gam00.dist + plot.INLA00.dist + plot.rt00.dist +
  plot.obs10.dist + plot.pde10.dist + plot.gam10.dist + plot.INLA10.dist + plot.rt10.dist +
  plot.obs16.dist + plot.pde16.dist + plot.gam16.dist + plot.INLA16.dist + plot.rt16.dist +
  plot.obs20.dist + plot.pde20.dist + plot.gam20.dist + plot.INLA20.dist + plot.rt20.dist +
  plot_layout(nrow = 4, ncol = 5, guides = "collect") & 
  scale_x_continuous(breaks = c(0, 3000, 6000, 8700),
                     labels = c("0", "3000", "6000", "x (m)"),
                     expand = expansion(add = c(300, 300))) &
  scale_y_continuous(breaks = c(0, 3000, 6000, 9000),
                     labels = c("0", "3000", "6000", "y (m)"),
                     expand = expansion(add = c(300, 300))) &
  theme(axis.text = element_text(size = 28),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 28),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

png(filename = file.path(plot_dir, "obs_PDE_GAM_INLA_RT_4yrs.png"), height = 10800, width = 13500, res = 300)
combined_plot_obs_PDE_GAM_INLA_RT_4yrs
dev.off()

# mu
common_fill_mu <- scale_fill_gradient(
  low="white", high="green4", limits=pic_limit,
  name = "tree cover\n(proportion)",
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

plot.estmu85 <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  common_fill_mu +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 30),
    legend.text  = element_text(size = 28),
    legend.margin = margin(20, 40, 20, 40))  + 
  coord_equal()

plot.estmu85a <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = )) +
  geom_tile() +
  common_fill_mu +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 30),
    legend.text  = element_text(size = 28),
    legend.margin = margin(20, 40, 20, 40))  + 
  coord_equal()

plot.estmu90 <- ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  common_fill_mu +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  )  + 
  coord_equal()

plot.estmu95 <- ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  common_fill_mu +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  )  + 
  coord_equal()

plot.estmu00 <- ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  common_fill_mu +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  )  + 
  coord_equal()

plot.estmu05 <- ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  common_fill_mu +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  )  + 
  coord_equal()

plot.estmu08 <- ggplot(data = dat08, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  common_fill_mu +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  )  + 
  coord_equal()

plot.estmu10 <- ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  common_fill_mu +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  )  + 
  coord_equal()

plot.estmu15 <- ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  common_fill_mu +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  )  + 
  coord_equal()

plot.estmu16 <- ggplot(data = dat16, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  common_fill_mu +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  )  + 
  coord_equal()

plot.estmu18 <- ggplot(data = dat18, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  common_fill_mu +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  )  + 
  coord_equal()

plot.estmu20 <- ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = mu.pred)) +
  geom_tile() +
  common_fill_mu +  # use the common scale
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # hide the legend for non-first subplots
  )  + 
  coord_equal()

title.plot.mu <- theme(plot.title = element_text(size = 32, 
                                                 color = "black", 
                                                 # face = "bold", 
                                                 hjust = 0.5, vjust = 1))

plot.estmu85.tl <- plot.estmu85 + labs(title = expression(bold(mu * ", 1985 (initial condition)"))) + title.plot.mu
plot.estmu90.tl <- plot.estmu90 + labs(title = expression(bold(paste("Fitted " * mu * " (1990)")))) + title.plot.mu
plot.estmu95.tl <- plot.estmu95 + labs(title = expression(bold(paste("Fitted " * mu * " (1995)")))) + title.plot.mu
plot.estmu00.tl <- plot.estmu00 + labs(title = expression(bold(paste("Fitted " * mu * " (2000)")))) + title.plot.mu
plot.estmu05.tl <- plot.estmu05 + labs(title = expression(bold(paste("Fitted " * mu * " (2005)")))) + title.plot.mu
plot.estmu08.tl <- plot.estmu08 + labs(title = expression(bold(paste("Fitted " * mu * " (2008)")))) + title.plot.mu
plot.estmu10.tl <- plot.estmu10 + labs(title = expression(bold(paste("Fitted " * mu * " (2010)")))) + title.plot.mu
plot.estmu15.tl <- plot.estmu15 + labs(title = expression(bold(paste("Fitted " * mu * " (2015)")))) + title.plot.mu
plot.estmu16.tl <- plot.estmu16 + labs(title = expression(bold(paste("Predicted " * mu * " (2016)")))) + title.plot.mu
plot.estmu18.tl <- plot.estmu18 + labs(title = expression(bold(paste("Predicted " * mu * " (2018)")))) + title.plot.mu
plot.estmu20.tl <- plot.estmu20 + labs(title = expression(bold(paste("Predicted " * mu * " (2020)")))) + title.plot.mu

combined_plot_mu_tl <- plot.estmu85.tl + plot.estmu90.tl + plot.estmu95.tl + 
  plot.estmu00.tl + plot.estmu08.tl + plot.estmu15.tl + 
  plot.estmu16.tl + plot.estmu18.tl + plot.estmu20.tl + 
  plot_layout(nrow = 3, ncol = 3, guides = "collect") & 
  scale_x_continuous(breaks = c(0, 3000, 6000, 8650),
                     labels = c("0", "3000", "6000", "x (m)"),
                     expand = expansion(add = c(300, 300))) &
  scale_y_continuous(breaks = c(0, 3000, 6000, 9000),
                     labels = c("0", "3000", "6000", " y (m)"),
                     expand = expansion(add = c(300, 300))) &
  theme(axis.text = element_text(size = 28),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 28),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

png(filename = file.path(plot_dir, "mu_inside_title.png"), height = 7000, width = 7500, res = 300)
combined_plot_mu_tl
dev.off()

combined_plot_mu_sl <- plot.estmu85.tl + plot.estmu90.tl + 
  plot.estmu00.tl + plot.estmu10.tl + plot.estmu15.tl + 
  plot.estmu20.tl + 
  plot_layout(nrow = 2, ncol = 3, guides = "collect") & 
  theme(axis.text = element_text(size = 24),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

png(filename = file.path(plot_dir, "mu_inside_slides.png"), height = 5000, width = 7500, res = 300)
combined_plot_mu_sl
dev.off()

# diffusion coefs
title.plot.coef <- theme(plot.title = element_text(size = 20, 
                                                   color = "black", 
                                                   face = "bold", 
                                                   hjust = 0.5, vjust = 1))

plot.estdiffu <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = diffu.prediction)) +
  geom_tile() +
  scale_fill_gradient(
    low="white", high="blue", limit = c(0,2000), 
    name = expression(m^2/yr),             
    breaks = seq(0, 2000, by = 500),
    labels = seq(0, 2000, by = 500),
    guide = guide_colorbar(
      barheight = unit(13.5, "cm"),    
      barwidth  = unit(1, "cm"),
      ticks = TRUE, 
      ticks.colour = "black",
      title.position = "top",
      title.hjust = 0.5,
      title.vjust = 2,
      label = TRUE
      )
  ) +
  theme_bw() +
  theme(plot.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14)) +
  coord_equal() 

png(filename = file.path(plot_dir, "diffusion_coefficients_no_title.png"), height = 2000, width = 2400, res = 300)
plot.estdiffu
dev.off()

plot.estdiffu.tl <- plot.estdiffu + labs(title = "Estimated Diffusion Coefficients") + title.plot.coef +
  scale_x_continuous(breaks = c(0, 3000, 6000, 9000),
                     labels = c("0", "3000", "6000", "dist. (m)"),
                     expand = expansion(add = c(200, 200))) +
  scale_y_continuous(breaks = c(0, 3000, 6000, 9000),
                     labels = c("0", "3000", "6000", "dist.\n(m)"),
                     expand = expansion(add = c(200, 200))) +
  theme(axis.text.x = element_text(margin = margin(t = 10)),
        axis.text.y = element_text(margin = margin(r = 10)),
        axis.ticks = element_blank())

png(filename = file.path(plot_dir, "diffusion_coefficients.png"), height = 2200, width = 2400, res = 300)
plot.estdiffu.tl
dev.off()

# growth coefficients
plot.estgrowth <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = growth.prediction)) +
  geom_tile() +
  scale_fill_gradient(
    low="white", high="blue", limit = c(0, 0.08), 
    name = "1/yr",
    breaks = seq(0, 0.08, by = 0.02),
    labels = seq(0, 0.08, by = 0.02),           
    guide = guide_colorbar(
      barheight = unit(13.5, "cm"),    
      barwidth  = unit(1, "cm"),
      ticks = TRUE, 
      ticks.colour = "black",
      title.position = "top",
      title.hjust = 0.5,
      title.vjust = 2,
      label = TRUE
      )
  ) +
  theme_bw() +
  theme(plot.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14)) +
  coord_equal() 

png(filename = file.path(plot_dir, "growth_coefficients_no_title.png"), height = 2000, width = 2400, res = 300)
plot.estgrowth
dev.off()

plot.estgrowth.tl <- plot.estgrowth + labs(title = "Estimated Growth Coefficients") + title.plot.coef +
  scale_x_continuous(breaks = c(0, 3000, 6000, 9000),
                     labels = c("0", "3000", "6000", "dist. (m)"),
                     expand = expansion(add = c(200, 200))) +
  scale_y_continuous(breaks = c(0, 3000, 6000, 9000),
                     labels = c("0", "3000", "6000", "dist.\n(m)"),
                     expand = expansion(add = c(200, 200))) +
  theme(axis.text.x = element_text(margin = margin(t = 10)),
        axis.text.y = element_text(margin = margin(r = 10)),
        axis.ticks = element_blank())

png(filename = file.path(plot_dir, "growth_coefficients.png"), height = 2200, width = 2400, res = 300)
plot.estgrowth.tl
dev.off()

# gamma
plot.estgamma <- ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = gamma.prediction)) +
  geom_tile() +
  scale_fill_gradient(
    low="white", high="blue", limit = c(0, 1), 
    name = "probability",             
    breaks = seq(0, 1, by = 0.2),
    labels = seq(0, 1, by = 0.2),
    guide = guide_colorbar(
      barheight = unit(13.5, "cm"),    
      barwidth  = unit(1, "cm"),
      ticks = TRUE, 
      ticks.colour = "black",
      title.position = "top",
      title.hjust = 0.5,
      title.vjust = 2,
      label = TRUE
      )
  ) +
  theme_bw() +
  theme(plot.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14)) +
  coord_equal() 

png(filename = file.path(plot_dir, "gamma_no_title.png"), height = 2000, width = 2400, res = 300)
plot.estgamma
dev.off()

plot.estgamma.tl <- plot.estgamma + 
  labs(title = expression(bold(paste("Estimated ",gamma)))) +
  theme(plot.title = element_text(size = 20, 
                                  color = "black", 
                                  # face = "bold", 
                                  hjust = 0.5, vjust = 1)) +
  scale_x_continuous(breaks = c(0, 3000, 6000, 9000),
                     labels = c("0", "3000", "6000", "dist. (m)"),
                     expand = expansion(add = c(200, 200))) +
  scale_y_continuous(breaks = c(0, 3000, 6000, 9000),
                     labels = c("0", "3000", "6000", "dist.\n(m)"),
                     expand = expansion(add = c(200, 200))) +
  theme(axis.text.x = element_text(margin = margin(t = 10)),
        axis.text.y = element_text(margin = margin(r = 10)),
        axis.ticks = element_blank())

png(filename = file.path(plot_dir, "gamma.png"), height = 2200, width = 2400, res = 300)
plot.estgamma.tl
dev.off()

# histrogram
histdat.train <- dat %>% filter(type == "train" & edge == 0 & Stream == 0)
histdat.test <- dat %>% filter(type == "test" & edge == 0 & Stream == 0)
histdat <- dat %>% filter(edge == 0 & Stream == 0)

png(filename = file.path(plot_dir, "histogram_mu.png"), height = 2500, width = 2000, res = 300)
par(mfrow = c(1, 1), 
    cex.lab = 2,
    cex.axis = 1.8,
    cex.main = 2,
    mar = c(5, 5, 10, 5),
    oma = c(1, 1, 1, 1))
hist(histdat[ ,"mu.pred"], main = expression(bold(paste("Distribution of Estimated " * mu))), 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2500000), 
     breaks = 20)
dev.off()

png(filename = file.path(plot_dir, "histogram_training_no_title.png"), height = 2500, width = 4000, res = 300)
par(mfrow = c(1, 2), 
    cex.lab = 2,
    cex.axis = 1.8, 
    mar = c(5, 5, 2, 2),
    oma = c(1, 1, 1, 1))
hist(histdat.train[ ,"trees"], main = "", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2000000),
     breaks = 20)
hist(histdat.train[ ,"pred.expectation"], main = "", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2000000), 
     breaks = 20)
dev.off()

png(filename = file.path(plot_dir, "histogram_training.png"), height = 2500, width = 4000, res = 300)
par(mfrow = c(1, 2), 
    cex.lab = 2,
    cex.axis = 1.8,
    cex.main = 2.5,
    mar = c(5, 5, 10, 2),
    oma = c(1, 1, 1, 1))
hist(histdat.train[ ,"trees"], main = "Distribution of Tree Coverage\n(Training Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2000000), 
     breaks = 20)
hist(histdat.train[ ,"pred.expectation"], main = "Distribution of Estimated Tree Coverage\nwith PDE-Based Method\n(Training Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2000000), 
     breaks = 20)
dev.off()

png(filename = file.path(plot_dir, "histogram_obs_PDE_GAM_RT_training.png"), height = 2500, width = 8000, res = 300)
par(mfrow = c(1, 4), 
    cex.lab = 2,
    cex.axis = 1.8,
    cex.main = 2.5,
    mar = c(5, 5, 10, 2),
    oma = c(1, 1, 1, 1))
hist(histdat.train[ ,"trees"], main = "Distribution of Tree Coverage\n(Training Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2000000), 
     breaks = 20)
hist(histdat.train[ ,"pred.expectation"], main = "Distribution of Estimated Tree Coverage\nwith PDE-Based Method\n(Training Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2000000), 
     breaks = 20)
hist(histdat.train[ ,"GAM.pred"], main = "Distribution of Estimated Tree Coverage\nwith GAM\n(Training Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2000000), 
     breaks = 10)
hist(histdat.train[ ,"RT.pruned.pred"], main = "Distribution of Estimated Tree Coverage\nwith RT\n(Training Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2000000), 
     breaks = 20)
dev.off()

png(filename = file.path(plot_dir, "histogram_test_no_title.png"), height = 2500, width = 4000, res = 300)
par(mfrow = c(1, 2), 
    cex.lab = 2,
    cex.axis = 1.8, 
    mar = c(5, 5, 2, 2),
    oma = c(1, 1, 1, 1))
hist(histdat.test[ ,"trees"], main = "", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 300000), 
     breaks = 20)
hist(histdat.test[ ,"pred.expectation"], main = "", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 300000), 
     breaks = 20)
dev.off()

png(filename = file.path(plot_dir, "histogram_test.png"), height = 2500, width = 4000, res = 300)
par(mfrow = c(1, 2), 
    cex.lab = 2,
    cex.axis = 1.8,
    cex.main = 2,
    mar = c(5, 5, 10, 2),
    oma = c(1, 1, 1, 1))
hist(histdat.test[ ,"trees"], main = "Distribution of Tree Coverage\n(Test Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 300000), 
     breaks = 20)
hist(histdat.test[ ,"pred.expectation"], main = "Distribution of Predicted Tree Coverage\nwith PDE-Based Method\n(Test Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 300000), 
     breaks = 20)
dev.off()

png(filename = file.path(plot_dir, "histogram_obs_PDE_GAM_RT_test.png"), height = 2500, width = 8000, res = 300)
par(mfrow = c(1, 4), 
    cex.lab = 2,
    cex.axis = 1.8,
    cex.main = 2,
    mar = c(5, 5, 10, 2),
    oma = c(1, 1, 1, 1))
hist(histdat.test[ ,"trees"], main = "Distribution of Tree Coverage\n(Test Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 300000), 
     breaks = 20)
hist(histdat.test[ ,"pred.expectation"], main = "Distribution of Predicted Tree Coverage\nwith PDE-Based Method\n(Test Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 300000), 
     breaks = 20)
hist(histdat.test[ ,"GAM.pred"], main = "Distribution of Predicted Tree Coverage\nwith GAM\n(Test Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 300000), 
     breaks = 10)
hist(histdat.test[ ,"RT.pruned.pred"], main = "Distribution of Predicted Tree Coverage\nwith RT\n(Test Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 300000), 
     breaks = 20)
dev.off()

png(filename = file.path(plot_dir, "histogram_obs_PDE_GAM_INLA_RT.png"), height = 5000, width = 10000, res = 300)
par(mfrow = c(2, 5), 
    cex.lab = 2,
    cex.axis = 1.8,
    cex.main = 2.5,
    mar = c(10, 5, 10, 2),
    oma = c(1, 1, 1, 1))
hist(histdat.train[ ,"trees"], main = "Distribution of Tree Coverage\n(Training Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2000000), 
     breaks = 20)
hist(histdat.train[ ,"pred.expectation"], main = "Distribution of Estimated Tree Coverage\nwith PDE-Based Method\n(Training Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2000000), 
     breaks = 20)
hist(histdat.train[ ,"GAM.pred"], main = "Distribution of Estimated Tree Coverage\nwith GAM\n(Training Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2000000), 
     breaks = 10)
hist(histdat.train[ ,"INLA.pred"], main = "Distribution of Estimated Tree Coverage\nwith INLA\n(Training Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2000000), 
     breaks = 20)
hist(histdat.train[ ,"RT.pruned.pred"], main = "Distribution of Estimated Tree Coverage\nwith RT\n(Training Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 2000000), 
     breaks = 20)
hist(histdat.test[ ,"trees"], main = "Distribution of Tree Coverage\n(Test Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 300000), 
     breaks = 20)
hist(histdat.test[ ,"pred.expectation"], main = "Distribution of Predicted Tree Coverage\nwith PDE-Based Method\n(Test Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 300000), 
     breaks = 20)
hist(histdat.test[ ,"GAM.pred"], main = "Distribution of Predicted Tree Coverage\nwith GAM\n(Test Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 300000), 
     breaks = 10)
hist(histdat.test[ ,"INLA.pred"], main = "Distribution of Predicted Tree Coverage\nwith INLA\n(Test Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 300000), 
     breaks = 20)
hist(histdat.test[ ,"RT.pruned.pred"], main = "Distribution of Predicted Tree Coverage\nwith RT\n(Test Data)", 
     xlab = "tree coverage", ylab = "frequency",
     xlim = c(0,1), ylim = c(0, 300000), 
     breaks = 20)
dev.off()

# MSE
dat.without.edge.train <- dat %>% filter(type == "train" & edge == 0 & Stream == 0)
dat.without.edge.test <- dat %>% filter(type == "test" & edge == 0 & Stream == 0)

MSE.train.PDE <- mean((dat.without.edge.train$trees-dat.without.edge.train$pred.expectation)^2, na.rm = TRUE)
MSE.train.GAM <- mean((dat.without.edge.train$trees-dat.without.edge.train$GAM.pred)^2, na.rm = TRUE)
MSE.train.INLA <- mean((dat.without.edge.train$trees-dat.without.edge.train$INLA.pred)^2, na.rm = TRUE)
MSE.train.pruned.rt <- mean((dat.without.edge.train$trees-dat.without.edge.train$RT.pruned.pred)^2, na.rm = TRUE)

MSE.test.PDE <- mean((dat.without.edge.test$trees-dat.without.edge.test$pred.expectation)^2, na.rm = TRUE)
MSE.test.GAM <- mean((dat.without.edge.test$trees-dat.without.edge.test$GAM.pred)^2, na.rm = TRUE)
MSE.test.INLA <- mean((dat.without.edge.test$trees-dat.without.edge.test$INLA.pred)^2, na.rm = TRUE)
MSE.test.pruned.rt <- mean((dat.without.edge.test$trees-dat.without.edge.test$RT.pruned.pred)^2, na.rm = TRUE)

sink(file = file.path(uncertainty_dir, "MSE.txt"))
cat("training: year = 1985-2015; test: time = 2016-2020.\n")
cat("MSE PDE training data: ")
MSE.train.PDE
cat("MSE GAM training data: ")
MSE.train.GAM
cat("MSE INLA training data: ")
MSE.train.INLA
cat("MSE RT training data: ")
MSE.train.pruned.rt
cat("MSE PDE test data: ")
MSE.test.PDE
cat("MSE GAM test data: ")
MSE.test.GAM
cat("MSE INLA test data: ")
MSE.test.INLA
cat("MSE RT test data: ")
MSE.test.pruned.rt
sink(file = NULL)

# MSE in terms of homogenization
dat.without.edge.homog.train <- dat %>% filter(type == "train" & edge.homog == 0 & stream.homog == 0)
dat.without.edge.homog.test <- dat %>% filter(type == "test" & edge.homog == 0 & stream.homog == 0)

MSE.train.PDE.homog <- mean((dat.without.edge.homog.train$trees-dat.without.edge.homog.train$pred.expectation)^2, na.rm = TRUE)
MSE.train.GAM.homog <- mean((dat.without.edge.homog.train$trees-dat.without.edge.homog.train$GAM.pred)^2, na.rm = TRUE)
MSE.train.INLA.homog <- mean((dat.without.edge.homog.train$trees-dat.without.edge.homog.train$INLA.pred)^2, na.rm = TRUE)
MSE.train.pruned.rt.homog <- mean((dat.without.edge.homog.train$trees-dat.without.edge.homog.train$RT.pruned.pred)^2, na.rm = TRUE)

MSE.test.PDE.homog <- mean((dat.without.edge.homog.test$trees-dat.without.edge.homog.test$pred.expectation)^2, na.rm = TRUE)
MSE.test.GAM.homog <- mean((dat.without.edge.homog.test$trees-dat.without.edge.homog.test$GAM.pred)^2, na.rm = TRUE)
MSE.test.INLA.homog <- mean((dat.without.edge.homog.test$trees-dat.without.edge.homog.test$INLA.pred)^2, na.rm = TRUE)
MSE.test.pruned.rt.homog <- mean((dat.without.edge.homog.test$trees-dat.without.edge.homog.test$RT.pruned.pred)^2, na.rm = TRUE)

sink(file = file.path(uncertainty_dir, "MSE_homog.txt"))
cat("training: year = 1985-2015; test: time = 2016-2020.\n")
cat("MSE PDE training data: ")
MSE.train.PDE.homog
cat("MSE GAM training data: ")
MSE.train.GAM.homog
cat("MSE INLA training data: ")
MSE.train.INLA.homog
cat("MSE RT training data: ")
MSE.train.pruned.rt.homog
cat("MSE PDE test data: ")
MSE.test.PDE.homog
cat("MSE GAM test data: ")
MSE.test.GAM.homog
cat("MSE INLA test data: ")
MSE.test.INLA.homog
cat("MSE RT test data: ")
MSE.test.pruned.rt.homog
sink(file = NULL)