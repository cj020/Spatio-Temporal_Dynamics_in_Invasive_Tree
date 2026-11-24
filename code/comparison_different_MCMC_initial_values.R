packages <- c("raster", "fields", "RColorBrewer", "rasterVis", "mvtnorm",
              "gridExtra", "fBasics", "coda", "tidyverse", "egg", "paletteer")
missing_packages <- packages[!packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

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
library(paletteer)


# Set the path to the local project directory
working_dir <- "C:/Users/YourName/Documents/Spatio-Temporal_Dynamics_in_Invasive_Tree"

# Set working directory
setwd(working_dir)

# Get the working directory
wd <- getwd()

source(file.path(wd, "code", "functions.R"))

mcmc_dir <- file.path(wd, "real_data", "PDE", "MCMC")
plot_dir <- file.path(mcmc_dir, "MCMC_paths")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

result.table1 <- read.csv(file.path(mcmc_dir, "MCMC_result_20000_real_data_300_by_300_zibinom_homo_a.csv"))
result.table2 <- read.csv(file.path(mcmc_dir, "MCMC_result_20000_real_data_300_by_300_zibinom_homo_b.csv"))
result.table3 <- read.csv(file.path(mcmc_dir, "MCMC_result_20000_real_data_300_by_300_zibinom_homo_c.csv"))
N_iteration <- nrow(result.table1)
START <- 5002
burnin <- burnin_1 <- burnin_2 <- c(START:N_iteration)

# Gelman-Rubin and ESS
beta1.int.mcmc <- mcmc.list(mcmc(result.table1[2], start = START), 
                            mcmc(result.table2[2], start = START),
                            mcmc(result.table3[2], start = START))
r.beta1.int <- gelman.diag(beta1.int.mcmc)$psrf[, "Point est."]
ess.beta1.int <- as.numeric(effectiveSize(beta1.int.mcmc))

beta3.int.mcmc <- mcmc.list(mcmc(result.table1[3], start = START), 
                            mcmc(result.table2[3], start = START),
                            mcmc(result.table3[3], start = START))
r.beta3.int <- gelman.diag(beta3.int.mcmc)$psrf[, "Point est."]
ess.beta3.int <- as.numeric(effectiveSize(beta3.int.mcmc))

betag.int.mcmc <- mcmc.list(mcmc(result.table1[4], start = START), 
                            mcmc(result.table2[4], start = START),
                            mcmc(result.table3[4], start = START))
r.betag.int <- gelman.diag(betag.int.mcmc)$psrf[, "Point est."]
ess.betag.int <- as.numeric(effectiveSize(betag.int.mcmc))

beta1.dist.mcmc <- mcmc.list(mcmc(result.table1[5], start = START), 
                             mcmc(result.table2[5], start = START),
                             mcmc(result.table3[5], start = START))
r.beta1.dist <- gelman.diag(beta1.dist.mcmc)$psrf[, "Point est."]
ess.beta1.dist <- as.numeric(effectiveSize(beta1.dist.mcmc))

beta3.dist.mcmc <- mcmc.list(mcmc(result.table1[6], start = START), 
                             mcmc(result.table2[6], start = START),
                             mcmc(result.table3[6], start = START))
r.beta3.dist <- gelman.diag(beta3.dist.mcmc)$psrf[, "Point est."]
ess.beta3.dist <- as.numeric(effectiveSize(beta3.dist.mcmc))

betag.dist.mcmc <- mcmc.list(mcmc(result.table1[7], start = START), 
                             mcmc(result.table2[7], start = START),
                             mcmc(result.table3[7], start = START))
r.betag.dist <- gelman.diag(betag.dist.mcmc)$psrf[, "Point est."]
ess.betag.dist <- as.numeric(effectiveSize(betag.dist.mcmc))

beta1.elev.mcmc <- mcmc.list(mcmc(result.table1[8], start = START), 
                             mcmc(result.table2[8], start = START),
                             mcmc(result.table3[8], start = START))
r.beta1.elev <- gelman.diag(beta1.elev.mcmc)$psrf[, "Point est."]
ess.beta1.elev <- as.numeric(effectiveSize(beta1.elev.mcmc))

beta3.elev.mcmc <- mcmc.list(mcmc(result.table1[9], start = START), 
                             mcmc(result.table2[9], start = START),
                             mcmc(result.table3[9], start = START))
r.beta3.elev <- gelman.diag(beta3.elev.mcmc)$psrf[, "Point est."]
ess.beta3.elev <- as.numeric(effectiveSize(beta3.elev.mcmc))

betag.elev.mcmc <- mcmc.list(mcmc(result.table1[10], start = START), 
                             mcmc(result.table2[10], start = START),
                             mcmc(result.table3[10], start = START))
r.betag.elev <- gelman.diag(betag.elev.mcmc)$psrf[, "Point est."]
ess.betag.elev <- as.numeric(effectiveSize(betag.elev.mcmc))

sink(file = file.path(plot_dir, "GelmanRubin and ESS.txt"))
cat("beta1.int: ")
c(r.beta1.int, ess.beta1.int)
cat("beta3.int: ")
c(r.beta3.int, ess.beta3.int)
cat("betag.int: ")
c(r.betag.int, ess.betag.int)
cat("beta1.dist: ")
c(r.beta1.dist, ess.beta1.dist)
cat("beta3.dist: ")
c(r.beta3.dist, ess.beta3.dist)
cat("betag.dist: ")
c(r.betag.dist, ess.betag.dist)
cat("beta1.elev: ")
c(r.beta1.elev, ess.beta1.elev)
cat("beta3.elev: ")
c(r.beta3.elev, ess.beta3.elev)
cat("betag.elev: ")
c(r.betag.elev, ess.betag.elev)
sink(file = NULL)

# acceptance rate result 1
# accpt.rate1 <- c(mean(result.table1[1:N_iteration,11]),
#                  mean(result.table1[1:N_iteration,12]),
#                  mean(result.table1[1:N_iteration,13]),
#                  mean(result.table1[1:N_iteration,14]),
#                  mean(result.table1[1:N_iteration,15]),
#                  mean(result.table1[1:N_iteration,16]),
#                  mean(result.table1[1:N_iteration,17]),
#                  mean(result.table1[1:N_iteration,18]),
#                  mean(result.table1[1:N_iteration,19]))

# beta1.int_est1 <- mean(result.table1[burnin,2])
# beta3.int_est1 <- mean(result.table1[burnin,3])
# betag.int_est1 <- mean(result.table1[burnin,4])
# beta1.dist_est1 <- mean(result.table1[burnin,5])
# beta3.dist_est1 <- mean(result.table1[burnin,6])
# betag.dist_est1 <- mean(result.table1[burnin,7])
# beta1.elev_est1 <- mean(result.table1[burnin,8])
# beta3.elev_est1 <- mean(result.table1[burnin,9])
# betag.elev_est1 <- mean(result.table1[burnin,10])

# 95% confidence interval
# beta1.int_sort1 <- sort(result.table1[burnin,2])
# beta1.int_lb1 <- beta1.int_sort1[floor(length(burnin)*0.025)]
# beta1.int_ub1 <- beta1.int_sort1[ceiling(length(burnin)*0.975)]
# 
# beta3.int_sort1 <- sort(result.table1[burnin,3])
# beta3.int_lb1 <- beta3.int_sort1[floor(length(burnin)*0.025)]
# beta3.int_ub1 <- beta3.int_sort1[ceiling(length(burnin)*0.975)]
# 
# betag.int_sort1 <- sort(result.table1[burnin,4])
# betag.int_lb1 <- betag.int_sort1[floor(length(burnin)*0.025)]
# betag.int_ub1 <- betag.int_sort1[ceiling(length(burnin)*0.975)]
# 
# beta1.dist_sort1 <- sort(result.table1[burnin,5])
# beta1.dist_lb1 <- beta1.dist_sort1[floor(length(burnin)*0.025)]
# beta1.dist_ub1 <- beta1.dist_sort1[ceiling(length(burnin)*0.975)]
# 
# beta3.dist_sort1 <- sort(result.table1[burnin,6])
# beta3.dist_lb1 <- beta3.dist_sort1[floor(length(burnin)*0.025)]
# beta3.dist_ub1 <- beta3.dist_sort1[ceiling(length(burnin)*0.975)]
# 
# betag.dist_sort1 <- sort(result.table1[burnin,7])
# betag.dist_lb1 <- betag.dist_sort1[floor(length(burnin)*0.025)]
# betag.dist_ub1 <- betag.dist_sort1[ceiling(length(burnin)*0.975)]
# 
# beta1.elev_sort1 <- sort(result.table1[burnin,8])
# beta1.elev_lb1 <- beta1.elev_sort1[floor(length(burnin)*0.025)]
# beta1.elev_ub1 <- beta1.elev_sort1[ceiling(length(burnin)*0.975)]
# 
# beta3.elev_sort1 <- sort(result.table1[burnin,9])
# beta3.elev_lb1 <- beta3.elev_sort1[floor(length(burnin)*0.025)]
# beta3.elev_ub1 <- beta3.elev_sort1[ceiling(length(burnin)*0.975)]
# 
# betag.elev_sort1 <- sort(result.table1[burnin,10])
# betag.elev_lb1 <- betag.elev_sort1[floor(length(burnin)*0.025)]
# betag.elev_ub1 <- betag.elev_sort1[ceiling(length(burnin)*0.975)]

# 95% confidence interval
# beta1.int_sort2 <- sort(result.table2[burnin,2])
# beta1.int_lb2 <- beta1.int_sort2[floor(length(burnin)*0.025)]
# beta1.int_ub2 <- beta1.int_sort2[ceiling(length(burnin)*0.975)]
# 
# beta3.int_sort2 <- sort(result.table2[burnin,3])
# beta3.int_lb2 <- beta3.int_sort2[floor(length(burnin)*0.025)]
# beta3.int_ub2 <- beta3.int_sort2[ceiling(length(burnin)*0.975)]
# 
# betag.int_sort2 <- sort(result.table2[burnin,4])
# betag.int_lb2 <- betag.int_sort2[floor(length(burnin)*0.025)]
# betag.int_ub2 <- betag.int_sort2[ceiling(length(burnin)*0.975)]
# 
# beta1.dist_sort2 <- sort(result.table2[burnin,5])
# beta1.dist_lb2 <- beta1.dist_sort2[floor(length(burnin)*0.025)]
# beta1.dist_ub2 <- beta1.dist_sort2[ceiling(length(burnin)*0.975)]
# 
# beta3.dist_sort2 <- sort(result.table2[burnin,6])
# beta3.dist_lb2 <- beta3.dist_sort2[floor(length(burnin)*0.025)]
# beta3.dist_ub2 <- beta3.dist_sort2[ceiling(length(burnin)*0.975)]
# 
# betag.dist_sort2 <- sort(result.table2[burnin,7])
# betag.dist_lb2 <- betag.dist_sort2[floor(length(burnin)*0.025)]
# betag.dist_ub2 <- betag.dist_sort2[ceiling(length(burnin)*0.975)]
# 
# beta1.dm_sort2 <- sort(result.table2[burnin,8])
# beta1.dm_lb2 <- beta1.dm_sort2[floor(length(burnin)*0.025)]
# beta1.dm_ub2 <- beta1.dm_sort2[ceiling(length(burnin)*0.975)]
# 
# beta3.dm_sort2 <- sort(result.table2[burnin,9])
# beta3.dm_lb2 <- beta3.dm_sort2[floor(length(burnin)*0.025)]
# beta3.dm_ub2 <- beta3.dm_sort2[ceiling(length(burnin)*0.975)]
# 
# betag.dm_sort2 <- sort(result.table2[burnin,10])
# betag.dm_lb2 <- betag.dm_sort2[floor(length(burnin)*0.025)]
# betag.dm_ub2 <- betag.dm_sort2[ceiling(length(burnin)*0.975)]

# initial1 <- unlist(result.table1[1,2:10])
# initial2 <- unlist(result.table2[1,2:10])
# last1 <- unlist(result.table1[N_iteration,2:10])
# last2 <- unlist(result.table2[N_iteration,2:10])
# 
# est1 <- c(beta1.int_est1, beta3.int_est1, betag.int_est1, 
#           beta1.dist_est1, beta3.dist_est1, betag.dist_est1,
#           beta1.elev_est1, beta3.elev_est1, betag.elev_est1)
# lb1 <- c(beta1.int_lb1, beta3.int_lb1, betag.int_lb1,
#          beta1.dist_lb1, beta3.dist_lb1, betag.dist_lb1,
#          beta1.elev_lb1, beta3.elev_lb1, betag.elev_lb1)
# ub1 <- c(beta1.int_ub1, beta3.int_ub1, betag.int_ub1,
#          beta1.dist_ub1, beta3.dist_ub1, betag.dist_ub1,
#          beta1.elev_ub1, beta3.elev_ub1, betag.elev_ub1)

# est2 <- c(beta1.int_est2, beta3.int_est2, betag.int_est2,
#           beta1.dist_est2, beta3.dist_est2, betag.dist_est2,
#           beta1.dm_est2, beta3.dm_est2, betag.dm_est2)
# lb2 <- c(beta1.int_lb2, beta3.int_lb2, betag.int_lb2,
#          beta1.dist_lb2, beta3.dist_lb2, betag.dist_lb2,
#          beta1.dm_lb2, beta3.dm_lb2, betag.dist_lb2)
# ub2 <- c(beta1.int_ub2, beta3.int_ub2, betag.int_ub2,
#          beta1.dist_ub2, beta3.dist_ub2, betag.dist_ub2,
#          beta1.dm_ub2, beta3.dm_ub2, betag.dm_ub2)

# name <- c("beta1.int", "beta3.int", "betag.int", 
#           "beta1.dist", "beta3.dist", "betag.dist",
#           "beta1.elev", "beta3.elev", "betag.elev")
# 
# table <- data.frame(name = name, ini_1 = initial1, last_1 = last1, est_1 = est1, 
#                     accpt.rate_1 = accpt.rate1, lb_1 = lb1, ub_1 = ub1)

par(mfrow=c(1,1))
# comparison
# png(filename = file.path(plot_dir, "legend.png"))
par(mfrow=c(1,1))
plot(burnin, result.table1[burnin,2], type = 'l',
     main = expression(paste(beta["diffu,int"])),
     xlab = "iterations", ylab = "value", col = "white")
lines(burnin, result.table2[burnin,2], col = "white")
lines(burnin, result.table3[burnin,2], col = "white")
abline(h = -4, col = "white")
legend("topright", legend = c("true value", "est. with ini. value a", "est. with ini. value b", "est. with ini. value c"), col = c("red", "black", "blue","green"),lty=c(1,1,1,1))
# dev.off()

cols <- paletteer::paletteer_d("yarrr::google")

png(filename = file.path(plot_dir, "MCMC_diffu_int.png"))
plot(burnin, result.table1[burnin,2], type = 'l',
     main = expression(paste(beta["diffu,int"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim = c(-5.25, -5.20))
lines(burnin, result.table2[burnin,2], col = cols[3])
lines(burnin, result.table3[burnin,2], col = cols[1])
dev.off()

png(filename = file.path(plot_dir, "MCMC_growth_int.png"))
plot(burnin, result.table1[burnin,3], type = 'l',
     main = expression(paste(beta["growth,int"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim = c(2.785, 2.795))
lines(burnin, result.table2[burnin,3], col = cols[3])
lines(burnin, result.table3[burnin,3], col = cols[1])
dev.off()

png(filename = file.path(plot_dir, "MCMC_gamma_int.png"))
plot(burnin, result.table1[burnin,4], type = 'l',
     main = expression(paste(beta[paste(gamma,",int")])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(-1.84, -1.76))
lines(burnin, result.table2[burnin,4], col = cols[3])
lines(burnin, result.table3[burnin,4], col = cols[1])
dev.off()

png(filename = file.path(plot_dir, "MCMC_diffu_dist.png"))
plot(burnin, result.table1[burnin,5], type = 'l',
     main = expression(paste(beta["diffu,dist"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(-0.75, -0.65))
lines(burnin, result.table2[burnin,5], col = cols[3])
lines(burnin, result.table3[burnin,5], col = cols[1])
dev.off()

png(filename = file.path(plot_dir, "MCMC_growth_dist.png"))
plot(burnin, result.table1[burnin,6], type = 'l',
     main = expression(paste(beta["growth,dist"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(0.41, 0.44))
lines(burnin, result.table2[burnin,6], col = cols[3])
lines(burnin, result.table3[burnin,6], col = cols[1])
dev.off()

png(filename = file.path(plot_dir, "MCMC_gamma_dist.png"))
plot(burnin, result.table1[burnin,7], type = 'l', 
     main = expression(paste(beta[paste(gamma,",dist")])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(-0.85, -0.70))
lines(burnin, result.table2[burnin,7], col = cols[3])
lines(burnin, result.table3[burnin,7], col = cols[1])
dev.off()

png(filename = file.path(plot_dir, "MCMC_diffu_elev.png"))
plot(burnin, result.table1[burnin,8], type = 'l', 
     main = expression(paste(beta["diffu,elev"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(-1.60, -1.58))
lines(burnin, result.table2[burnin,8], col = cols[3])
lines(burnin, result.table3[burnin,8], col = cols[1])
dev.off()

png(filename = file.path(plot_dir, "MCMC_growth_elev.png"))
plot(burnin, result.table1[burnin,9], type = 'l',
     main = expression(paste(beta["growth,elev"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(0.035,0.055))
lines(burnin, result.table2[burnin,9], col = cols[3])
lines(burnin, result.table3[burnin,9], col = cols[1])
dev.off()

png(filename = file.path(plot_dir, "MCMC_gamma_elev.png"))
plot(burnin, result.table1[burnin,10], type = 'l',
     main = expression(paste(beta[paste(gamma,",elev")])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(1.20,1.35))
lines(burnin, result.table2[burnin,10], col = cols[3])
lines(burnin, result.table3[burnin,10], col = cols[1])
abline(h = 1, col = cols[2], lwd = 3)
dev.off()

png(filename = file.path(plot_dir, "MCMC_paths.png"), height = 6000, width = 6000, res = 300)
par(mfrow=c(3,3))
plot(burnin, result.table1[burnin,2], type = 'l',
     main = expression(paste(beta["diffu,int"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim = c(-5.25, -5.20))
lines(burnin, result.table2[burnin,2], col = cols[3])
lines(burnin, result.table3[burnin,2], col = cols[1])

plot(burnin, result.table1[burnin,3], type = 'l',
     main = expression(paste(beta["growth,int"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim = c(2.785, 2.795))
lines(burnin, result.table2[burnin,3], col = cols[3])
lines(burnin, result.table3[burnin,3], col = cols[1])

plot(burnin, result.table1[burnin,4], type = 'l',
     main = expression(paste(beta[paste(gamma,",int")])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(-1.84, -1.76))
lines(burnin, result.table2[burnin,4], col = cols[3])
lines(burnin, result.table3[burnin,4], col = cols[1])

plot(burnin, result.table1[burnin,5], type = 'l',
     main = expression(paste(beta["diffu,dist"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(-0.75, -0.65))
lines(burnin, result.table2[burnin,5], col = cols[3])
lines(burnin, result.table3[burnin,5], col = cols[1])

plot(burnin, result.table1[burnin,6], type = 'l',
     main = expression(paste(beta["growth,dist"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(0.41, 0.44))
lines(burnin, result.table2[burnin,6], col = cols[3])
lines(burnin, result.table3[burnin,6], col = cols[1])

plot(burnin, result.table1[burnin,7], type = 'l', 
     main = expression(paste(beta[paste(gamma,",dist")])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(-0.85, -0.70))
lines(burnin, result.table2[burnin,7], col = cols[3])
lines(burnin, result.table3[burnin,7], col = cols[1])

plot(burnin, result.table1[burnin,8], type = 'l', 
     main = expression(paste(beta["diffu,elev"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(-1.60, -1.58))
lines(burnin, result.table2[burnin,8], col = cols[3])
lines(burnin, result.table3[burnin,8], col = cols[1])

plot(burnin, result.table1[burnin,9], type = 'l',
     main = expression(paste(beta["growth,elev"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(0.035,0.055))
lines(burnin, result.table2[burnin,9], col = cols[3])
lines(burnin, result.table3[burnin,9], col = cols[1])

plot(burnin, result.table1[burnin,10], type = 'l',
     main = expression(paste(beta[paste(gamma,",elev")])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(1.20,1.35))
lines(burnin, result.table2[burnin,10], col = cols[3])
lines(burnin, result.table3[burnin,10], col = cols[1])
dev.off()

png(filename = file.path(plot_dir, "MCMC_paths_one_initial.png"), height = 6000, width = 6000, res = 300)
par(mfrow=c(3,3))
plot(burnin, result.table1[burnin,2], type = 'l',
     main = expression(paste(beta["diffu,int"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = "black", ylim = c(-5.25, -5.20))

plot(burnin, result.table1[burnin,3], type = 'l',
     main = expression(paste(beta["growth,int"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = "black", ylim = c(2.785, 2.795))

plot(burnin, result.table1[burnin,4], type = 'l',
     main = expression(paste(beta[paste(gamma,",int")])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = "black", ylim=c(-1.84, -1.76))

plot(burnin, result.table1[burnin,5], type = 'l',
     main = expression(paste(beta["diffu,dist"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = "black", ylim=c(-0.75, -0.65))

plot(burnin, result.table1[burnin,6], type = 'l',
     main = expression(paste(beta["growth,dist"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = "black", ylim=c(0.41, 0.44))

plot(burnin, result.table1[burnin,7], type = 'l', 
     main = expression(paste(beta[paste(gamma,",dist")])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = "black", ylim=c(-0.85, -0.70))

plot(burnin, result.table1[burnin,8], type = 'l', 
     main = expression(paste(beta["diffu,elev"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = "black", ylim=c(-1.60, -1.58))

plot(burnin, result.table1[burnin,9], type = 'l',
     main = expression(paste(beta["growth,elev"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = "black", ylim=c(0.035,0.055))

plot(burnin, result.table1[burnin,10], type = 'l',
     main = expression(paste(beta[paste(gamma,",elev")])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = "black", ylim=c(1.20,1.35))
dev.off()