# Compare different MCMC initial values

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

setwd("/small_region")
source("functions/raster_and_homogenization_functions.R")

data_dir <- "/small_region/data"
output_dir1 <- "/small_region/output/PDE_1"
output_dir2 <- "/small_region/output/PDE_2"
output_dir3 <- "/small_region/output/PDE_3"
result_dir <- "/small_region/result/PDE"

result.table1 <- read.csv(file.path(output_dir1, "MCMC_result_20000_real_data_300_by_300_zibinom_homo.csv"))
result.table2 <- read.csv(file.path(output_dir2, "MCMC_result_20000_real_data_300_by_300_zibinom_homo.csv"))
result.table3 <- read.csv(file.path(output_dir3, "MCMC_result_20000_real_data_300_by_300_zibinom_homo.csv"))
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

sink(file = file.path(result_dir, "GelmanRubin and ESS.txt"))
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

par(mfrow=c(1,1))
# comparison
# png(filename = "legend.png")
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

png(filename = file.path(result_dir, "MCMC_diffu_int.png"))
plot(burnin, result.table1[burnin,2], type = 'l',
     main = expression(paste(beta["diffu,int"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim = c(-5.25, -5.20))
lines(burnin, result.table2[burnin,2], col = cols[3])
lines(burnin, result.table3[burnin,2], col = cols[1])
dev.off()

png(filename = file.path(result_dir, "MCMC_growth_int.png"))
plot(burnin, result.table1[burnin,3], type = 'l',
     main = expression(paste(beta["growth,int"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim = c(2.785, 2.795))
lines(burnin, result.table2[burnin,3], col = cols[3])
lines(burnin, result.table3[burnin,3], col = cols[1])
dev.off()

png(filename = file.path(result_dir, "MCMC_gamma_int.png"))
plot(burnin, result.table1[burnin,4], type = 'l',
     main = expression(paste(beta[paste(gamma,",int")])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(-1.84, -1.76))
lines(burnin, result.table2[burnin,4], col = cols[3])
lines(burnin, result.table3[burnin,4], col = cols[1])
dev.off()

png(filename = file.path(result_dir, "MCMC_diffu_dist.png"))
plot(burnin, result.table1[burnin,5], type = 'l',
     main = expression(paste(beta["diffu,dist"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(-0.75, -0.65))
lines(burnin, result.table2[burnin,5], col = cols[3])
lines(burnin, result.table3[burnin,5], col = cols[1])
dev.off()

png(filename = file.path(result_dir, "MCMC_growth_dist.png"))
plot(burnin, result.table1[burnin,6], type = 'l',
     main = expression(paste(beta["growth,dist"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(0.41, 0.44))
lines(burnin, result.table2[burnin,6], col = cols[3])
lines(burnin, result.table3[burnin,6], col = cols[1])
dev.off()

png(filename = file.path(result_dir, "MCMC_gamma_dist.png"))
plot(burnin, result.table1[burnin,7], type = 'l', 
     main = expression(paste(beta[paste(gamma,",dist")])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(-0.85, -0.70))
lines(burnin, result.table2[burnin,7], col = cols[3])
lines(burnin, result.table3[burnin,7], col = cols[1])
dev.off()

png(filename = file.path(result_dir, "MCMC_diffu_elev.png"))
plot(burnin, result.table1[burnin,8], type = 'l', 
     main = expression(paste(beta["diffu,elev"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(-1.60, -1.58))
lines(burnin, result.table2[burnin,8], col = cols[3])
lines(burnin, result.table3[burnin,8], col = cols[1])
dev.off()

png(filename = file.path(result_dir, "MCMC_growth_elev.png"))
plot(burnin, result.table1[burnin,9], type = 'l',
     main = expression(paste(beta["growth,elev"])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(0.035,0.055))
lines(burnin, result.table2[burnin,9], col = cols[3])
lines(burnin, result.table3[burnin,9], col = cols[1])
dev.off()

png(filename = file.path(result_dir, "MCMC_gamma_elev.png"))
plot(burnin, result.table1[burnin,10], type = 'l',
     main = expression(paste(beta[paste(gamma,",elev")])), cex.main = 3,
     xlab = "iterations", ylab = "", cex.lab = 2.5,
     cex.axis = 1.75,
     col = cols[4], ylim=c(1.20,1.35))
lines(burnin, result.table2[burnin,10], col = cols[3])
lines(burnin, result.table3[burnin,10], col = cols[1])
abline(h = 1, col = cols[2], lwd = 3)
dev.off()

png(filename = file.path(result_dir, "MCMC_paths.png"), height = 6000, width = 6000, res = 300)
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

png(filename = file.path(result_dir, "MCMC_paths_one_initial.png"), height = 6000, width = 6000, res = 300)
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