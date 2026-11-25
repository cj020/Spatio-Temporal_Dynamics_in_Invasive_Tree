# 12/21/25 INLA
# 1/6/26 several INLA comparison
# 1/11/26 INLA with spatial dependent, dist+elevation+time+s1+s2

# Set the path to the local project directory
working_dir <- "C:/Users/YourName/Documents/Spatio-Temporal_Dynamics_in_Invasive_Tree"

# Set working directory
setwd(working_dir)

# Get the working directory
wd <- getwd()

data_dir <- file.path(wd, "real_data", "data", "N300")
output_dir <- file.path(wd, "real_data", "INLA")
uncertainty_dir <- file.path(wd, "real_data", "prediction_summary", "MSE_and_prediction_uncertainty")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
if (!dir.exists(uncertainty_dir)) {
  dir.create(uncertainty_dir, recursive = TRUE)
}

setwd(output_dir)

required_packages <- c(
  "spdep",
  "tidyverse",
  "ggplot2"
)
missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
if (!requireNamespace("INLA", quietly = TRUE)) {
  install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
}

library(spdep)
library(INLA)
library(tidyverse)
library(ggplot2)

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

# INLA with binomial (iid)
# formula_iid <- response100 ~ distance.str + DEM_30m + time + s1 + s2 
# INLA_fit_iid <- inla(formula_iid, data = dat,
#                      family = "binomial",
#                      Ntrials = size,
#                      control.fixed = list(
#                        mean = 0,
#                        prec = 1 / 10^2,          # precision = 1/variance
#                        mean.intercept = 0,
#                        prec.intercept = 1 / 10^2),
#                      control.predictor = list(compute = TRUE, link = 1),
#                      control.compute = list(dic = TRUE, return.marginals.predictor = TRUE))
# 
# sink(file = file.path(output_dir, "summary_(iid).txt"))
# INLA_fit_iid$summary.fixed
# sink(NULL)

# prediction_iid <- INLA_fit_iid$summary.fitted.values$mean
# lower_CI_iid <- INLA_fit_iid$summary.fitted.values$`0.025quant`
# upper_CI_iid <- INLA_fit_iid$summary.fitted.values$`0.975quant`
# lower_CI_3_iid <- round(lower_CI_iid, 3)
# upper_CI_3_iid <- round(upper_CI_iid, 3)
# 
# dat <- cbind(dat, prediction_iid, lower_CI_3_iid, upper_CI_3_iid)

# INLA with binomial (sp)
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

sink(file = file.path(output_dir, "summary_sp.txt"))
INLA_fit_sp$summary.fixed
sink(NULL)

prediction_sp <- INLA_fit_sp$summary.fitted.values$mean
lower_CI_sp <- INLA_fit_sp$summary.fitted.values$`0.025quant`
upper_CI_sp <- INLA_fit_sp$summary.fitted.values$`0.975quant`
lower_CI_3_sp <- round(lower_CI_sp, 3)
upper_CI_3_sp <- round(upper_CI_sp, 3)

dat <- cbind(dat, prediction_sp, lower_CI_3_sp, upper_CI_3_sp)

# INLA with binomial (sp1)
# formula_sp1 <- response100 ~ distance.str + DEM_30m + time +
# f(space_id, model = "besag", graph="grid_adj.graph")
# INLA_fit_sp1 <- inla(formula_sp1, data = dat,
#                      family = "binomial",
#                      Ntrials = size,
#                      control.fixed = list(
#                        mean = 0,
#                        prec = 1 / 10^2,          # precision = 1/variance
#                        mean.intercept = 0,
#                        prec.intercept = 1 / 10^2),
#                      control.predictor = list(compute = TRUE, link = 1),
#                      control.compute = list(dic = TRUE, return.marginals.predictor = TRUE))
# 
# sink(file = file.path(output_dir, "summary_sp1.txt"))
# INLA_fit_sp1$summary.fixed
# sink(NULL)

# prediction_sp1 <- INLA_fit_sp1$summary.fitted.values$mean
# lower_CI_sp1 <- INLA_fit_sp1$summary.fitted.values$`0.025quant`
# upper_CI_sp1 <- INLA_fit_sp1$summary.fitted.values$`0.975quant`
# lower_CI_3_sp1 <- round(lower_CI_sp1, 3)
# upper_CI_3_sp1 <- round(upper_CI_sp1, 3)
# 
# dat <- cbind(dat, prediction_sp1, lower_CI_3_sp1, upper_CI_3_sp1)

# INLA with binomial (st)
# formula_st <- response100 ~ 1 + distance.str + DEM_30m +
#   f(space_id, model = "besag", graph="grid_adj.graph",
#     group = t_id, control.group = list(model = "ar1"),
#     scale.model = TRUE)
# 
# INLA_fit_st <- inla(formula_st, data = dat,
#                     family = "binomial",
#                     Ntrials = size,
#                     control.fixed = list(
#                       mean = 0,
#                       prec = 1 / 10^2,          # precision = 1/variance
#                       mean.intercept = 0,
#                       prec.intercept = 1 / 10^2),
#                     control.predictor = list(compute = TRUE, link = 1),
#                     control.compute = list(dic = TRUE, waic = TRUE,
#                                            return.marginals.predictor = TRUE))
# 
# sink(file = file.path(output_dir, "summary_st.txt"))
# INLA_fit_st$summary.fixed
# sink(NULL)

# prediction_st <- INLA_fit_st$summary.fitted.values$mean
# lower_CI_st <- INLA_fit_st$summary.fitted.values$`0.025quant`
# upper_CI_st <- INLA_fit_st$summary.fitted.values$`0.975quant`
# lower_CI_3_st <- round(lower_CI_st, 3)
# upper_CI_3_st <- round(upper_CI_st, 3)
# 
# dat <- cbind(dat, prediction_st, lower_CI_3_st, upper_CI_3_st)

# # INLA with zero-inflated binomial
# # 1) Inflation covariates (can be same or different from mean part)
# X_zi <- cbind(1, dat$distance.str, dat$DEM_30m, dat$time, dat$s1, dat$s2)
# Kzi <- ncol(X_zi)
# 
# # priors for inflation coefficients beta1..betaKzi
# hyper_zi <- setNames(
#   lapply(1:Kzi, function(i) list(param = c(0, 1))),
#   paste0("beta", 1:Kzi)
# )
# 
# # 2) ZIB formula: this is for the BINOMIAL MEAN (pi_i)
# formula_zib <- inla.mdata(cbind(dat$response100, dat$size), X_zi) ~
#   1 + distance.str + DEM_30m + time + s1 + s2 +
#   f(space_id, model="besag", graph="grid_adj.graph", scale.model=TRUE)
# 
# # 3) Fit
# INLA_fit_zib <- inla(
#   formula_zib,
#   data = dat,
#   family = "0binomial",
#   control.family = list(
#     link.simple = "logit",   # link for inflation prob p_i
#     hyper = hyper_zi
#   ),
#   control.predictor = list(compute = TRUE),
#   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
#                          return.marginals.predictor = TRUE),
#   control.inla = list(num.threads = 1)
# )
# 
# sink(file = file.path(output_dir, "summary_zib.txt"))
# INLA_fit_zib$summary.fixed
# sink(NULL)
# 
# eta2 <- INLA_fit_zib$summary.linear.predictor$mean
# pi_hat <- plogis(eta2)
# # extract posterior means of inflation coefficients
# b_zi <- INLA_fit_zib$summary.hyperpar$mean[grep("^beta", rownames(INLA_fit_zib$summary.hyperpar))]
# 
# # inflation linear predictor and prob
# eta1 <- as.numeric(X_zi %*% b_zi)
# p_hat <- plogis(eta1)
# 
# # expected count
# mu_hat <- (1 - p_hat) * dat$size * pi_hat
# range(mu_hat)
# prediction_zib <- mu_hat
# dat <- cbind(dat, prediction_zib)

write.csv(dat, file.path(output_dir, "prediction_with_INLA_result_real_data_300_by_300.csv"))

pic_limit <- c(0,1)

dat90 <- dat %>% filter(year == 1990)
png(filename = file.path(output_dir, "observation_in_1990.png"))
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

# png(filename = file.path(output_dir, "estimation_in_1990_iid.png"))
# ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = prediction_iid)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 1990 (iid)") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

png(filename = file.path(output_dir, "estimation_in_1990_sp.png"))
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

# png(filename = file.path(output_dir, "estimation_in_1990_sp1.png"))
# ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = prediction_sp1)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 1990 (sp1)") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = file.path(output_dir, "estimation_in_1990_st.png"))
# ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = prediction_st)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 1990 (st)") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal()
# dev.off()

# png(filename = file.path(output_dir, "estimation_in_1990_with_zib.png"))
# ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = prediction_zib)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 1990 with zib") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal()
# dev.off()

dat00 <- dat %>% filter(year == 2000)
png(filename = file.path(output_dir, "observation_in_2000.png"))
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

# png(filename = file.path(output_dir, "estimation_in_2000_iid.png"))
# ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = prediction_iid)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 2000 (iid)") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

png(filename = file.path(output_dir, "estimation_in_2000_sp.png"))
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

# png(filename = file.path(output_dir, "estimation_in_2000_sp1.png"))
# ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = prediction_sp1)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 2000 (sp1)") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = file.path(output_dir, "estimation_in_2000_st.png"))
# ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = prediction_st)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 2000 (st)") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal()
# dev.off()

# png(filename = file.path(output_dir, "estimation_in_2000_with_zib.png"))
# ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = prediction_zib)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 2000 with zib") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal()
# dev.off()

dat10 <- dat %>% filter(year == 2010)
png(filename = file.path(output_dir, "observation_in_2010.png"))
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

# png(filename = file.path(output_dir, "estimation_in_2010_iid.png"))
# ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = prediction_iid)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 2010 (iid)") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

png(filename = file.path(output_dir, "estimation_in_2010_sp.png"))
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

# png(filename = file.path(output_dir, "estimation_in_2010_sp1.png"))
# ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = prediction_sp1)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 2010 (sp1)") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = file.path(output_dir, "estimation_in_2010_st.png"))
# ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = prediction_st)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 2010 (st)") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal()
# dev.off()

# png(filename = file.path(output_dir, "estimation_in_2010_with_zib.png"))
# ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = prediction_zib)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 2010 with zib") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal()
# dev.off()

dat20 <- dat %>% filter(year == 2020)
png(filename = file.path(output_dir, "observation_in_2020.png"))
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

# png(filename = file.path(output_dir, "estimation_in_2020_iid.png"))
# ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = prediction_iid)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 2020 (iid)") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

png(filename = file.path(output_dir, "estimation_in_2020_sp.png"))
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

# png(filename = file.path(output_dir, "estimation_in_2020_sp1.png"))
# ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = prediction_sp1)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 2020 (sp1)") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() 
# dev.off()

# png(filename = file.path(output_dir, "estimation_in_2020_st.png"))
# ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = prediction_st)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 2020 (st)") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal()
# dev.off()

# png(filename = file.path(output_dir, "estimation_in_2020_with_zib.png"))
# ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = prediction_zib)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
#   labs(title = "Estimated Tree Coverage at year 2020 with zib") +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal()
# dev.off()

# histrogram
# png(filename = file.path(output_dir, "histogram_comparison_training_dataset_obs_vs_est.png"),
#     height = 1000, width = 4000, res = 100)
# par(mfrow=c(1,5),
#     cex.lab = 2,
#     cex.axis = 1.8, 
#     mar = c(5, 5, 2, 2),
#     oma = c(1, 1, 1, 1))
# histdat.train <- dat %>% filter(type == "train" & edge == 0 & Stream == 0)
# hist(histdat.train[ ,"trees"], main = "Histogram of Observations \n(training)")
# hist(histdat.train[ ,"prediction_iid"], main = "Histogram of Estimation (iid)\n(training)")
# hist(histdat.train[ ,"prediction_sp"], main = "Histogram of Estimation (sp)\n(training)")
# hist(histdat.train[ ,"prediction_sp1"], main = "Histogram of Estimation (sp1)\n(training)")
# hist(histdat.train[ ,"prediction_st"], main = "Histogram of Estimation (st)\n(training)")
# # hist(histdat.train[ ,"prediction_zib"], main = "Histogram of Estimation \n(training)")
# dev.off()

# png(filename = file.path(output_dir, "histogram_comparison_test_dataset_obs_vs_est.png"),
#     height = 1000, width = 4000, res = 100)
# par(mfrow=c(1,5),
#     cex.lab = 2,
#     cex.axis = 1.8, 
#     mar = c(5, 5, 2, 2),
#     oma = c(1, 1, 1, 1))
# histdat.test <- dat %>% filter(type == "test" & edge == 0 & Stream == 0)
# hist(histdat.test[ ,"trees"], main = "Histogram of Observations \n(test)")
# hist(histdat.test[ ,"prediction_iid"], main = "Histogram of Prediction (iid)\n(test)")
# hist(histdat.test[ ,"prediction_sp"], main = "Histogram of Prediction (sp)\n(test)")
# hist(histdat.test[ ,"prediction_sp1"], main = "Histogram of Prediction (sp)\n(test)")
# hist(histdat.test[ ,"prediction_st"], main = "Histogram of Prediction (st)\n(test)")
# # hist(histdat.test[ ,"prediction_zib"], main = "Histogram of Prediction \n(test)")
# dev.off()

png(filename = file.path(output_dir, "histogram_comparison_training_dataset_obs_vs_est.png"),
    height = 1000, width = 1600, res = 100)
par(mfrow=c(1,2),
    cex.lab = 2,
    cex.axis = 1.8, 
    mar = c(5, 5, 2, 2),
    oma = c(1, 1, 1, 1))
histdat.train <- dat %>% filter(type == "train" & edge == 0 & Stream == 0)
hist(histdat.train[ ,"trees"], main = "Histogram of Observations \n(training)")
# hist(histdat.train[ ,"prediction_iid"], main = "Histogram of Estimation (iid)\n(training)")
hist(histdat.train[ ,"prediction_sp"], main = "Histogram of Estimation (sp)\n(training)")
# hist(histdat.train[ ,"prediction_sp1"], main = "Histogram of Estimation (sp1)\n(training)")
# hist(histdat.train[ ,"prediction_st"], main = "Histogram of Estimation (st)\n(training)")
# hist(histdat.train[ ,"prediction_zib"], main = "Histogram of Estimation \n(training)")
dev.off()

png(filename = file.path(output_dir, "histogram_comparison_test_dataset_obs_vs_est.png"),
    height = 1000, width = 1600, res = 100)
par(mfrow=c(1,2),
    cex.lab = 2,
    cex.axis = 1.8, 
    mar = c(5, 5, 2, 2),
    oma = c(1, 1, 1, 1))
histdat.test <- dat %>% filter(type == "test" & edge == 0 & Stream == 0)
hist(histdat.test[ ,"trees"], main = "Histogram of Observations \n(test)")
# hist(histdat.test[ ,"prediction_iid"], main = "Histogram of Prediction (iid)\n(test)")
hist(histdat.test[ ,"prediction_sp"], main = "Histogram of Prediction (sp)\n(test)")
# hist(histdat.test[ ,"prediction_sp1"], main = "Histogram of Prediction (sp)\n(test)")
# hist(histdat.test[ ,"prediction_st"], main = "Histogram of Prediction (st)\n(test)")
# hist(histdat.test[ ,"prediction_zib"], main = "Histogram of Prediction \n(test)")
dev.off()

dat.without.edge.train <- dat %>% filter(type == "train" & edge == 0 & Stream == 0)
dat.without.edge.test <- dat %>% filter(type == "test" & edge == 0 & Stream == 0)
# MSE.train.iid <- mean((dat.without.edge.train$trees-dat.without.edge.train$prediction_iid)^2, na.rm = TRUE)
# MSE.test.iid <- mean((dat.without.edge.test$trees-dat.without.edge.test$prediction_iid)^2, na.rm = TRUE)
MSE.train.sp <- mean((dat.without.edge.train$trees-dat.without.edge.train$prediction_sp)^2, na.rm = TRUE)
MSE.test.sp <- mean((dat.without.edge.test$trees-dat.without.edge.test$prediction_sp)^2, na.rm = TRUE)
# MSE.train.sp1 <- mean((dat.without.edge.train$trees-dat.without.edge.train$prediction_sp1)^2, na.rm = TRUE)
# MSE.test.sp1 <- mean((dat.without.edge.test$trees-dat.without.edge.test$prediction_sp1)^2, na.rm = TRUE)
# MSE.train.st <- mean((dat.without.edge.train$trees-dat.without.edge.train$prediction_st)^2, na.rm = TRUE)
# MSE.test.st <- mean((dat.without.edge.test$trees-dat.without.edge.test$prediction_st)^2, na.rm = TRUE)
# MSE.train.zib <- mean((dat.without.edge.train$trees-dat.without.edge.train$prediction_zib)^2, na.rm = TRUE)
# MSE.test.zib <- mean((dat.without.edge.test$trees-dat.without.edge.test$prediction_zib)^2, na.rm = TRUE)

sink(file = file.path(output_dir, "MSE_INLA.txt"))
cat("training: time 0-8; test: time 9-10.\n")
# cat("MSE of training dataset (iid): ")
# MSE.train.iid
# cat("MSE of test dataset (iid): ")
# MSE.test.iid
cat("MSE of training dataset (sp): ")
MSE.train.sp
cat("MSE of test dataset (sp): ")
MSE.test.sp
# cat("MSE of training dataset (sp1): ")
# MSE.train.sp1
# cat("MSE of test dataset (sp1): ")
# MSE.test.sp1
# cat("MSE of training dataset (st): ")
# MSE.train.st
# cat("MSE of test dataset (st): ")
# MSE.test.st
# cat("MSE of training dataset with zib: ")
# MSE.train.zib
# cat("MSE of test dataset with zib: ")
# MSE.test.zib
sink(file = NULL)

# sum_withinCI_iid <- rep(0, nrow(dat))
# for (j in 1:nrow(dat)) {
#   if(dat[j,"trees"] <= dat[j,"upper_CI_3_iid"] & dat[j,"trees"] >= dat[j,"lower_CI_3_iid"])
#     sum_withinCI_iid[j] <- 1
# }

sum_withinCI_sp <- rep(0, nrow(dat))
for (j in 1:nrow(dat)) {
  if(dat[j,"trees"] <= dat[j,"upper_CI_3_sp"] & dat[j,"trees"] >= dat[j,"lower_CI_3_sp"])
    sum_withinCI_sp[j] <- 1
}

# sum_withinCI_sp1 <- rep(0, nrow(dat))
# for (j in 1:nrow(dat)) {
#   if(dat[j,"trees"] <= dat[j,"upper_CI_3_sp1"] & dat[j,"trees"] >= dat[j,"lower_CI_3_sp1"])
#     sum_withinCI_sp1[j] <- 1
# }
# 
# sum_withinCI_st <- rep(0, nrow(dat))
# for (j in 1:nrow(dat)) {
#   if(dat[j,"trees"] <= dat[j,"upper_CI_3_st"] & dat[j,"trees"] >= dat[j,"lower_CI_3_st"])
#     sum_withinCI_st[j] <- 1
# }

sink(file = file.path(uncertainty_dir, "prediction_coverage_INLA_method.txt"))
# cat("prediction coverage in training data (iid): ")
# mean(sum_withinCI_iid[1:(N*N*(TRAIN+1))])
# cat("prediction coverage in test data (iid): ")
# mean(sum_withinCI_iid[(N*N*(TRAIN+1)+1):nrow(dat)])
cat("prediction coverage in training data (sp): ")
mean(sum_withinCI_sp[1:(N*N*(TRAIN+1))])
cat("prediction coverage in test data (sp): ")
mean(sum_withinCI_sp[(N*N*(TRAIN+1)+1):nrow(dat)])
# cat("prediction coverage in training data (sp1): ")
# mean(sum_withinCI_sp1[1:(N*N*(TRAIN+1))])
# cat("prediction coverage in test data (sp1): ")
# mean(sum_withinCI_sp1[(N*N*(TRAIN+1)+1):nrow(dat)])
# cat("prediction coverage in training data (st): ")
# mean(sum_withinCI_st[1:(N*N*(TRAIN+1))])
# cat("prediction coverage in test data (st): ")
# mean(sum_withinCI_st[(N*N*(TRAIN+1)+1):nrow(dat)])
sink(file = NULL)