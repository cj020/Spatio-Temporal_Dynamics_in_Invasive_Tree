# combine data set to generate small region data.

# Set the path to the local project directory
working_dir <- "C:/Users/YourName/Documents/Spatio-Temporal_Dynamics_in_Invasive_Tree"

# Set working directory
setwd(working_dir)

# Get the working directory
wd <- getwd()

data_dir <- file.path(wd, "real_data", "data", "N300")

dat_part1 <- read.csv(file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized_part1.csv"))
dat_part2 <- read.csv(file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized_part2.csv"))
dat_part3 <- read.csv(file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized_part3.csv"))
dat_part4 <- read.csv(file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized_part4.csv"))
dat_part5 <- read.csv(file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized_part5.csv"))
dat_part6 <- read.csv(file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized_part6.csv"))

dat <- rbind(dat_part1, dat_part2, dat_part3, dat_part4, dat_part5, dat_part6)

write.csv(dat, file.path(data_dir, "dat_8520_small_area_300_by_300_smoothed_integered_standardized.csv"))
