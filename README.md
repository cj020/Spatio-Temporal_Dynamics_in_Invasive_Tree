# A Bayesian Mechanistic Model of Spatiotemporal Dynamics in Invasive Tree Species

**Jiaqi Chen, Huijing Du, and Yawen Guan**

This repository fits a Bayesian ecological-diffusion and logistic-growth PDE, with a zero-inflated binomial observation model, to satellite-derived invasive tree cover. It then predicts held-out cover and compares that mechanistic model with a generalized additive model, a regression tree, and a spatial model fitted by INLA.

## Example data

The analysis uses a 300 × 300 grid. GitHub’s file-size limit required splitting that table into six CSV files under `real_data/data/N300/`. Merge those six parts into one CSV before any model is fit.

## Requirements

- R 4.5
- On Windows, [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching that R version, so `Rcpp` can compile `PDE_fun.cpp` and `H_matrix.cpp`

R packages used by the scripts:

`raster`, `fields`, `RColorBrewer`, `rasterVis`, `mvtnorm`, `gridExtra`, `fBasics`, `coda`, `tidyverse`, `egg`, `Rcpp`, `RcppArmadillo`, `Matrix`, `mgcv`, `ggplot2`, `caret`, `rpart`, `rpart.plot`, `spdep`, `patchwork`, `cowplot`, `paletteer`, `INLA`

`INLA` is installed from its own repository, not CRAN:

```r
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
```

Each script checks for missing packages and installs them when it starts. `parallel` and `grid` come with R.

## Folder guide

- `code` — R scripts and C++ sources that merge the data, fit the models, predict, and compare results.
- `real_data` — the split 300 × 300 example data and the folders where MCMC output, predictions, and comparison plots are written.
- `manuscript` — the paper PDF and the supplementary PDF.

## Script guide

- `Dataset_combining_small_region.R` — stacks the six part files into `dat_8520_small_area_300_by_300_smoothed_integered_standardized.csv`.
- `Functions.R` — shared likelihood, neighborhood, and PDE solution functions, and compiles the C++ sources.
- `PDE_fun.cpp` — numerical solution of the ecological-diffusion PDE.
- `H_matrix.cpp` — builds the homogenization matrix used by that solver.
- `MCMC_real_data.R` — fits the PDE model by MCMC on the merged 300 × 300 data.
- `Comparison_different_MCMC_initial_values.R` — reruns the MCMC from different starting values and compares the chains.
- `Prediction_PDE_real_data.R` — draws posterior predictions, MSE, and predictive coverage from the PDE fit.
- `GAM_real_data.R` — fits the generalized additive model benchmark and writes its predictions.
- `Prediction_uncertainty_GAM_real_data.R` — computes predictive uncertainty for the GAM.
- `RT_real_data.R` — fits the regression-tree benchmark and writes its predictions.
- `Prediction_uncertainty_RT_real_data.R` — computes predictive uncertainty for the regression tree.
- `INLA_real_data.R` — fits the spatial hierarchical benchmark with INLA and writes its predictions.
- `Plot_prediction_real_data.R` — maps and compares PDE, GAM, regression-tree, and INLA predictions.

## Run order

1. In every script, set `working_dir` to this repository on your machine, then run the script. The placeholder is `C:/Users/YourName/Documents/Spatio-Temporal_Dynamics_in_Invasive_Tree`.
2. Merge the data with `code/Dataset_combining_small_region.R`.
3. Fit the PDE model with `code/MCMC_real_data.R`.
4. Predict from the PDE fit with `code/Prediction_PDE_real_data.R`.
5. Fit the benchmarks with `code/GAM_real_data.R`, `code/RT_real_data.R`, and `code/INLA_real_data.R`.
6. Compare the models with `code/Plot_prediction_real_data.R`.

`Comparison_different_MCMC_initial_values.R` is optional and can be run after the main MCMC fit. The two prediction-uncertainty scripts can be run after the matching benchmark fit and before the comparison plots.
