
# A Bayesian Mechanistic Model of Spatiotemporal Dynamics in Invasive Tree Species

## Summary
This repository contains the code and data used in the manuscript **"A Bayesian Mechanistic Model of Spatiotemporal Dynamics in Invasive Tree Species"** by Jiaqi Chen, Huijing Du and Yawen Guan. We propose a Bayesian mechanistic framework for invasive tree cover that combines an ecological diffusion–logistic growth partial differential equation (PDE) with a zero-inflated binomial observation model. Then compare our model with GAM, RT and INLA.

## Example data note. 
The repository ships a 100 × 100 example for speed, the paper uses 300 × 300, and results will differ slightly.

## Requirements
R version: 4.5.

R packages: raster, fields, RColorBrewer, rasterVis, mvtnorm, gridExtra, fBasics, coda, tidyverse, egg, Rcpp, RcppArmadillo, Matrix, mgcv, caret, rpart, rpart.plot, spdep, INLA.

## Folder guide
`raw_data`: The `examples` folder contains sample datasets for the years 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2019, and 2020.

`dataset_generating`:

`functions`: 

`small_region` (with its data, output, result).

## Script guide
One line per script in small_region saying what it does.

## Function guide.
One line per function in functions with inputs and output.

## Run order. 
Clean data, fit the PDE model, predict, fit the benchmarks, compare, with rough run times.

## Paper map. Which script produces each table and figure.

## Citation, license, contact.

## Dataset and Spatial Resolution

The data cover a 9 km × 9 km subregion near the Niobrara River in the Nebraska Sandhills.
The response is annual tree cover (integer percentage, 0 to 100) from the Rangeland
Analysis Platform for 1985–2020. The covariates are distance to the nearest stream and
elevation.

The data are provided at two spatial resolutions. The 300 × 300 grid (30 m pixels) is
used for the analysis in the manuscript. The 100 × 100 grid (90 m pixels) is obtained by
averaging and is provided for faster runs to test the code.

## Code and Analysis

The following scripts correspond to the main steps of data preparation and statistical analysis:

1. **Create the dataset for the study region from the raw CSV files:**
   `dataset_generating/dataset_generating_small_area.R`

2. **Prepare the data for model fitting:**  
   `small_region/data_cleaning_real_data_small_area.R`

3. **PDE-based method:**  
   `small_region/PDE_based_method_small_region.R`

4. **Generalized additive model (GAM):**  
   `small_region/GAM_small_region.R`

5. **Regression tree (RT):**  
   `small_region/RT_small_region.R`

6. **Integrated Nested Laplace Approximation (INLA):**  
   `small_region/INLA_small_region.R`

The scripts can be run individually to reproduce the corresponding data preparation steps and model analyses.
