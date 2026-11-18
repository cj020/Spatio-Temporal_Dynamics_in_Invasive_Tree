
# A Bayesian Mechanistic Model of Spatiotemporal Dynamics in Invasive Tree Species

This repository contains the code and data used in the manuscript **"A Bayesian Mechanistic Model of Spatiotemporal Dynamics in Invasive Tree Species"** by Jiaqi Chen, Huijing Du and Yawen Guan.

## Dataset and Spatial Resolution

The datasets at resolutions of 300 × 300 and 100 × 100 cover the same geographic region but differ in spatial resolution. The primary analysis in the manuscript uses the 300 × 300 resolution. However, to reduce computational time, the code can also be run using the 100 × 100 resolution. The results obtained at the two resolutions are expected to be similar.

## Code and Analysis

The following scripts correspond to the main steps of data preparation and statistical analysis:

1. **Generate the study region:**  
   `dataset_generating/dataset_generating_small_area.R`

2. **Clean the data:**  
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
