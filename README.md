
# A Bayesian Mechanistic Model of Spatiotemporal Dynamics in Invasive Tree Species

This repository contains the code and data used in the manuscript **"A Bayesian Mechanistic Model of Spatiotemporal Dynamics in Invasive Tree Species"** by Jiaqi Chen, Huijing Du and Yawen Guan.

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
