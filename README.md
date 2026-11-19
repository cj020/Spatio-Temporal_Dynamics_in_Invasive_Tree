# A Bayesian Mechanistic Model of Spatiotemporal Dynamics in Invasive Tree Species

## 1. Overview

This repository contains the code and data used in the manuscript **"A Bayesian Mechanistic Model of Spatiotemporal Dynamics in Invasive Tree Species"** by Jiaqi Chen, Huijing Du, and Yawen Guan.

We propose a Bayesian mechanistic framework for modeling the spatiotemporal dynamics of invasive tree cover. The framework combines a diffusion–logistic growth partial differential equation (PDE) with a zero-inflated binomial observation model to account for both the ecological processes governing tree cover dynamics and the excess zeros in the observed data.

We compare the proposed PDE-based model with three benchmark approaches:

- Generalized Additive Models (GAM)
- Regression Trees (RT)
- Integrated Nested Laplace Approximation (INLA)

The repository provides the data-processing scripts, model implementation, prediction procedures, and code for fitting the benchmark models.

## 2. Example Data

For computational efficiency and ease of demonstration, this repository includes a **100 × 100 grid example dataset**. The analysis presented in the manuscript uses the full **300 × 300 grid dataset**.

Because of the difference in spatial resolution, results obtained using the example dataset may differ slightly from those reported in the manuscript.

## 3. Requirements

**R version:** 4.5

**Required R packages:**

`raster`, `fields`, `RColorBrewer`, `rasterVis`, `mvtnorm`, `gridExtra`, `fBasics`, `coda`, `tidyverse`, `egg`, `Rcpp`, `RcppArmadillo`, `Matrix`, `mgcv`, `caret`, `rpart`, `rpart.plot`, `spdep`, and `INLA`.

Note: INLA may require installation from its official repository rather than directly from CRAN.

## 4. Repository Structure

The repository is organized into the following main folders:

| Folder | Description |
|---|---|
| `raw_data` | Annual tree cover estimates from 1985 to 2020 obtained from the Rangeland Analysis Platform (RAP). The `examples` subfolder contains sample datasets for 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2019, and 2020. |
| `dataset_generating` | R scripts for generating the study-region dataset from the raw data. |
| `functions` | General-purpose R functions used in the analysis, including numerical PDE solvers and homogenization procedures. |
| `small_region` | Working directory for the analysis of the study region, including the processed datasets and scripts for model fitting, prediction, and benchmark comparisons. |

## 5. Main Scripts

The following scripts implement the main steps of the analysis.

| Script | Description |
|---|---|
| `dataset_generating/dataset_generating_small_area.R` | Generates the dataset for the study region from the raw tree cover data. |
| `small_region/data_cleaning_real_data_small_area.R` | Performs data cleaning and prepares the dataset for model fitting. |
| `small_region/PDE_based_method_small_region.R` | Implements the proposed Bayesian PDE-based mechanistic model, including model fitting and prediction. |
| `small_region/GAM_small_region.R` | Fits the Generalized Additive Model (GAM) benchmark. |
| `small_region/RT_small_region.R` | Fits the Regression Tree (RT) benchmark. |
| `small_region/INLA_small_region.R` | Fits the INLA benchmark. |

## 6. General Workflow

To reproduce the analysis, the scripts should generally be executed in the following order:

1. **Data generation:** Run `dataset_generating_small_area.R` to generate the study-region dataset from the raw data.
2. **Data preprocessing:** Run `data_cleaning_real_data_small_area.R` to clean and prepare the data.
3. **PDE-based model:** Run `PDE_based_method_small_region.R` to fit the proposed model and obtain predictions.
4. **Benchmark models:** Run the GAM, RT, and INLA scripts to fit the corresponding benchmark models and generate predictions.
5. **Model comparison:** Compare the predictive performance of the proposed model with the benchmark approaches.

The general-purpose functions required by these scripts are provided in the `functions` folder.

## 7. Citation

If you use this code or data in your research, please cite the associated manuscript:

**Chen, J., Du, H., and Guan, Y.** *A Bayesian Mechanistic Model of Spatiotemporal Dynamics in Invasive Tree Species.*

[Add the journal or preprint information and DOI when available.]

## 8. License

[Specify the repository license here.]

## 9. Contact

For questions about the code or the analysis, please contact:

**Jiaqi Chen**
