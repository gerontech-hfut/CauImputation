# CauImputation: Causal Imputation Strategy for Missing Data

This repository contains the implementation of the Causal Imputation (CauImputation) strategy, designed for the imputation of missing data in biomedical datasets. This MATLAB codebase applies various kernel methods to improve imputation accuracy by leveraging both instance and feature information.

## Overview

Missing data is a significant challenge in data analysis, particularly in biomedical research. Traditional imputation methods often apply a uniform strategy across all instances, neglecting the essential role of the feature with the missing value. CauImputation addresses this by selecting neighbor instances for imputation based on both instance and missing feature information using a Structural Causal Model (SCM).

## Key Features

- Implements Polynomial, Linear, and RBF kernel methods for neighbor selection.
- Integrates Local Least Squares (LLS) imputation methods.
- Evaluates performance on real-world biomedical datasets.
- Provides both RMSE and MAE metrics for performance evaluation.

## Dataset Information

The repository includes the following datasets under the `data` directory:

- `completegds1761.mat`: Complete data for GDS1761.
- `completegds38.mat`: Complete data for GDS38.
- `completegds3835.mat`: Complete data for GDS3835.
- `missingdata_gds1761.mat`: Missing data for GDS1761.
- `missingdata_gds38.mat`: Missing data for GDS38.
- `missingdata_gds3835.mat`: Missing data for GDS3835.

These datasets correspond to the gene microarray datasets mentioned in the paper.

## Installation

This codebase is compatible with MATLAB 2021b. Ensure you have this version or later installed on your system.


### Running the Imputation

1. **Load the dataset**:
    ```matlab
    load('data/completegds1761.mat');
    load('data/missingdata_gds1761.mat');
    ```
2. **Choose an imputation method**:
    ```matlab
    imputed_data = CauPLLS_continuity(missing_data);
    ```
3. **Evaluate the imputed data** using RMSE or MAE:
    ```matlab
    rmse = sqrt(mean((complete_data(:) - imputed_data(:)).^2));
    mae = mean(abs(complete_data(:) - imputed_data(:)));
    ```
