<br />
<div align="center">
  <a href="https://github.com/cfjimenezv07/FANOVA_FFM">
    <img src="MQ.png" alt="MQ Logo" height="150">
    <img src="UoY.png" alt="York Logo" height="150">
  </a>

<h3 align="center">Interpretable models for forecasting high-dimensional functional time series</h3>
</div>

## Abstract
<p align="justify">
We study the modeling and forecasting of high-dimensional functional time series, which can be temporally dependent and cross-sectionally correlated. We implement a functional analysis of variance (FANOVA) to decompose high-dimensional functional time series, such as subnational age- and sex-specific mortality observed over years, into two distinct components: a deterministic mean structure and a residual process varying over time. From the residual process, we implement a functional factor model to capture the remaining stochastic trends. Illustrated by the age-specific Japanese subnational mortality rates from 1975 to 2023, we evaluate and compare the accuracy of the point and interval forecasts across various forecast horizons. The results demonstrate that leveraging these interpretable components not only clarifies the underlying drivers of the data but also improves forecast accuracy, providing more transparent insights for evidence-based policy decisions.
</p>

### Main Results
The R script files in the `R Code` folder should be used in the following order:

#### 1. Setup and auxiliary functions
* **1_load_packages.R**: Installs and loads the necessary R libraries.
* **2_read_data.R**: Loads the raw datasets and formats them for functional data analysis.
* **choice_K.R**: Selection of the optimal number of components ($K$) for the functional models.
* **interval_score.R**: Computes prediction interval error metrics.

#### 2. Functional ANOVA (FANOVA) with MFTS
* **One_way_FANOVA.R**: Performs the baseline One-Way Functional ANOVA (OWA).
* **FANOVA_MFTS_Arima.R**: FANOVA with MFTS using ARIMA.
* **FANOVA_MFTS_ets.R**: FANOVA with MFTS using ETS.
* **FANOVA_MFTS_sequential.R**: Sequential estimation approach for FANOVA-MFTS.
* **FANOVA_MFTS_Split.R**: Split-sample estimation approach for FANOVA-MFTS.

#### 3. Functional Factor Model (FFM)
* **FFM.R**: Functional Factor Model implementation.

#### 4. Estimation Procedures for point and interval forecast
* **TWA-FFM.R**: Two-Way FANOVA (TWA) integrated with FFM.
* **TWA-OWA-FFM.R**: Hybrid Two-Way and One-Way FANOVA plus FFM framework.
* **TWA_FFM_sequential.R** & **TWA_FFM_split.R**: Sequential and split estimation based on TWA-FFM.
* **TWA_OWA_FFM_sequential.R** & **TWA_OWA_FFM_split.R**: Sequential and split estimation for the combined TWA-OWA-FFM.

## How to Cite
If you use this code or methodology in your research, please cite our paper:

> Shang, H. L., & Jiménez-Varón, C. F. (2026). Interpretable models for forecasting high-dimensional functional time series. *Statistics and Computing*, 36(6), Article 236. https://doi.org/10.1007/s11222-026-10993-8

**BibTeX:**
```bibtex
@article{shang2026interpretable,
  author   = {Shang, Han Lin and Jim{\'e}nez-Var{\'o}n, Cristian F.},
  title    = {Interpretable models for forecasting high-dimensional functional time series},
  journal  = {Statistics and Computing},
  volume   = {36},
  number   = {6},
  pages    = {236},
  year     = {2026},
  doi      = {10.1007/s11222-026-10993-8},
  url      = {[https://doi.org/10.1007/s11222-026-10993-8](https://doi.org/10.1007/s11222-026-10993-8)}
}