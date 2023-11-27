# SPM: Statistical Process Monitoring Package

[![Build Status](https://github.com/DedZago/SPM.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DedZago/SPM.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/DedZago/SPM.jl/graph/badge.svg?token=F1KFUFLD9A)](https://codecov.io/gh/DedZago/SPM.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DedZago.github.io/SPM.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DedZago.github.io/SPM.jl/dev/)


## Overview

**SPM.jl** is a comprehensive package for Statistical Process Monitoring, which provides users with tools for monitoring the stability of sequential processes.

## Package Features

### Control Charts

- **Univariate Control Charts:**
  - Shewhart
  - Exponentially Weighted Moving Average (EWMA)
  - Adaptive EWMA (AEWMA)
  - Cumulative Sum (CUSUM)
  - Adaptive CUSUM

- **Multivariate Control Charts:**
  - Shewhart
  - Multivariate EWMA (MEWMA)
  - Multivariate Adaptive EWMA (MAEWMA)
  - Multivariate CUSUM (MCUSUM)
  - Adaptive Multivariate CUSUM (AMCUSUM)
  
- **Other Control Charts:**
  - Variance-Covariance Matrix Monitoring
  - Data Categorization-based Mean Monitoring
  - Partially-Observed Data Streams Monitoring
  - Profile Monitoring
  - Support for Multi-Chart Monitoring Schemes

### Control Limits

- **Fixed Control Limits:**
  - Classical one-sided and two-sided limits
  
- **Time-Varying Control Limits:**
  - Deterministic time-varying limits
  
- **Dynamic Control Limits:**
  - Bootstrap-based limits with constant false-alarm rate

### Control Limit Design

- **Design Methods:**
  - Classical and Improved Bisection Methods
  - Stochastic Approximation Algorithms
  - Support for In-Control Average Run Length (ARL) and Run Length Quantiles
  - Support for Multi-Chart Control Limit Design

### Hyperparameter Tuning

- **Optimization:**
  - Control Chart Parameter Optimization for User-Defined Out-of-Control Scenarios
  - Algorithms based on Grid Search and Nonlinear Optimizers

### Extensibility

- **User-Defined Monitoring Statistics:**
  - Easily extend SPM to incorporate custom monitoring statistics

