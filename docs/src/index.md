```@meta
CurrentModule = StatisticalProcessMonitoring
```

# StatisticalProcessMonitoring.jl

**StatisticalProcessMonitoring.jl**  is a Julia package for testing the stability of sequential data streams using statistical process monitoring (SPM) techniques. The package provides a flexible framework for implementing control charts, monitoring statistics, and control limit calibration methods. The package is designed to be extensible, allowing users to define custom control charts and monitoring statistics.

For an extensive review of SPM and control charts, a recommended resource is the book [Introduction to Statistical Process Control](https://www.taylorfrancis.com/books/mono/10.1201/b15016/introduction-statistical-process-control-peihua-qiu).

## Installation
StatisticalProcessMonitoring.jl requires Julia version 1.8 or above. To install StatisticalProcessMonitoring.jl, press the `]` key inside the Julia REPL to access the interactive package manager model and run the following command
```julia
pkg> add StatisticalProcessMonitoring
```

Otherwise, in the standard REPL run the following command

```julia
julia> using Pkg
julia> Pkg.add("StatisticalProcessMonitoring")
```

Once installed, the package can be loaded by running the following command

```julia
julia> using StatisticalProcessMonitoring
```

## Features overview

- The package focuses on Monte-Carlo-based calibration of control limits and hyperparameter optimization for control charts.
- The package implements a general control chart interface which can be easily extended to accommodate user-defined monitoring statistics and customized control limit behaviour.
    
1. Classical control charts
    * Univariate and multivariate Shewhart, EWMA, CUSUM charts, as well as their adaptive generalization
    * Control charts for monitoring the variance-covariance matrix.
    * Control charts based on data categorization for monitoring the process mean.
    * Control charts for partially-observed data streams.
    * Control charts for profile monitoring.
2. Multi-chart monitoring schemes
    * Support for arbitrary combination of control charts
    * Joint control limit calibration for control charts
3. Metrics
    * Metrics based on the Average Run Length
    * Metrics based on Run Length quantiles
4. Control limit calibration
    * Control limit estimation based on data distributions using [Distributions.jl](https://juliastats.org/Distributions.jl/stable/starting/)
    * Control limit estimation based on bootstrap and block bootstrap of initial data
    * Dynamic control limits with constant false-alarm rate via bootstrap resampling
5. Hyperparameter tuning
    * Optimization of control chart parameters for user-defined out-of-control scenarios, using grid search and nonlinear optimizers ([NLopt.jl](https://github.com/JuliaOpt/NLopt.jl))


# Overview of other software packages 

## Traditional control charts in SPM software

Most SPM software packages primarily focus on implementing traditional control charts for analyzing univariate or multivariate data. However, these packages may lack support for advanced methodologies like multi-chart schemes and user-defined control charts, requiring practitioners to develop their routines to address these limitations.

### Julia packages
Other packages implementing SPM tools in Julia are:

- [MultivariateAnomalies.jl](https://github.com/milanflach/MultivariateAnomalies.jl) for multivariate anomaly detection

### R packages

Some notable SPM software packages in R include:
- [qcc](https://luca-scr.github.io/qcc/) for Shewhart, CUSUM, and EWMA charts
- [MSQC](https://rdrr.io/cran/MSQC/) for multivariate control charts
- [spc](https://cran.r-project.org/web/packages/spc/index.html) for in-control average run length and run length quantile calculations
- [qicharts](https://cran.r-project.org/web/packages/qicharts/index.html) and [qicharts2](https://cran.r-project.org/web/packages/qicharts2/index.html) for run charts, Shewhart and Pareto control charts
- [edcc](https://rdrr.io/cran/edcc/) for control charts designed using the optimal economic design framework

Other specialized packages offer functionalities for specific applications, such as:
- [spcadjust](https://cran.r-project.org/web/packages/spcadjust/index.html) for adjusting control limits based on estimated parameters
- [funcharts](https://cran.r-project.org/web/packages/funcharts/index.html) for multivariate functional data, function-on-scalar regression and function-on-function regression
- [surveillance](https://cran.r-project.org/web/packages/surveillance/index.html) for change detection with a focus on public health surveillance
- [cpm](https://cran.r-project.org/web/packages/cpm/index.html), [strucchange](https://cran.r-project.org/web/packages/strucchange/index.html), [bcp](https://www.jstatsoft.org/article/view/v023i03), and [changepoint](https://cran.r-project.org/web/packages/changepoint/index.html) for process monitoring based on change-point models
- [dfphase1](https://cran.r-project.org/web/packages/dfphase1/index.html) for change detection in retrospective samples using distribution-free control charts

### Packages in other programming languages

Limited support for SPM methodologies exists in other programming languages:
- [PySpc](https://pypi.org/project/pyspc/) in Python for classical control charts
- [Pre-Screen](https://www.cpact.com/About/Software/PreScreen) in MATLAB for various control charts

For a detailed review of existing SPM packages across different programming languages, including their features and functionalities, refer to the complete article.


