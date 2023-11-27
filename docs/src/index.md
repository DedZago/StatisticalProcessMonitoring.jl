```@meta
CurrentModule = SPM
```

# SPM

Documentation for [SPM](https://github.com/DedZago/SPM.jl), a package for Statistical Process Monitoring.

## Package features
#### Control charts
* Univariate Shewhart, EWMA, AEWMA, CUSUM, adaptive CUSUM control charts.
* Multivariate Shewhart, MEWMA, MAEWMA, MCUSUM, AMCUSUM control charts.
* Control charts for monitoring the variance-covariance matrix.
* Control charts based on data categorization for monitoring the process mean.
* Control charts for partially-observed data streams.
* Control charts for profile monitoring.
* Support for multi-chart monitoring schemes.
#### Control limits
* Classical one-sided and two-sided fixed control limits.
* Deterministic time-varying control limits.
* Dynamic control limits based on bootstrap with constat false-alarm rate.
#### Control limit design
* Classical and improved bisection methods.
* Stochastic approximation algorithms.
* Support for in-control ARL and in-control run length quantiles.
* Support for multi-chart control limit design.
#### Hyperparameter tuning
* Optimization of control chart parameters for user-defined out-of-control scenarios.
* Algorithms based on grid search and nonlinear optimizers.
#### Extensibility to user-defined monitoring statistics

