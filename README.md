# StatisticalProcessMonitoring.jl: Statistical Process Monitoring in Julia

| Build | Coverage | Documentation |
|-------|----------|---------------|
| [![Build status](https://github.com/DedZago/StatisticalProcessMonitoring.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DedZago/StatisticalProcessMonitoring.jl/actions/workflows/CI.yml?query=branch%3Amain)| [![codecov](https://codecov.io/gh/DedZago/StatisticalProcessMonitoring.jl/graph/badge.svg?token=F1KFUFLD9A)](https://codecov.io/gh/DedZago/StatisticalProcessMonitoring.jl)| [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DedZago.github.io/StatisticalProcessMonitoring.jl/stable/)|



## Overview

**StatisticalProcessMonitoring.jl** is a comprehensive package for Statistical Process Monitoring, which provides users with tools for monitoring the stability of sequential processes.

## Package Features

This package implements a number of univariate and multivariate control charts, as well as control charts for monitoring partially-observed processes, profiles, and support for multi-chart designs.

Various type of control limits are implemented, with dedicated algorithms for estimating their values based on common metrics such as in-control average run length and in-control run length quantiles.

Hyperparameter tuning is supported with a native grid search implementation, as well as black-box optimization.

To ensure the package loads quickly, two extensions are [conditionally loaded](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)) only if needed:

* **PlottingSPM.jl**: functionalities for plotting the results of applying a control chart to a dataset. This extension can be loaded by running `using Plots` after loading the **StatisticalProcessMonitoring.jl** package.
* **NLoptSPM.jl**: functionalities for hyperparameter tuning using the [NLopt](https://github.com/JuliaOpt/NLopt.jl) package. This extension can be loaded by running `using NLopt` after loading the **StatisticalProcessMonitoring.jl** package.

The package is highly extensible and can incorporate custom monitoring statistics.
