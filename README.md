# SPM.jl: Statistical Process Monitoring in Julia

| Build | Coverage | Documentation |
|-------|----------|---------------|
| [![Build status](https://github.com/DedZago/SPM.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DedZago/SPM.jl/actions/workflows/CI.yml?query=branch%3Amain)| [![codecov](https://codecov.io/gh/DedZago/SPM.jl/graph/badge.svg?token=F1KFUFLD9A)](https://codecov.io/gh/DedZago/SPM.jl)| [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DedZago.github.io/SPM.jl/dev/)|

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DedZago.github.io/SPM.jl/stable/) -->



## Overview

**SPM.jl** is a comprehensive package for Statistical Process Monitoring, which provides users with tools for monitoring the stability of sequential processes.

## Package Features

This package provides a number of univariate and multivariate control charts, as well as control charts for monitoring partially-observed processes, profiles, and support for multi-chart designs.

Various type of control limits are implemented, with dedicated algorithms for estimating their values based on common metrics such as in-control average run length and in-control run length quantiles.

Hyperparameter tuning is supported with a native grid search implementation, as well as black-box optimization using the [NLopt](https://github.com/JuliaOpt/NLopt.jl) package.

The package is highly extensible and can incorporate custom monitoring statistics.

