```@meta
CurrentModule = SPM
```

# SPM

Documentation for [SPM](https://github.com/DedZago/SPM.jl), a package for Statistical Process Monitoring.

## Package features
1. Standard control charts
    * Univariate Shewhart, EWMA, AEWMA, and CUSUM control charts;
    * MEWMA and MAEWMA control charts;
2. Control limits
    * Classical one-sided and two-sided control limits;
    * Support for multi-chart combinations;
    * Dynamic control limits based on bootstrap and permutation methods;
3. Optimization methods
    * State-of-the-art methods for estimating control limits;
    * Optimization of control chart parameters against user-defined out-of-control scenarios;
4. Extensibility to user-made control statistics
    * Users only need to define the behaviour of the control statistic (`struct` and sequential update function), everything else is taken care of by the package.

```@index
```

```@autodocs
Modules = [SPM]
Order   = [:function, :type]
```
