# Getting started

## Type hierarchy

The `StatisticalProcessMonitoring.jl` package introduces a flexible definition of a control chart through the `ControlChart` type. The attributes of `ControlChart` determine its main properties.

### ControlChart type
```julia
mutable struct ControlChart{STAT, LIM, NOM, PH2} <: AbstractChart{STAT, LIM, NOM, PH2}
  stat::STAT
  limit::LIM
  nominal::NOM
  phase2::PH2
  t::Int
end
```
- **Attributes:**
  - `stat`: Monitoring statistic (see [Monitoring statistics](#Monitoring-statistics)).
  - `limit`: Control limit (see [Control limits](#Control-limits)).
  - `nominal`: Nominal property (see Section [Nominal properties](#Nominal-properties)).
  - `phase2`: Phase II data simulator (see [Simulating new observations](#Simulating-new-observations)).
  - `t`: Current time point indicator.

This type is defined as `mutable`, allowing updates via functions such as `update_chart!`.

```julia
function update_chart!(CH::AbstractChart, x)
  CH.t += 1
  update_statistic!(get_statistic(CH), x)
end
```

Then, a suitable implementation of `update_statistic!` will produce the desired behaviour of the monitoring statistic.

#### Example: EWMA and AEWMA statistics

To demonstrate the flexibility of the interface, consider the EWMA and AEWMA control charts that share control limits, nominal property, and phase II data, differing only in their statistics.

```julia
mutable struct EWMA <: AbstractStatistic 
  λ::Float64
  value::Float64
end
```
```julia
mutable struct AEWMA <: AbstractStatistic 
  λ::Float64
  k::Float64
  value::Float64
end
```
The `update_statistic!` function is specialized for each statistic:
```julia
function update_statistic!(stat::EWMA, x::Real)
  stat.value = (1.0 - stat.λ) * stat.value + stat.λ * x
end

function update_statistic!(stat::AEWMA, x::Real)
  stat.value = stat.value + huber(x - stat.value, stat.λ, stat.k)
end
```
Here, the `huber` function implements Huber's score function.
Then, as an example, we can define two control charts as:

```julia
CH_E  = ControlChart(EWMA(λ = 0.2),
                     TwoSidedFixedLimit(0.25),
                     ARL(200),
                     Phase2Distribution(Normal(0,1)
                     )

CH_AE = ControlChart(AEWMA(λ = 0.2, k=3.0),
                     TwoSidedFixedLimit(0.25),
                     ARL(200),
                     Phase2Distribution(Normal(0,1)
                     )
```

The two charts are defined with the same control limit, nominal properties, and Phase 2 object. However, their behaviour will be that of the EWMA and AEWMA control charts, respectively.

#### Implementation of multi-chart monitoring schemes

Multi-chart monitoring is supported by the `MultipleControlChart` type alias.
```julia
const MultipleControlChart{S,L,N,P} = ControlChart{S,L,N,P} where {S<:Tuple,L<:Tuple,N,P}
```

For example, a multi-chart statistic composed of two EWMA charts run simultaneously for monitoring normal data can be defined as

```julia
ControlChart(
    (EWMA(λ = 0.05), EWMA(λ=0.2)),
    (TwoSidedFixedLimit(0.25), TwoSidedFixedLimit(0.5)),
    ARL(200),
    Phase2Distribution(Normal(0,1))
    )
```

## Control limits

#### Types of control limits

Control limits $(\text{LCL}_t, \text{UCL}_t)$ can be fixed or dynamic (time-varying).

##### Fixed control limits

Defined using constant boundaries. For $h > 0$,

$\text{LCL}_t = -h, \quad \text{UCL}_t = h \quad \text{for all  } t.$


Implemented as `TwoSidedFixedLimit` and `OneSidedFixedLimit`. The `OneSidedFixedLimit` allows choosing limits of the form $(0, h)$, if its `upw` attribute is set to `true`, or $(-h, 0)$ if set to `false`.

##### Deterministic time-varying control limits

For example, time-varying limits for the EWMA chart with fast initial response:

$\text{LCL}_t = h \cdot g_l(t), \quad \text{UCL}_t = h \cdot g_u(t)$

Implemented as `TwoSidedCurvedLimit` and `OneSidedCurvedLimit`, which also require specification of a deterministic function $g(t)$ .

Alternatively, the `TwoSidedBootstrapLimit` and `OneSidedBootstrapLimit` allow defining control limits based on bootstrap resampling, which provide time-varying control limits with approximately constant false alarm rate

$\mathbb{P}(C_{t} \not\in (\text{LCL}_t, \text{UCL}_t) | C_{s} \in (\text{LCL}_s, \text{UCL}_s)  \text{ for all } s < t) = \alpha.$ 


## Nominal properties

Subtypes of `NominalProperties`, they define IC run length constraints for the control chart.
- `ARL` specifies the nominal average run length, $\mathbb{E}_0[\text{RL}]$.
- `QRL` specifies the nominal $p$-level quantile of the IC run length.

### Simulating new observations

Implemented in subtypes of `AbstractPhase2`.
- `Phase2Distribution` samples data from distributions using the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) package.
- `Phase2` resamples from an IC reference dataset (either a vector, a matrix or a data frame), supporting techniques like bootstrap.

## Monitoring statistics

Monitoring statistics, subtyped from `AbstractStatistic`, implement `update_statistic!` to define the update behaviour as new data is sequentially observed.

| Monitoring statistic             | Type name            | Hyperparameters                                                                 | Reference Paper       |
|----------------------------------|----------------------|--------------------------------------------------------------------------------|-----------------------|
| **Mean (univariate)**            |                      |                                                                                |                       |
|  Shewhart            | `Shewhart`           | --                                                                             | [Shewhart (1931)](https://archive.org/details/in.ernet.dli.2015.150272)       |
|  CUSUM               | `CUSUM`              | $ k \in \mathbb{R}^+ $                                                         | [Page (1954)](https://www.jstor.org/stable/2333009)           |
|  EWMA                | `EWMA`               | $ \lambda \in (0,1)$                                                           | [Roberts (1959)](https://www.tandfonline.com/doi/abs/10.1080/00401706.1959.10489860)        |
|  One-sided EWMA      | `OneSidedEWMA`       | $ \lambda\in (0,1)$                                                            | [Champ et Al. (1991)](https://www.sciencedirect.com/science/article/abs/pii/016771529190145H?via%3Dihub)          |
|  Adaptive EWMA       | `AEWMA`              | $ \lambda\in (0,1), k\in \mathbb{R}^+ $                                        | [Capizzi and Masarotto (2003)](https://www.tandfonline.com/doi/abs/10.1198/004017003000000023)        |
|  Weighted CUSUM      | `WCUSUM`             | $ \lambda \in (0,1)$, $ k \in \mathbb{R}^+ $                                   | [Shu et Al. (2008)](https://www.tandfonline.com/doi/abs/10.1080/00224065.2008.11917725)            |
| **Mean (multivariate)**          |                      |                                                                                |                       |
|  MShewhart           | `MShewhart`          | --                                                                             | [Shewhart (1931)](https://archive.org/details/in.ernet.dli.2015.150272)       |
|  MEWMA               | `DiagMEWMA`          | $ \bm{\lambda} \in (0,1)^{p}$                                                  | [Lowry et Al. (1992)](https://www.jstor.org/stable/1269551)          |
|  MCUSUM              | `MCUSUM`             | $ k\in \mathbb{R}^+ $                                                          | [Crosier (1988)](https://www.tandfonline.com/doi/abs/10.1080/00401706.1988.10488402)        |
|  Adaptive MEWMA      | `MAEWMA`             | $ \lambda\in (0,1)$, $ k \in \mathbb{R}^+ $                                    | [Mahmoud and Zahran (2010)](https://doi.org/10.1080/03610920902755813)        |
|  Adaptive MCUSUM     | `AMCUSUM`            | $ \lambda\in (0,1)$                                                            | [Dai et Al. (2011)](https://onlinelibrary.wiley.com/doi/abs/10.1002/qre.1177)            |
|  LLCUSUM             | `LLCUSUM`            | $ k \in \mathbb{R}^+ $                                                         | [Qiu (2008)](https://doi.org/10.1080/07408170701744843)            |
|  LLD                 | `LLD`                | $ \lambda\in (0,1)$                                                            | [Li et Al. (2012)](https://doi.org/10.1080/00224065.2012.11917889)             |
|  MOC                 | `MOC`                | $ \lambda\in (0,1)$                                                            | [Wang et Al. (2017)](https://doi.org/10.1080/00224065.2017.11917983)           |
| **Variance/covariance**          |                      |                                                                                |                       |
|  GLR-based statistic | `ALT`                | --                                                                             | [Alt (1984)](https://books.google.it/books/about/ENCYCLOPEDIA_OF_STATISTICAL_SCIENCES_VOL.html?id=dG7mzQEACAAJ&redir_esc=y)            |
|  MEWMS               | `MEWMS`              | $ \lambda\in (0,1)$                                                            | [Huwang et Al. (2007)](https://doi.org/10.1080/00224065.2007.11917692)         |
|  MEWMC               | `MEWMC`              | $ \lambda\in (0,1)$                                                            | [Hawkins and Maboudou-Tchao (2008)](https://www.jstor.org/stable/25471456)        |
| **Partially-observed data**      |                      |                                                                                |                       |
|  R-SADA              | `RSADA`              | $ k, \mu_\text{min}\in \mathbb{R}^+ $                       | [Xian et Al. (2019)](https://doi.org/10.1080/00224065.2019.1681924)           |
|  TRAS                | `TRAS`               | $ k, \mu_\text{min}, \Delta \in \mathbb{R}^+ $, $ r\in \{1, \ldots,q\} $ | [Liu et Al. (2015)](https://doi.org/10.1080/00401706.2014.947005)            |
| **Profile monitoring**           |                      |                                                                                |                       |
|  NEWMA               | `NEWMA`              | $ \lambda\in (0,1)$                                                            | [Zou et Al. (2008)](https://doi.org/10.1198/004017008000000433)            |
|  Risk-adjusted CUSUM | `RiskAdjustedCUSUM`  | $ \Delta\in \mathbb{R} $                                                       | [Steiner et Al. (2000)](https://academic.oup.com/biostatistics/article/1/4/441/238348)        |

*Table: List of available monitoring statistics in the **StatisticalProcessMonitoring.jl** package. Here, $p$ is the number of quality variables under monitoring. For partially-observed data, $q$ is the number of observed variables.*


### Statistics with estimated parameters

Separation of monitoring statistic and parameter estimation promotes code compartmentalization, facilitated by subtypes like `ResidualStatistic`.

For example, the definition of a monitoring statistic that behaves for $k = 1.0$ as

$C_{t} = \max\left\{ 0, \left( \frac{X_{t} - 0.5}{2} \right) - k \right\},$

can be defined using the `LocationScaleStatistic` subclass of `ResidualStatistic`.

```julia
STAT = LocationScaleStatistic(CUSUM(k=1.0), 0.5, 2.0)
```

## Control limit design

| Functions                              | Symbol         | Description                                                                                                                     | Single charts | Multi-charts |
|----------------------------------------|----------------|---------------------------------------------------------------------------------------------------------------------------------|---------------|--------------|
| `bisectionCL`, `bisectionCL!`          | `:Bisection`   | The standard bisection algorithm [Qiu (2013)](https://www.taylorfrancis.com/books/mono/10.1201/b15016/introduction-statistical-process-control-peihua-qiu). Requires an initial interval to search for the control limit.               | ✔             |              |
| `saCL`, `saCL!`                        | `:SA`          | Algorithm based on stochastic approximations. See [Capizzi & Masarotto (2016)](https://www.tandfonline.com/doi/full/10.1080/0740817X.2015.1055392) for a complete description of the tuning parameters.           | ✔             | ✔            |
| `combinedCL`, `combinedCL!`            | `:Combined`    | The standard bisection algorithm combined with a preliminary low-accuracy application of the SA algorithm to automatically select the initial interval. | ✔             |              |
| `bootstrapCL`, `bootstrapCL!`          | `:Bootstrap`   | Approximates the distribution of the monitoring statistic at each time $t = 1, 2, ..., T$ for $T > 0$, and then applies bisection search. Does not require an initial interval to begin the search. | ✔             | ✔            |

*Table: List of available algorithms for designing fixed and deterministic time-varying control limits in **StatisticalProcessMonitoring.jl**.*

### Bisection search

A common method for determining control limits such that the nominal run length characteristic (e.g. `ARL` or `QRL`) equals the nominal value.
Starting from an interval $[0, h_\text{max}]$, the method finds the appropriate control limit value via bisection, by approximating the run length characteristic at each iteration with a large number of simulated run lengths from the Phase 2 chart attribute.

This algorithm is implemented in the `bisectionCL!` and `bisectionCL` functions.

### Stochastic approximations

For single- and multi-chart schemes, the method simulates one run length at a time and applies a stochastic gradient descent algorithm,

$\bm{h}_{k+1} = \Psi\left( \bm{h}_{k} - \frac{1}{(k+1)^{q}} D \bm{s}(\bm{h}_k) \right),$

which converges to the required control limit value. 
This algorithm is implemented in the `saCL!` and `saCL` functions.

```julia
using Distributions, Random
Random.seed!(12345)

chart = ControlChart(
    (CUSUM(k = 0.25), CUSUM(k=0.5)),
    (OneSidedFixedLimit(5.0, true), OneSidedFixedLimit(5.0, true)),
    ARL(370),
    Phase2Distribution(Normal(0,1))
    )
  
saCL(chart)
(h = [7.290634106470891, 4.403721760449869], iter = 32847, status = "Convergence")
```




## Hyperparameter tuning

Hyperparameter tuning for optimal performance with respect to a specific out-of-control scenario is implemented via:

- Grid search (slow for multidimensional parameters)
- Nonlinear numerical solvers available from the [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl) package (e.g., BOBYQA)

The function `optimize_design!` provides a high-level interface for optimization and requires being able to simulate from the out-of-control scenario of interest.
An example of hyperparameter optimization can be found in the [Residual-Based Monitoring of Autocorrelated Data](monitoring_autoregressive.md) tutorial.
