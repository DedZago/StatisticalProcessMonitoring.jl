# Joint monitoring of mean and covariance matrix

```julia
using DrWatson
@quickactivate
using SPM, LinearAlgebra, Random, Distributions, CSV, DataFrames, Plots
```

```julia
seed = 54397858713
Random.seed!(seed)
# Normal distribution given by example 7.7 from Qiu (2013)
μ = [0, 0, 0]
Σ = [1.0 0.8 0.5; 0.8 1.0 0.8; 0.5 0.8 1.0]
Σ_sqm1 = inv(sqrt(Σ))
DIST = MultivariateNormal(μ, Σ)
NM = QRL(200, 0.5)
STATS = (MCUSUM(k=0.25, p=3), MCUSUM(k=0.5, p=3), MCUSUM(k=1.0, p=3), MEWMC(λ=0.2, p = 3))
EST_STATS = Tuple(LocationScaleStatistic(s, μ, Σ_sqm1) for s in STATS)
LIMS = Tuple(OneSidedFixedLimit(1.0, true) for _ in 1:4)
PH2 = Phase2Distribution(DIST)
CH = ControlChart(EST_STATS, LIMS, NM, PH2)
bootstrapCL!(CH)
```

```julia
data = CSV.read(datadir("example77.csv"), DataFrame)
xnew = Matrix(data)[:, 1:3]
proc = apply_chart(CH, xnew)
subtitles = ["MCUSUM k = 0.25", "MCUSUM k = 0.5", "MCUSUM k = 1", "MEWMC λ = 0.2"]
plt = plot_series(proc, dpi=300, subtitles=subtitles)
```

```julia
save(plotsdir("example-multiple-mean-covariance.png"), plt)
```

![example-multiple-mean-covariance](./assets/img/example-multiple-mean-covariance.png)
