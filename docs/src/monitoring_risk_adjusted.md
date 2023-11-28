# Monitoring risk-adjusted surgical outcomes

```julia
using StatisticalProcessMonitoring, Distributions, Random, Parameters, CSV, DataFrames, CategoricalArrays, MixedModels, Plots
```

```julia
dat = CSV.read("cardiacsurgery.csv", DataFrame)
dat.surgeon = categorical(dat.surgeon)
dat_ic = dat[dat.date .<= 730, :]
dat_oc = dat[730 .< dat.date .<= 1095, :]
mod = fit(MixedModel, @formula(status ~ Parsonnet + (1|surgeon)), dat_ic, Bernoulli())
```

```julia
Random.seed!(239184367)
STAT = RiskAdjustedCUSUM(Δ = 0.75, model = mod, response=:status)
LIM = OneSidedFixedLimit(1.0, true)
NOM = ARL(1000)
PH2 = Phase2(Bootstrap(), dat_ic)

CH = ControlChart(STAT, LIM, NOM, PH2)
saCL!(CH, verbose=true, gamma=0.05)
```

```julia
proc = apply_chart(CH, dat_oc)
```

```julia
plt = plot_series(proc, dpi=300, label="")
```

![example-risk-adjusted](./assets/img/example-risk-adjusted.png)
