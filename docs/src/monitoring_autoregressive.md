# Monitoring residuals of an autoregressive process

```julia
using DrWatson
@quickactivate
using SPM, Distributions, Random, NLopt, Plots, Parameters, LaTeXStrings
```

```julia
import SPM.residual!, SPM.new_data!
mutable struct AR1Statistic{S} <: ResidualStatistic
    stat::S
    phi::Float64
    ym1::Float64
end

function residual!(x, S::AR1Statistic)
    yhat = x - S.phi * S.ym1   
    S.ym1 = x
    return yhat
end

@with_kw mutable struct Phase2AR1<:SPM.AbstractPhase2
    phi::Float64
    y::Float64 = 0.0
    init::Bool = false
end

function new_data!(PH2::Phase2AR1)
    if !PH2.init
        PH2.y = randn()       
        PH2.init = true
    end
    yhat = PH2.phi * PH2.y + randn()
    PH2.y = yhat
    return yhat
end
```

```julia
phi = 0.5
PH2 = Phase2AR1(phi = phi)
STAT = AR1Statistic(EWMA(λ=0.1), phi, 0.0)
NOM = ARL(500)
LIM = TwoSidedFixedLimit(1.0)
CH = ControlChart(STAT, LIM, NOM, PH2)

seed = 4398354798
Random.seed!(seed)
delta = 2.0
rlsim_oc = x -> run_sim_oc(x, shift = delta)
settings = OptSettings(verbose = false, minpar = [0.001], maxpar = [0.99])
optimize_design!(CH, rlsim_oc, settings, optimizer=:LN_BOBYQA)
```

```julia
n = 100
tau = 50
DGP = Phase2AR1(phi=phi)
y = [new_data!(DGP) for _ in 1:n]
y[(tau+1):n] .+= delta
proc = apply_chart(CH, y)
plt = plot_series(proc, dpi=300, label="")
vline!(plt, [tau], label=L"\tau", linestyle=:dot, colour="black")
```

```julia
save(plotsdir("example-autocorrelated.png"), plt)
```

![example-autocorrelated](./assets/img/example-autocorrelated.png)
