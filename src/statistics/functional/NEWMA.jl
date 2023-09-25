using Distributions
using LinearAlgebra
using Parameters
using StatsBase

"""
    NEWMA

Nonparametric Exponentially Weighted Moving Average object for monitoring nonparametric regression estimates.

# Fields
- `λ::Float64`: The EWMA smoothing constant. Must be between 0 and 1.
- `value::Float64`: The current value of the NEWMA statistic.
- `σ::Float64`: The standard deviation of the error term.
- `cdf_σ::E`: A function that computes the (estimated) cdf of the error term `σ` with signature `cdf_σ(x::Real)`.
- `g::F`: The estimated regression function object. Must have a method of signature `predict(g::F, x::AbstractVector)`.
- `Ej::Vector{Float64}`: The current smoothed observations.

# References
Zou, C., Tsung, F., & Wang, Z. (2008). Monitoring Profiles Based on Nonparametric Regression Methods. Technometrics, 50(4), 512-526. https://doi.org/10.1198/004017008000000433
"""
@with_kw mutable struct NEWMA{E,F,M} <: UnivariateStatistic
    λ::Float64              # EWMA smoothing constant
    value::Float64 = 0.0
    σ::Float64              # Standard deviation of the error term
    cdf_σ::E                # Cumulative distribution function of σ
    V::M                    # Covariance matrix of the nonparametric regression
    g::F                    # Estimated regression function object, must have a method
                            # of signature `predict(g::F, x::AbstractVector)`
    Ej::Vector{Float64}     # Current smoothed observations
    @assert 0 < λ <= 1
end
export NEWMA

get_design(stat::NEWMA) = [stat.λ]
set_design!(stat::NEWMA, λ::AbstractVector) = stat.λ = first(k)
set_design!(stat::NEWMA, λ::Float64) = stat.λ = λ


function NEWMA(λ::Float64, g::F, σ::Float64, cdf_σ::E, n::Int)
    return NEWMA(λ=λ, value=0.0, g = g, σ = σ, cdf_σ = cdf_σ, Ej = zeros(n))
end

#FIXME: test
"""
    NEWMA(λ::Float64, g::F, dat::FunctionalData)

Construct a NEWMA control chart by calculating the standard errors for a functional regression model.

# Arguments
- `λ::Float64`: The EWMA smoothing constant. Must be between 0 and 1.
- `g::F`: The estimated regression function object. Must have a method of signature `predict(g::F, x::AbstractVector)`.
- `dat::FunctionalData`: A `FunctionalData` object containing observations of the regression curves.
- `Ej::Vector{Float64}`: The current smoothed observations.

# Returns
The constructed NEWMA control chart.

# Examples
    using Loess, Plots
    xs = 10 .* rand(100)
    ys = sin.(xs) .+ rand(100)
    g = loess(xs, ys)
    dat = FunctionalData(xs, ys)
    NEWMA(0.2, g, dat)
"""
function NEWMA(λ::Float64, g::F, dat::FunctionalData)
    error("Must specify matrix Sigma")
    σjs = zeros(length(dat))
    res = zeros(length(get_covariates(first(dat))))
    for i in 1:length(dat)
       res[:] = dat[i].y - predict(g, get_covariates(dat[i])) 
       σjs[i] = std(res)
    end
    n = length(get_covariates(first(dat)))
    ecdf_σ = ecdf(σjs)
    error("Add regression covariance matrix")
    return NEWMA(λ, g, σ, ecdf_σ, n + 1)
end

update_statistic(STAT::NEWMA, x::FunctionalObservation) = update_statistic!(deepcopy(STAT), x)

#FIXME: test
function update_statistic!(STAT::NEWMA, x::FunctionalObservation)
    yhat = predict(STAT.g, get_covariates(x))
    zj = (x.y .- yhat) ./ STAT.σ
    new_σ = quantile(Normal(0,1), STAT.cdf_σ(std(zj)))
    STAT.Ej[1:(end-1)] = (1 - STAT.λ) * STAT.Ej[1:(end-1)] + STAT.λ * zj
    STAT.Ej[end] = (1 - STAT.λ) * STAT.Ej[end] + STAT.λ * new_σ
    error("Add regression covariance matrix")
    STAT.value = dot(STAT.Ej, STAT.Sigma, STAT.Ej)
end