using Parameters
using StatsBase

"""
    NEWMA{F}(λ::Float64, σ::Float64, g::F) where {F <: AbstractRegressionFunction}

Construct a NEWMA (Nonparametric Exponentially Weighted Moving Average) object.

# Fields
- `λ::Float64`: The EWMA smoothing constant. Must be between 0 and 1.
- `value::Float64`: The current value of the NEWMA statistic.
- `σ::Float64`: The standard deviation of the error term.
- `g::F`: The estimated regression function object. Must have a method of signature `predict(g::F, x::AbstractVector)`.
- `Uj::Vector{Float64}`: The current smoothed observations.
"""
@with_kw mutable struct NEWMA{F} <: UnivariateStatistic
    λ::Float64              # EWMA smoothing constant
    value::Float64 = 0.0
    σ::Float64              # Standard deviation of the error term
    g::F                    # Estimated regression function object, must have a method
                            # of signature `predict(g::F, x::AbstractVector)`
    Uj::Vector{Float64}     # Current smoothed observations
    @assert 0 < λ <= 1
end
export NEWMA

get_design(stat::NEWMA) = [stat.λ]
set_design!(stat::NEWMA, λ::AbstractVector) = stat.λ = first(k)
set_design!(stat::NEWMA, λ::Float64) = stat.λ = λ


function NEWMA(λ::Float64, g::F, σ::Float64, n::Int)
    return NEWMA(λ=λ, value=0.0, g = g, σ = σ, Uj = zeros(n))
end

#FIXME: test
"""
    NEWMA(λ::Float64, g::F, dat::FunctionalData)

Construct a NEWMA control chart by calculating the standard errors for a functional regression model.

# Arguments
- `λ::Float64`: The EWMA smoothing constant. Must be between 0 and 1.
- `g::F`: The estimated regression function object. Must have a method of signature `predict(g::F, x::AbstractVector)`.
- `dat::FunctionalData`: A `FunctionalData` object containing observations of the regression curves.
- `Uj::Vector{Float64}`: The current smoothed observations.

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
    res = [d.y .- predict(g, d.x) for d in dat]     # Calculate residuals from estimated regression function
    n = first(dat).x
    σ = std(hcat(res...))                           # Assume random errors are the same across functions
    return NEWMA(λ, g, σ, n)
end

update_statistic(STAT::NEWMA, x::FunctionalObservation) = update_statistic!(deepcopy(STAT), x)

#FIXME: test
function update_statistic!(STAT::NEWMA, x::FunctionalObservation)
    yhat = predict(STAT.g, x.x)
    zj = (x.y .- yhat) ./ STAT.σ
    error("Not implemented yet.")
end