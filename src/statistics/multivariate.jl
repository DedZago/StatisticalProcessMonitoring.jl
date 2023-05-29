using LinearAlgebra

"""
    MEWMA(Λ, value)

Exponentially weighted moving average with smoothing matrix `Λ` and initial value `value`.

The update mechanism based on a new observation `x` is given by

``Z = (I-Λ)*Z + Λ * x``,

and the chart value is defined as

``value = Z' Λ^(-1) Z``.

### References 
Lowry, C. A., Woodall, W. H., Champ, C. W., & Rigdon, S. E. (1992). A Multivariate Exponentially Weighted Moving Average Control Chart. Technometrics, 34(1), 46-53. https://doi.org/10.1080/00401706.1992.10485232
"""
@with_kw mutable struct MEWMA{L,V} <: UnivariateStatistic 
    Λ::Vector{L}
    value::V = 0.0
    μ::Vector{L} = zeros(L, size(L)[1])
    z::Vector{L} = deepcopy(μ)
    Σ::Matrix{L} = diagm(ones(L, size(L)[1]))
    inv_Σz::Matrix{L} = inv(diagm(Λ)*Σ*diagm(Λ))
    @assert !isinf(value)
    @assert dim(Λ[1]) == dim(Λ[2])
end
export MEWMA


get_design(stat::MEWMA) = (Λ = stat.Λ,)

function set_design!(stat::MEWMA, Λ::AbstractVector)
    stat.inv_Σz = inv(diagm(Λ)*stat.Σ*diagm(Λ))
    stat.Λ = Λ
end

function update_statistic!(stat::MEWMA, x::AbstractVector)
    for j in eachindex(x)
        stat.z[j] = (1.0 - stat.Λ[j]) * stat.z[j] + stat.Λ[j] * (x[j] - stat.μ[j])
    end  
    error("TODO")
end

function update_statistic(stat::MEWMA, x::AbstractVector)
    error("Not implemented for the MEWMA control chart")
end

