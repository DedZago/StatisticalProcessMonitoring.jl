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
    z::Vector{L} = zeros(length(Λ))
    inv_Σz::Matrix{L} = inv(diagm(Λ)*diagm(Λ))
    @assert !isinf(value)
end
export MEWMA


get_design(stat::MEWMA) = (Λ = stat.Λ,)

function set_design!(stat::MEWMA, Λ::AbstractVector)
    stat.inv_Σz = inv(diagm(Λ)*diagm(Λ))
    stat.Λ = Λ
end

function update_statistic!(stat::MEWMA, x::AbstractVector)
    for j in eachindex(x)
        stat.z[j] = (1.0 - stat.Λ[j]) * stat.z[j] + stat.Λ[j] * x[j]
    end  
    stat.value = stat.z' * stat.inv_Σz * stat.z
end

update_statistic(stat::MEWMA, x::AbstractVector) = update_statistic!(deepcopy(stat), x)

