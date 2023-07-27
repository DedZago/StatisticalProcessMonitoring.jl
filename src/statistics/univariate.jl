abstract type UnivariateStatistic <: AbstractStatistic end
export UnivariateStatistic


"""
    Shewhart(value)

Shewhart control chart with initial value `value`.

The update mechanism based on a new observation `x` is given by

``value = x``.

### References 
* Shewhart, W. A. (1931). Economic Control of Quality Of Manufactured Product. D. Van Nostrand Company.
"""
@with_kw mutable struct Shewhart{V} <: UnivariateStatistic 
    value::V = 0.0
    @assert !isinf(value)
end
export Shewhart

get_design(stat::Shewhart) = NamedTuple()
set_design!(stat::Shewhart, ::Float64) = error("Cannot set a design for Shewhart chart.")

update_statistic(stat::Shewhart, x::Real) = x
update_statistic!(stat::Shewhart, x::Real) = stat.value = update_statistic(stat, x)



"""
    EWMA(λ, value)

Exponentially weighted moving average with design parameter `λ` and initial value `value`.

The update mechanism based on a new observation `x` is given by

``value = (1-λ)*value + λ * x``.

### References 
* Roberts, S. W. (1959). Control Chart Tests Based on Geometric Moving Averages. Technometrics, 1(3), 239-250. https://doi.org/10.1080/00401706.1959.10489860
"""
@with_kw mutable struct EWMA{L,V} <: UnivariateStatistic 
    λ::L = 0.1
    value::V = 0.0
    @assert 0.0 < λ <= 1.0
    @assert !isinf(value)
end
export EWMA


get_design(stat::EWMA) = (λ = stat.λ,)
set_design!(stat::EWMA, λ::Float64) = stat.λ = λ
set_design!(stat::EWMA, λ) = stat.λ = λ[1]


update_statistic(stat::EWMA, x::Real) = (1.0 - stat.λ) * stat.value + stat.λ * x
update_statistic!(stat::EWMA, x::Real) = stat.value = update_statistic(stat, x)


"""
    CUSUM(k, value, upw::Bool)

CUSUM statistic with design parameter `k` and initial value `value`.

The update mechanism based on a new observation `x` is given by:
* if `upw == true`, then ``value = \\max\\{0, value + x - k\\}``;
* if `upw == false`, then ``value = \\min\\{0, value + x + k\\}``.

### References 
* Page, E. S. (1954). Continuous Inspection Schemes. Biometrika, 41(1/2), 100. https://doi.org/10.2307/2333009
"""
@with_kw mutable struct CUSUM{K,V} <: UnivariateStatistic 
    k::K = 1.0
    value::V = 0.0
    upw::Bool = true
    @assert !isinf(value)
    @assert (k > 0.0) && !isinf(k)
end
export CUSUM


get_design(stat::CUSUM) = (k = stat.k,)
set_design!(stat::CUSUM, k::Float64) = stat.k = k
set_design!(stat::CUSUM, k) = stat.k = k[1]


function update_statistic(stat::CUSUM, x::Real)
    if stat.upw
        return max(0.0, stat.value + x - stat.k)
    else
        return min(0.0, stat.value + x + stat.k)
    end
end


update_statistic!(stat::CUSUM, x::Real) = stat.value = update_statistic(stat, x)


"""
    AEWMA(λ, k, value)

Adaptive exponentially weighted moving average with design parameters `λ`, `k`, and initial value `value`.

The update mechanism based on a new observation `x` is given by

``value = (1-phi(e))*value + phi(e) * x``,

where `phi(e)` is a forecast error function based on the Huber function.

### References 
Capizzi, G. & Masarotto, G. (2003). An Adaptive Exponentially Weighted Moving Average Control Chart. Technometrics, 45(3), 199-207.
"""
@with_kw mutable struct AEWMA{L,V} <: UnivariateStatistic 
    λ::L = 0.1
    k::L = 3.0
    value::V = 0.0
    @assert 0.0 < λ <= 1.0
    @assert 0.0 < k
    @assert !isinf(value)
end
export AEWMA
#TODO: test


get_design(stat::AEWMA) = (λ = stat.λ, k = stat.k)

function set_design!(stat::AEWMA, par)
    @assert length(par) == 2
    stat.λ = par[1]
    stat.k = par[2]   
    return par
end
 
function huber(e, λ, k)
    if abs(e) <= k
        return λ*e
    elseif e > k
        return e - (1-λ)*k
    else
        return e + (1-λ)*k
    end
end

update_statistic(stat::AEWMA, x::Real) = stat.value + huber(x - stat.value, stat.λ, stat.k)
update_statistic!(stat::AEWMA, x::Real) = stat.value = update_statistic(stat, x)