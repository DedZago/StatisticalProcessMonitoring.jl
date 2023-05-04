abstract type UnivariateStatistic <: AbstractStatistic end

"""
    EWMA(λ, value)

Exponentially weighted moving average with parameter `λ` and initial value `value`.

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


get_parameter(stat::EWMA) = (λ = stat.λ,)
set_parameter!(stat::EWMA, λ) = stat.λ = λ


function update_statistic!(stat::EWMA, x::Real)
    stat.value = (1.0 - stat.λ) * stat.value + stat.λ * x   
end


"""
    CUSUM(k, value, upw::Bool)

CUSUM statistic with parameter `k` and initial value `value`.

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


get_parameter(stat::CUSUM) = (k = stat.k,)
set_parameter!(stat::CUSUM, k) = stat.k = k


function update_statistic!(stat::CUSUM, x::Real)
    if stat.upw
        stat.value = max(0.0, stat.value + x - stat.k)
    else
        stat.value = min(0.0, stat.value + x + stat.k)
    end
end