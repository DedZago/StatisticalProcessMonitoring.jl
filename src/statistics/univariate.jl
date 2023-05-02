abstract type UnivariateStatistic <: AbstractStatistic end

@with_kw mutable struct EWMA{L,V} <: UnivariateStatistic 
    λ::L = 0.1
    value::V = 0.0
    @assert 0.0 < λ <= 1.0
    @assert !isinf(value)
end
export EWMA


get_param(stat::EWMA) = (λ = stat.λ,)
set_param!(stat::EWMA, λ) = stat.λ = λ


function update_statistic!(stat::EWMA, x::Real)
    stat.value = (1.0 - stat.λ) * stat.value + stat.λ * x   
end


@with_kw mutable struct CUSUM{K,V} <: UnivariateStatistic 
    k::K = 1.0
    value::V = 0.0
    upw::Bool = true
    @assert !isinf(value)
    @assert (k > 0.0) && !isinf(k)
end
export CUSUM


get_param(stat::CUSUM) = (k = stat.k,)
set_param!(stat::CUSUM, k) = stat.k = k


function update_statistic!(stat::CUSUM, x::Real)
    if stat.upw
        stat.value = max(0.0, stat.value + x - stat.k)
    else
        stat.value = min(0.0, stat.value + x + stat.k)
    end
end