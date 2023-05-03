using FunctionWrappers
import FunctionWrappers.FunctionWrapper

abstract type OneSidedLimit <: AbstractLimit end
abstract type TwoSidedLimit <: AbstractLimit end

function is_IC(L::AbstractLimit, stat::AbstractStatistic)
    val = get_value(stat)
    lim = get_value(L)
    @assert typeof(val) == typeof(lim)
    return compare_values(lim, val, L)
end

function is_IC(L::AbstractLimit, t, stat::AbstractStatistic)
    val = get_value(stat)
    lim = get_value(L, t, stat)
    @assert typeof(val) == typeof(lim)
    return compare_values(lim, val, L)
end

#FIXME: test
function compare_values(lim_val, stat_val, L::LIM) where LIM <: TwoSidedLimit
    if (stat_val > lim_val) || (stat_val < -lim_val)
        return false
    end
    return true
end

function compare_values(lim_val, stat_val, L::LIM) where LIM <: OneSidedLimit
    for i in eachindex(lim_val)
        if L.upw[i]
            stat_val[i] < lim_val[i] || return false
        else
            stat_val[i] > lim_val[i] || return false
        end
    end
    return true
end

"""
    OneSidedFixedLimit(value::Float64, upw::Bool)

Classical fixed one-sided limit, such that the run length ``RL`` of a control chart is the first time ``t`` in which the statistic ``C_t`` crosses the limit.

* if `upw == true`, ``RL = \\inf\\{t : C_t > value\\}``
* if `upw == false`, ``RL = \\inf\\{t : C_t < -value\\}``

Note that `value > 0` by the way it is defined.
"""
@with_kw mutable struct OneSidedFixedLimit{T} <: OneSidedLimit
    value::T
    upw::Bool = true
    @assert value > 0.0
end
export OneSidedFixedLimit

"""
    TwoSidedFixedLimit(value::Float64)

Classical fixed two-sided limit, such that the run length ``RL`` of a control chart is the first time ``t`` in which the statistic ``C_t`` crosses the limit:

``RL = \\inf\\{t > 0 : |C_t| > value\\}``.

Note that `value > 0` by the way it is defined.
"""
@with_kw mutable struct TwoSidedFixedLimit{T} <: TwoSidedLimit
    value::T
    @assert value > 0.0
end
export TwoSidedFixedLimit


"""
    OneSidedCurvedLimit(value::Float64, upw::Bool)
    OneSidedCurvedLimit(value::Vector{T}, upw::Vector{Bool})

Curved one-sided limit, such that the run length ``RL`` of a control chart is the first time ``t`` in which the statistic ``C_t`` crosses the limit.

* if `upw == true`, ``RL = \\inf\\{t : C_t > value\\cdot f(t)\\}``
* if `upw == false`, ``RL = \\inf\\{t : C_t < -value\\cdot f(t)\\}``

Note that `value > 0` by the way it is defined.
"""
#FIXME: test
@with_kw mutable struct OneSidedCurvedLimit{T, S} <: OneSidedLimit
    value::T
    upw::Bool = true
    fun::FunctionWrapper{Float64, Tuple{Float64, AbstractStatistic}}

    OneSidedCurvedLimit(value, upw::Bool, f::Function, stat::AbstractStatistic) = new{typeof(first(value)), typeof(stat)}(value, upw, FunctionWrapper{typeof(first(value)), Tuple{typeof(first(value)), typeof(stat)}}(f))
end
export OneSidedCurvedLimit

get_value(::OneSidedCurvedLimit) = error("Curved limit requires specifying a time and a statistic.")
get_value(L::OneSidedCurvedLimit, t, stat) = L.value * L.fun(t, stat)


"""
    TwoSidedCurvedLimit(value::Float64)
    TwoSidedCurvedLimit(value::Vector{T})

Curved one-sided limit, such that the run length ``RL`` of a control chart is the first time ``t`` in which the statistic ``C_t`` crosses the limit.

``RL = \\inf\\{t > 0 : |C_t| > value\\cdot f(t)\\}``.

Note that `value > 0` by the way it is defined.
"""
@with_kw mutable struct TwoSidedCurvedLimit{T, S} <: TwoSidedLimit
    value::T
    fun::FunctionWrapper{Float64, Tuple{Float64, AbstractStatistic}}

    TwoSidedCurvedLimit(value, f::Function, stat::AbstractStatistic) = new{typeof(first(value)), typeof(stat)}(value, FunctionWrapper{typeof(first(value)), Tuple{typeof(first(value)), typeof(stat)}}(f))
end
export TwoSidedCurvedLimit

get_value(::TwoSidedCurvedLimit) = error("Curved limit requires specifying a time and a statistic.")
get_value(L::TwoSidedCurvedLimit, t, stat) = L.value * L.fun(t, stat)