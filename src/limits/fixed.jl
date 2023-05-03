abstract type FixedLimit <: AbstractLimit


"""
    OneSidedFixedLimit(value::Float64, upw::Bool)
    OneSidedFixedLimit(value::Vector{T}, upw::Vector{Bool})

Classical fixed one-sided limit, such that the run length ``RL`` of a control chart is the first time ``t`` in which the statistic ``C_t`` crosses the limit.

* if `upw == true`, RL = \\inf\\{t : C_t > value\\}``
* if `upw == false`, RL = \\inf\\{t : C_t < -value\\}``

Note that `value > 0` by the way it is defined.
"""
@with_kw mutable struct OneSidedFixedLimit{T} <: AbstractLimit
    value::Vector{T}
    upw::Vector{Bool} = ones(length(value))
    @assert length(value) == length(upw)
end
export OneSidedFixedLimit

OneSidedFixedLimit(h::Float64; upw = true) = OneSidedFixedLimit([h], [upw])

function is_IC(L::OneSidedFixedLimit, stat::AbstractStatistic)
    val = get_value(stat)
    lim = get_value(L)
    @assert length(val) == length(lim)
    for i in eachindex(val)
        if L.upw[i]
            val[i] < lim[i] || return false
        else
            val[i] > lim[i] || return false
        end
    end
    return true
end


"""
    TwoSidedFixedLimit(value::Float64)
    TwoSidedFixedLimit(value::Vector{T})

Classical fixed two-sided limit, such that the run length ``RL`` of a control chart is the first time ``t`` in which the statistic ``C_t`` crosses the limit:

``RL = \\inf\\{t > 0 : |C_t| > value\\}``.

Note that `value > 0` by the way it is defined.
"""
@with_kw mutable struct TwoSidedFixedLimit{T} <: AbstractLimit
    value::Vector{T}
    @assert all(value .> 0.0)
end
export TwoSidedFixedLimit

TwoSidedFixedLimit(h::Float64) = TwoSidedFixedLimit([h])

function is_IC(L::TwoSidedFixedLimit, stat::AbstractStatistic)
    val = get_value(stat)
    lim = get_value(L)
    @assert length(val) == length(lim)
    for i in eachindex(val)
        if (val[i] > lim[i]) || (val[i] < -lim[i])
            return false
        end
    end
    return true
end

