"""
    OneSidedLimit(value::Float64, upw::Bool)
    OneSidedLimit(value::Vector{T}, upw::Vector{Bool})

Classical fixed one-sided limit, such that the run length ``RL`` of a control chart is the first time ``t`` in which the statistic ``C_t`` crosses the limit.

* if `upw == true`, RL = \\inf\\{t : C_t > value\\}``
* if `upw == false`, RL = \\inf\\{t : C_t < -value\\}``

Note that `value > 0` by the way it is defined.
"""
@with_kw mutable struct OneSidedLimit{T} <: AbstractLimit
    value::Vector{T}
    upw::Vector{Bool} = ones(length(value))
    @assert length(value) == length(upw)
end
export OneSidedLimit

OneSidedLimit(h::Float64; upw = true) = OneSidedLimit([h], [upw])

function is_IC(L::OneSidedLimit, stat::AbstractStatistic)
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
    TwoSidedLimit(value::Float64)
    TwoSidedLimit(value::Vector{T})

Classical fixed two-sided limit, such that the run length ``RL`` of a control chart is the first time ``t`` in which the statistic ``C_t`` crosses the limit:

``RL = \\inf\\{t > 0 : |C_t| > value\\}``.

Note that `value > 0` by the way it is defined.
"""
@with_kw mutable struct TwoSidedLimit{T} <: AbstractLimit
    value::Vector{T}
    @assert all(value .> 0.0)
end
export TwoSidedLimit

TwoSidedLimit(h::Float64) = TwoSidedLimit([h])

function is_IC(L::TwoSidedLimit, stat::AbstractStatistic)
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

