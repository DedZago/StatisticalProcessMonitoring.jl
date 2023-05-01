@with_kw mutable struct OneSidedLimit{T} <: AbstractLimit
    value::Vector{T}
    upw::Vector{Bool} = ones(length(value))
    @assert length(value) == length(upw)
end
export OneSidedLimit

OneSidedLimit(h::Float64; upw = true) = OneSidedLimit([h], [upw])

function is_IC(L::OneSidedLimit, stat::AbstractStatistic)
    #TODO: run tests to check if it is correct
    val = get_value(stat)
    lim = get_value(L)
    @assert length(val) == length(lim)
    for i in eachindex(val)
        if upw[i]
            val[i] < lim[i] || return false
        else
            val[i] > lim[i] || return false
        end
    end
    return true
end


@with_kw mutable struct TwoSidedLimit{T} <: AbstractLimit
    value::Vector{T}
    @assert all(value .> 0.0)
end
export TwoSidedLimit

TwoSidedLimit(h::Float64) = TwoSidedLimit([h])

function is_IC(L::TwoSidedLimit, stat::AbstractStatistic)
    #TODO: run tests to check if it is correct
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

