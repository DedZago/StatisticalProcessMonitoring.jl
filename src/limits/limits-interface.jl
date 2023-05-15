abstract type AbstractLimit{T} end
get_h(L::AbstractLimit) = L.h

function get_h(L::Vector{AbstractLimit{T}}) where T
    ret = Vector{T}(undef, length(L))
    for i in 1:length(L)
        ret[i] = get_h(L[i])
    end
    return ret
end
export get_h

get_value(L::AbstractLimit) = get_h(L)
get_value(L::Vector{AbstractLimit{T}}) where T = get_h(L)
export get_value
set_h!(L::AbstractLimit{T}, h::T) where T = L.h = h
set_h!(L::AbstractLimit{T}, h::Vector{T}) where T = L.h = h
export set_h!

function is_IC(L::AbstractLimit, stat::AbstractStatistic)
    val = get_value(stat)
    lim = get_value(L)
    return compare_values(lim, val, L)
end
export is_IC

is_OC(L::AbstractLimit, stat::AbstractStatistic) = !is_IC(L, stat)
export is_OC

is_IC_vec(L::Vector{AbstractLimit{T}}, stat::Vector{STAT}) where {T, STAT <: AbstractStatistic} = is_IC.(L, stat)
export is_IC_vec

export get_curved_value

include("fixed.jl")
include("dynamic.jl")