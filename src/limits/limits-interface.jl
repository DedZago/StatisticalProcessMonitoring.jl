abstract type AbstractLimit end
get_h(L::AbstractLimit) = L.h
get_h(L::Vector{LIM}) where LIM <: AbstractLimit = get_h.(L)
export get_h
get_value(L::AbstractLimit) = get_h(L)
get_value(L::Vector{LIM}) where LIM <: AbstractLimit = get_h(L)
export get_value
set_h!(L::AbstractLimit, h::Float64) = L.h = h
set_h!(L::AbstractLimit, h::Vector{Float64}) = L.h = h
export set_h!

function is_IC(L::AbstractLimit, stat::AbstractStatistic)
    val = get_value(stat)
    lim = get_value(L)
    return compare_values(lim, val, L)
end
export is_IC

is_OC(L::AbstractLimit, stat::AbstractStatistic) = !is_IC(L, stat)
export is_OC

is_IC_vec(L::Vector{LIM}, stat::Vector{STAT}) where {LIM <: AbstractLimit, STAT <: AbstractStatistic} = is_IC.(L, stat)
export is_IC_vec

export get_curved_value

include("fixed.jl")
include("dynamic.jl")