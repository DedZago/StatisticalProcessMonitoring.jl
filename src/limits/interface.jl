abstract type AbstractLimit end
get_value(L::AbstractLimit) = L.value
export get_value
set_value!(L::AbstractLimit, h::AbstractVector) = L.value .= h
set_value!(L::AbstractLimit, h::Float64) = L.value .= h
export set_value!
is_IC(L::AbstractLimit, stat::AbstractStatistic) = error("Not implemented for abstract interface.")
export is_IC
is_OC(L::AbstractLimit, stat::AbstractStatistic) = !is_IC(L, stat)
export is_OC

include("fixed.jl")