abstract type AbstractLimit end
get_value(L::AbstractLimit) = L.value
is_IC(L::AbstractLimit, stat::AbstractStatistic) = error("Not implemented for abstract interface.")
export is_IC
is_OC(L::AbstractLimit, stat::AbstractStatistic) = !is_IC(L, stat)
export is_OC

include("fixed.jl")