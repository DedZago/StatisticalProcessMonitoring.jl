abstract type AbstractLimit end
get_limit_value(L::AbstractLimit) = L.value
is_IC(L::AbstractLimit, stat::AbstractStatistic) = error("Not implemented for abstract interface.")
is_OC(L::AbstractLimit, stat::AbstractStatistic) = !is_IC(L, stat)

include("fixed.jl")