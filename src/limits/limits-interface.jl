abstract type AbstractLimit end
get_value(L::AbstractLimit) = L.value
export get_value
set_value!(L::AbstractLimit, h::Float64) = L.value = h
export set_value!

function is_IC(L::AbstractLimit, stat::AbstractStatistic)
    val = get_value(stat)
    lim = get_value(L)
    @assert typeof(val) == typeof(lim)
    return compare_values(lim, val, L)
end
export is_IC
is_OC(L::AbstractLimit, stat::AbstractStatistic) = !is_IC(L, stat)
export is_OC

export get_curved_value

include("fixed.jl")
include("dynamic.jl")