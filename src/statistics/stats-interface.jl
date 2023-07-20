abstract type AbstractStatistic end

get_design(::AbstractStatistic) = @NamedTuple{}
export get_design
set_design!(::AbstractStatistic) = error("Not implemented for abstract class.") 
export set_design!

"""
    get_value(stat::AbstractStatistic)

Get the current value of a statistic.
"""
get_value(stat::AbstractStatistic) = stat.value
export get_value

"""
    function set_value!(stat::AbstractStatistic, value)

Set the value of a statistic.
"""
function set_value!(stat::AbstractStatistic, value)
    stat.value = value
end
export set_value!

"""
    update_statistic!(stat::AbstractStatistic, x)

Update a statistic with a new observation x
"""
update_statistic(stat::AbstractStatistic, x) = error("Not implemented for abstract class.")
export update_statistic
update_statistic!(stat::AbstractStatistic, x) = error("Not implemented for abstract class.")
export update_statistic!
get_maxrl(::AbstractStatistic) = Inf
get_maxrl(stat::Vector{T}) where T <: AbstractStatistic = minimum(get_maxrl.(stat))
export get_maxrl

include("univariate.jl")
include("residual.jl")