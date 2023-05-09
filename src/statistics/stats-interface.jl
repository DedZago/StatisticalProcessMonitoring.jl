abstract type AbstractStatistic end

get_parameter(::AbstractStatistic) = @NamedTuple{}
export get_parameter
set_parameter!(::AbstractStatistic) = @NamedTuple{}
export set_parameter!
get_value(stat::AbstractStatistic) = stat.value
export get_value
update_statistic!(stat::AbstractStatistic, x) = error("Not implemented for abstract class.")
export update_statistic!
update_statistic(stat::AbstractStatistic, x) = error("Not implemented for abstract class.")
export update_statistic
get_maxrl(::AbstractStatistic) = Inf
get_maxrl(stat::Vector{T}) where T <: AbstractStatistic = minimum(get_maxrl.(stat))
export get_maxrl

include("univariate.jl")