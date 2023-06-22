abstract type AbstractStatistic end

get_design(::AbstractStatistic) = @NamedTuple{}
export get_design
set_design!(::AbstractStatistic) = error("Not implemented for abstract class.") 
export set_design!
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
include("estimated.jl")