abstract type AbstractStatistic end

get_param(::AbstractStatistic) = @NamedTuple{}
export get_param
set_param(::AbstractStatistic) = @NamedTuple{}
export set_param
get_value(stat::AbstractStatistic) = stat.value
export get_value
update_statistic!(stat::AbstractStatistic, x) = stat
export update_statistic!
get_maxrl(::AbstractStatistic) = Inf
export get_maxrl

include("univariate.jl")