abstract type AbstractStatistic end

get_parameters(::AbstractStatistic) = @NamedTuple{}
export get_parameters
get_value(stat::AbstractStatistic) = stat.value
export get_value
update_statistic!(stat::AbstractStatistic, x) = stat
export update_statistic!
get_maxrl(::AbstractStatistic) = Inf
export get_maxrl