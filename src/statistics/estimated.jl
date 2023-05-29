abstract type EstimatedStatistic <: AbstractStatistic end

################ Control charts with estimated parameter values ################

get_statistic(S::EstimatedStatistic) = S.stat
get_value(S::EstimatedStatistic) = get_value(get_statistic(S))
get_design(S::EstimatedStatistic) = get_design(S.stat)
set_design!(S::EstimatedStatistic, par) = set_design!(S.stat, par)

@with_kw mutable struct LocationScaleEstimatedStatistic{S, F} <: EstimatedStatistic
    stat::S
    μ::F = 0.0
    σ::F = 0.0
end
export LocationScaleEstimatedStatistic
#FIXME: test

update_statistic(S::EstimatedStatistic, x) = update_statistic(get_statistic(S), (x - S.μ)/S.σ)
update_statistic!(S::EstimatedStatistic, x) = update_statistic!(get_statistic(S), (x - S.μ)/S.σ)