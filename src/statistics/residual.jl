abstract type ResidualStatistic <: AbstractStatistic end
export ResidualStatistic

residual!(x, S::ResidualStatistic) = error("Not implemented for abstract class")
export residual!

################ Control charts based on residuals ################

get_statistic(S::ResidualStatistic) = S.stat
set_value!(S::ResidualStatistic, x) = set_value!(get_statistic(S), x)
get_value(S::ResidualStatistic) = get_value(get_statistic(S))
get_design(S::ResidualStatistic) = get_design(get_statistic(S))
set_design!(S::ResidualStatistic, par::Real) = set_design!(get_statistic(S), par)
set_design!(S::ResidualStatistic, par::AbstractVector) = set_design!(get_statistic(S), par)
get_maxrl(S::ResidualStatistic) = get_maxrl(get_statistic(S))

update_statistic(S::ResidualStatistic, x) = update_statistic(get_statistic(S), residual!(x, deepcopy(S)))
update_statistic!(S::ResidualStatistic, x) = update_statistic!(get_statistic(S), residual!(x, S))

"""
    LocationScaleStatistic{S, M, P}

A mutable struct representing a statistic applied to a location-scale family.

# Fields
- `stat::S`: The statistic.
- `μ::M`: The location parameter.
- `Ω::P`: The inverse square root of the variance.

# Examples
    STAT = EWMA(λ = 0.2)
    RSTAT = LocationScaleStatistic(STAT, 1.0, 2.5)
"""
@with_kw mutable struct LocationScaleStatistic{S, M, P} <: ResidualStatistic
    stat::S
    μ::M = 0.0
    Ω::P = 1.0
end
export LocationScaleStatistic

LocationScaleStatistic(stat, x::AbstractVector) = LocationScaleStatistic(stat, mean(x), 1.0/std(x))
LocationScaleStatistic(stat, x::AbstractMatrix) = LocationScaleStatistic(stat, mean.(eachcol(x)), inv(sqrt(cov(x))))
residual!(x, S::LocationScaleStatistic) = S.Ω * (x .- S.μ)

