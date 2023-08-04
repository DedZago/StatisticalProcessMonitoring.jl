abstract type ResidualStatistic <: AbstractStatistic end
export ResidualStatistic

residual!(x, S::ResidualStatistic) = error("Not implemented for abstract class")
export residual!

################ Control charts based on residuals ################

#FIXME: test general interface for residual statistic
get_statistic(S::ResidualStatistic) = S.stat
get_value(S::ResidualStatistic) = get_value(get_statistic(S))
get_design(S::ResidualStatistic) = get_design(get_statistic(S))
set_design!(S::ResidualStatistic, par) = set_design!(get_statistic(S), par)
get_maxrl(S::ResidualStatistic) = get_maxrl(get_statistic(S))

#FIXME: residual! in update_statistic may result in unwanted bugs
update_statistic(S::ResidualStatistic, x) = update_statistic(get_statistic(S), residual!(x, S))
update_statistic!(S::ResidualStatistic, x) = update_statistic!(get_statistic(S), residual!(x, S))

"""
    LocationScaleStatistic{S, M, P}

A mutable struct representing a statistic applied to a location-scale family.

# Fields
- `stat::S`: The statistic.
- `μ::M`: The location parameter.
- `Ω::P`: The precision parameter (inverse of the variance).

# Examples
    STAT = EWMA(λ = 0.2)
    RSTAT = LocationScaleStatistic(STAT, 1.0, 2.5)
"""
@with_kw mutable struct LocationScaleStatistic{S, M, P} <: ResidualStatistic
    stat::S
    μ::M = 0.0
    Ω::P = 0.0
end


export LocationScaleStatistic

#FIXME: test whether the location scale constructors and the residual function work on vector and matrix data
LocationScaleStatistic(stat, x::AbstractVector) = LocationScaleStatistic(stat, mean(x), 1.0/std(x))
LocationScaleStatistic(stat, x::AbstractMatrix) = LocationScaleStatistic(stat, mean.(eachcol(x)), inv(sqrt(cov(x))))
residual!(x, S::LocationScaleStatistic) = S.Ω * (x .- S.μ)