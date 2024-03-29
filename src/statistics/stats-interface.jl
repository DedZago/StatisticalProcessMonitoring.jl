abstract type AbstractStatistic end

"""
    get_design(stat::AbstractStatistic)

Get the vector of hyperparameters of a statistic.
"""
get_design(::AbstractStatistic) = Vector{Float64}()
export get_design

"""
    set_design!(stat::AbstractStatistic, par::AbstractVector)

Set the vector of hyperparameters of a statistic.
"""
set_design!(::AbstractStatistic, ::AbstractVector) = error("Not implemented for abstract class.") 
export set_design!

"""
    get_value(stat::AbstractStatistic)

Get the current value of a statistic.
"""
get_value(stat::AbstractStatistic) = stat.value
get_value(S::NTuple{N, AbstractStatistic}) where N = collect(get_value.(S))
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
update_statistic!(stat::AbstractStatistic, x) = error("Not implemented for abstract class.")
export update_statistic!
update_statistic(stat::AbstractStatistic, x) = update_statistic!(deepcopy(stat), x)
export update_statistic

"""
    get_maxrl(stat::AbstractStatistic)
    get_maxrl(stat::Vector{T <: AbstractStatistics})

Get the maximum value of the run length for a statistic `stat`.
"""
get_maxrl(::AbstractStatistic) = Inf
get_maxrl(stat::Vector{T}) where T <: AbstractStatistic = minimum(get_maxrl.(stat))
get_maxrl(stat::Tuple) = minimum(get_maxrl.(stat))
export get_maxrl

include("univariate.jl")
include("multivariate.jl")
include("residual.jl")
include("categorization/categorization.jl")
include("categorization/LLCUSUM.jl")
include("categorization/LLD.jl")
include("categorization/MOC.jl")
include("functional/functional.jl")
include("functional/NEWMA.jl")
include("partially-observed/adaptive-sampling.jl")
include("partially-observed/RSADA.jl")
include("partially-observed/TRAS.jl")