using SPM
using Parameters

abstract type AbstractChart{STAT, LIM, NOM, PH2} end

#################################################################
#               Generic control chart interface                 #
#################################################################
mutable struct ControlChart{STAT, LIM, NOM, PH2} <: AbstractChart{STAT, LIM, NOM, PH2}
    stat::STAT
    limit::LIM
    nominal::NOM
    phase2::PH2
    t::Int

    ControlChart(stat::S, limit::L, nominal::N, phase2::P, t::Int) where {S <: AbstractStatistic, L <: AbstractLimit, N <: NominalProperties, P <: AbstractPhase2} = new{S,L,N,P}(deepcopy(stat), deepcopy(limit), deepcopy(nominal), deepcopy(phase2), t)

    ControlChart(stat::S, limit::L, nominal::N, phase2::P) where {S <: AbstractStatistic, L <: AbstractLimit, N <: NominalProperties, P <: AbstractPhase2} = new{S,L,N,P}(deepcopy(stat), deepcopy(limit), deepcopy(nominal), deepcopy(phase2), 0)

    ControlChart(stat::Vector{S}, limit::Vector{L}, nominal::N, phase2::P) where {S <: AbstractStatistic, L <: AbstractLimit, N <: NominalProperties, P <: AbstractPhase2} = new{Vector{S},Vector{L},N,P}(deepcopy(stat), deepcopy(limit), deepcopy(nominal), deepcopy(phase2), 0)

    ControlChart(stat::Vector{S}, limit::Vector{L}, nominal::N, phase2::P, t::Int) where {S <: AbstractStatistic, L <: AbstractLimit, N <: NominalProperties, P <: AbstractPhase2} = new{Vector{S},Vector{L},N,P}(deepcopy(stat), deepcopy(limit), deepcopy(nominal), deepcopy(phase2), t)

    ControlChart(stat::S, limit::L, nominal::N, phase2::P, t::Int) where {S <: Tuple, L <: Tuple, N <: NominalProperties, P <: AbstractPhase2} = new{S,L,N,P}(deepcopy(stat), deepcopy(limit), deepcopy(nominal), deepcopy(phase2), t)

    ControlChart(stat::S, limit::L, nominal::N, phase2::P) where {S <: Tuple, L <: Tuple, N <: NominalProperties, P <: AbstractPhase2} = new{S,L,N,P}(deepcopy(stat), deepcopy(limit), deepcopy(nominal), deepcopy(phase2), 0)
end
export ControlChart

Base.show(io::IO, CH::ControlChart) = print(io, "stat: $(get_statistic(CH))\nlimit: $(get_limit(CH))\nnominal: $(get_nominal(CH))\nphase2: $(typeof(get_phase2(CH)))\n\nt: $(get_t(CH))")

const MultipleControlChart{S,L,N,P} = ControlChart{S,L,N,P} where {S<:Tuple,L<:Tuple,N,P}
export MultipleControlChart


"""
    shallow_copy_sim(CH::ControlChart)

Create a shallow copy of a control chart, so that only the statistic and the control limit are copied.
This is done to prevent copying eventual Phase 2 data multiple times and thus reduce computational effort when optimizing the control limit and the chart tuning parameters.
"""
shallow_copy_sim(CH::ControlChart) = ControlChart(deepcopy(get_statistic(CH)), deepcopy(get_limit(CH)), deepcopy(get_nominal(CH)), shallow_copy_sim(get_phase2(CH)), deepcopy(get_t(CH)))
export shallow_copy_sim


"""
    get_limit(CH::AbstractChart)

Get the control limit of a control chart.
"""
get_limit(CH::AbstractChart) = CH.limit
export get_limit


"""
    get_limit_value(CH::AbstractChart)

Get the control limit value of a control chart.
"""
get_limit_value(CH::AbstractChart) = get_value(get_limit(CH))

get_limit_value(CH::MultipleControlChart) = collect(get_value.(get_limit(CH)))

# get_limit_value(CH::AbstractChart{STAT,LIM,NOM,PH2}) where {STAT, LIM <: OneSidedCurvedLimit, NOM, PH2} = get_value(get_limit(CH), get_t(CH), get_statistic(CH))

# get_limit_value(CH::AbstractChart{STAT,LIM,NOM,PH2}) where {STAT, LIM <: TwoSidedCurvedLimit, NOM, PH2} = get_value(get_limit(CH), get_t(CH), get_statistic(CH))
export get_limit_value


"""
    get_statistic(CH::AbstractChart)

Get the statistic of a control chart.
"""
get_statistic(CH::AbstractChart) = CH.stat
export get_statistic


"""
    get_value(CH::AbstractChart)
    
Get the current value of the control chart statistic.
"""
get_value(CH::AbstractChart) = get_value(get_statistic(CH))
get_value(CH::MultipleControlChart) = collect(get_value.(get_statistic(CH)))
export get_value


"""
    get_nominal(CH::AbstractChart)

Get the nominal properties of a control chart.
"""
get_nominal(CH::AbstractChart) = CH.nominal
export get_nominal

"""
    get_nominal(CH::AbstractChart)

Get the nominal value of a control chart.
"""
get_nominal_value(CH::AbstractChart) = get_value(get_nominal(CH))
export get_nominal_value


"""
    get_phase2(CH::AbstractChart)

Get the phase 2 information of a control chart.
"""
get_phase2(CH::AbstractChart) = CH.phase2
export get_phase2

"""
    get_t(CH::AbstractChart)

Get the current time point of a control chart.
"""
get_t(CH::AbstractChart) = CH.t
export get_t



"""
    get_design(CH::AbstractChart)
    
Get the designs of the control chart statistic.
"""
get_design(CH::AbstractChart) = get_design(get_statistic(CH))

get_design(CH::MultipleControlChart) = collect(get_design.(get_statistic(CH)))
export get_design


"""
    get_design(CH::AbstractChart)
    
Set the designs of the control chart statistic.
"""
set_design!(CH::AbstractChart, par) = set_design!(get_statistic(CH), par)

set_design!(CH::MultipleControlChart, par) = set_design!.(get_statistic(CH), par)

function set_design!(CH::MultipleControlChart, par::Vector{T}) where T
    @assert length(get_statistic(CH)) == length(par)
    for i in 1:length(get_statistic(CH))
        set_design!(get_statistic(CH)[i], par[i])
    end
end
export set_design!



"""
    get_maxrl(CH::AbstractChart)
    
Get the maximum run length of the control chart.
"""
get_maxrl(CH::AbstractChart) = get_maxrl(get_statistic(CH))
export get_maxrl


"""
    is_IC(CH::AbstractChart)
    is_OC(CH::AbstractChart)
    
Check whether the control chart is in control or out of control.
"""
is_IC(CH::AbstractChart) = is_IC(get_limit(CH), get_statistic(CH))

is_IC(CH::MultipleControlChart) = all(is_IC_vec(get_limit(CH), get_statistic(CH)))

is_IC(CH::AbstractChart{STAT,LIM,NOM,PH2}) where {STAT, LIM <: OneSidedCurvedLimit, NOM, PH2} = is_IC(get_limit(CH), get_t(CH), get_statistic(CH))

is_IC(CH::AbstractChart{STAT,LIM,NOM,PH2}) where {STAT, LIM <: TwoSidedCurvedLimit, NOM, PH2} = is_IC(get_limit(CH), get_t(CH), get_statistic(CH))

is_OC(CH::AbstractChart) = !is_IC(CH)
export is_IC
export is_OC

"""
    is_IC_vec(CH::MultipleControlChart)
    is_OC_vec(CH::MultipleControlChart)

Check whether each individual control chart that makes up a multiple control chart is in control or out of control.

### Returns
A vector of Bool, whose length is equal to the number of individual statistics.
"""
is_IC_vec(CH::MultipleControlChart) = is_IC_vec(get_limit(CH), get_statistic(CH))
is_OC_vec(CH::MultipleControlChart) = .!(is_IC_vec(get_limit(CH), get_statistic(CH)))
export is_IC_vec
export is_OC_vec

"""
    function set_statistic!(CH::AbstractChart, statistic::AbstractStatistic)

Set the statistic of a control chart.

### Returns
The new value of the statistic.
"""
function set_statistic!(CH::C, statistic::STAT) where C <: AbstractChart where STAT <: AbstractStatistic
    CH.stat = statistic
    return statistic
end
export set_statistic!


"""
    function set_t!(CH::AbstractChart, t)

Set the current time value of a control chart.

### Returns
The new time value of the control chart.
"""
function set_t!(CH::C, t::Int) where C <: AbstractChart
    CH.t = t
end
export set_t!


"""
    function set_value!(CH::AbstractChart, value)

Set the value of the statistic of a control chart.

### Returns
The new value of the control chart's statistic.
"""
function set_value!(CH::C, value) where C <: AbstractChart
    set_value!(get_statistic(CH), value)
end
function set_value!(CH::MultipleControlChart, value::AbstractVector)
    set_value!.(get_statistic(CH), value)
end
function set_value!(CH::MultipleControlChart, value::Real, j::Int)
    set_value!(get_statistic(CH)[j], value)
end
export set_value!


"""
    function set_limit!(CH::AbstractChart, limit::AbstractLimit)
    function set_limit!(CH::AbstractChart, h::Float64)
    function set_limit!(CH::MultipleControlChart, h::Vector{Float64})
    function set_limit!(CH::MultipleControlChart, h::Float64)
    function set_limit!(CH::MultipleControlChart, h::Float64, j::Int)

Set the control limit of a control chart.

### Returns
The new control limit.
"""
function set_limit!(CH::AbstractChart, limit::AbstractLimit)
    CH.limit = limit
    return limit
end

function set_limit!(CH::AbstractChart, h::Float64)
    set_h!(get_limit(CH), h)
    return get_limit(CH)
end
export set_limit! 

set_limit!(CH::MultipleControlChart, h::Float64) = set_h!.(get_limit(CH), h)
set_limit!(CH::MultipleControlChart, h::Float64, j::Int) = set_h!(get_limit(CH)[j], h)

function set_limit!(CH::MultipleControlChart, h::Vector{Float64})
    @assert length(get_limit(CH)) == length(h)
    for i in 1:length(get_limit(CH))
        set_h!(get_limit(CH)[i], h[i])
    end
end

"""
    set_phase2!(CH::AbstractChart, phase2::AbstractPhase2)

Set the Phase 2 information of a control chart to simulate run lenghts.
"""
function set_phase2!(CH::C, phase2::PH2) where C <: AbstractChart where PH2 <: AbstractPhase2
    CH.phase2 = phase2
    return phase2
end
export set_phase2!


"""
    new_data(CH::AbstractChart)

Simulate a new observation based on the control chart's Phase II object.
"""
new_data(CH::AbstractChart) = new_data(get_phase2(CH))
export new_data

"""
    new_data(CH::AbstractChart)

Simulate a new observation for the control chart from the phase 2 data, eventually modifying the underlying phase 2 object.
"""
new_data!(CH::AbstractChart) = new_data!(get_phase2(CH))
export new_data
    

"""
    set_nominal!(CH::AbstractChart, nominal::NominalProperties)

Set the nominal properties of a control chart.
"""
function set_nominal!(CH::C, nominal::N) where C <: AbstractChart where N <: NominalProperties
    CH.nominal = nominal
    return nominal
end
export set_nominal!


"""
    update_chart(CH::AbstractChart, x)
    
Update the control chart without modifying it using a new observation `x`.
"""
update_chart(CH::AbstractChart, x) = update_chart!(shallow_copy_sim(CH), x)
export update_chart


"""
    update_chart!(CH::AbstractChart, x)
    
Update the control chart inplace using a new observation `x`.
"""
function update_chart!(CH::AbstractChart, x)
    CH.t += 1
    update_statistic!(get_statistic(CH), x)
end

function update_chart!(CH::MultipleControlChart, x)
    CH.t += 1
    for i in eachindex(get_statistic(CH))
        update_statistic!(get_statistic(CH)[i], x)
    end
end

function update_chart!(CH::AbstractChart{STAT, LIM, NOM, PH2}, x) where {STAT, LIM <: DynamicLimit, NOM, PH2}
    CH.t += 1
    update_limit!(CH)
    update_statistic!(get_statistic(CH), x)
end
export update_chart!


"""
    update_limit!(CH::AbstractChart, x)
    
Update the dynamic control limit of a control chart inplace.
"""
function update_limit!(CH::AbstractChart{S,L,N,P1}) where {S, L <: DynamicLimit, N, P1}
    update_value!(get_limit(CH))
end

function update_limit!(CH::AbstractChart{S,L,N,P1}) where {S, L <: BootstrapLimit, N, P1}
    for b in 1:length(get_limit(CH).sim)
        get_limit(CH).sim[b] = update_statistic(get_statistic(CH), new_data(CH))
    end
    update_value!(get_limit(CH), get_nominal(CH))
    resample_sims!(get_limit(CH))
end
export update_limit!

include("simulate.jl")