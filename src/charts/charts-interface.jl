using SPM
using Parameters

abstract type AbstractChart{STAT, LIM, NOM, PH1} end

#################################################################
#               Generic control chart interface                 #
#################################################################
mutable struct ControlChart{STAT, LIM, NOM, PH1} <: AbstractChart{STAT, LIM, NOM, PH1}
    stat::STAT
    limit::LIM
    nominal::NOM
    phase1::PH1
    t::Int

    ControlChart(stat::S, limit::L, nominal::N, phase1::P, t::Int) where {S <: AbstractStatistic, L <: AbstractLimit, N <: NominalProperties, P <: AbstractPhase1} = new{S,L,N,P}(stat, limit, nominal, phase1, t)
    ControlChart(stat::S, limit::L, nominal::N, phase1::P) where {S <: AbstractStatistic, L <: AbstractLimit, N <: NominalProperties, P <: AbstractPhase1} = new{S,L,N,P}(stat, limit, nominal, phase1, 0)
    ControlChart(stat::Vector{S}, limit::Vector{L}, nominal::N, phase1::P) where {S <: AbstractStatistic, L <: AbstractLimit, N <: NominalProperties, P <: AbstractPhase1} = new{Vector{S},Vector{L},N,P}(stat, limit, nominal, phase1, 0)
    ControlChart(stat::Vector{S}, limit::Vector{L}, nominal::N, phase1::P, t::Int) where {S <: AbstractStatistic, L <: AbstractLimit, N <: NominalProperties, P <: AbstractPhase1} = new{Vector{S},Vector{L},N,P}(stat, limit, nominal, phase1, t)
end
export ControlChart

Base.show(io::IO, CH::ControlChart) = print(io, "stat: $(get_statistic(CH))\nlimit: $(get_limit(CH))\nnominal: $(get_nominal(CH))\nphase1: $(typeof(get_phase1(CH)))\n\nt: $(get_t(CH))")

const MultipleControlChart{S,L,N,P} = ControlChart{Vector{S}, Vector{L},N,P} where {S,L,N,P}
export MultipleControlChart

shallow_copy_sim(CH::ControlChart) = ControlChart(deepcopy(get_statistic(CH)), deepcopy(get_limit(CH)), get_nominal(CH), get_phase1(CH), get_t(CH))
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

get_limit_value(CH::MultipleControlChart) = get_value.(get_limit(CH))

get_limit_value(CH::AbstractChart{STAT,LIM,NOM,PH1}) where {STAT, LIM <: OneSidedCurvedLimit, NOM, PH1} = get_curved_value(get_limit(CH), get_t(CH) + 1, get_statistic(CH))

get_limit_value(CH::AbstractChart{STAT,LIM,NOM,PH1}) where {STAT, LIM <: TwoSidedCurvedLimit, NOM, PH1} = get_curved_value(get_limit(CH), get_t(CH) + 1, get_statistic(CH))

get_limit_value(CH::MultipleControlChart{STAT,LIM,NOM,PH1}) where {STAT, LIM <: OneSidedCurvedLimit, NOM, PH1} = get_curved_value.(get_limit(CH), get_t(CH) + 1, get_statistic(CH))

get_limit_value(CH::MultipleControlChart{STAT,LIM,NOM,PH1}) where {STAT, LIM <: TwoSidedCurvedLimit, NOM, PH1} = get_curved_value.(get_limit(CH), get_t(CH) + 1, get_statistic(CH))
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
get_value(CH::MultipleControlChart) = get_value.(get_statistic(CH))
export get_value


"""
    get_nominal(CH::AbstractChart)
    get_nominal_value(CH::AbstractChart)

Get the nominal properties of a control chart.
"""
get_nominal(CH::AbstractChart) = CH.nominal
export get_nominal
get_nominal_value(CH::AbstractChart) = get_value(get_nominal(CH))
export get_nominal_value


"""
    get_phase1(CH::AbstractChart)

Get the Phase 1 information of a control chart.
"""
get_phase1(CH::AbstractChart) = CH.phase1
export get_phase1

"""
    get_t(CH::AbstractChart)

Get the current time point from a control chart.
"""
get_t(CH::AbstractChart) = CH.t
export get_t



"""
    get_parameter(CH::AbstractChart)
    set_parameter!(CH::AbstractChart, par)
    
Get and set the parameters of the control chart statistic.
"""
get_parameter(CH::AbstractChart) = get_parameter(get_statistic(CH))

get_parameter(CH::MultipleControlChart) = get_parameter.(get_statistic(CH))
export get_parameter


set_parameter!(CH::AbstractChart, par) = set_parameter!(get_statistic(CH), par)

set_parameter!(CH::MultipleControlChart, par) = set_parameter!.(get_statistic(CH), par)

function set_parameter!(CH::MultipleControlChart, par::AbstractVector)
    @assert length(get_statistic(CH)) == length(par)
    for i in 1:length(get_statistic(CH))
        set_parameter!(get_statistic(CH)[i], par[i])
    end
end
export set_parameter!



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

is_IC(CH::AbstractChart{STAT,LIM,NOM,PH1}) where {STAT, LIM <: OneSidedCurvedLimit, NOM, PH1} = is_IC(get_limit(CH), get_t(CH) + 1, get_statistic(CH))

is_IC(CH::AbstractChart{STAT,LIM,NOM,PH1}) where {STAT, LIM <: TwoSidedCurvedLimit, NOM, PH1} = is_IC(get_limit(CH), get_t(CH) + 1, get_statistic(CH))

is_OC(CH::AbstractChart) = !is_IC(CH)
export is_IC
export is_OC

"""
    is_IC_vec(CH::MultipleControlChart)
    is_OC_vec(CH::MultipleControlChart)


"""
is_IC_vec(CH::MultipleControlChart) = is_IC_vec(get_limit(CH), get_statistic(CH))
is_OC_vec(CH::MultipleControlChart) = .!(is_IC_vec(get_limit(CH), get_statistic(CH)))
export is_IC_vec
export is_OC_vec

"""
    function set_statistic!(CH::AbstractChart, statistic::AbstractStatistic)

Set the statistic of a control chart.
"""
function set_statistic!(CH::C, statistic::STAT) where C <: AbstractChart where STAT <: AbstractStatistic
    CH.stat = statistic
    return statistic
end
export set_statistic!


"""
    function set_limit!(CH::AbstractChart, limit::AbstractLimit)

Set the control limit of a control chart.
"""
function set_limit!(CH::AbstractChart, limit::AbstractLimit)
    CH.limit = limit
    return limit
end

function set_limit!(CH::AbstractChart, limit::Float64)
    set_value!(get_limit(CH), limit)
    return get_limit(CH)
end
export set_limit! 

set_limit!(CH::MultipleControlChart, limit::Float64) = set_value!.(get_limit(CH), limit)

function set_limit!(CH::MultipleControlChart, limit::Vector{Float64})
    @assert length(get_limit(CH)) == length(limit)
    for i in 1:length(get_limit(CH))
        set_value!(get_limit(CH)[i], limit[i])
    end
end

"""
    set_phase1!(CH::AbstractChart, phase1::AbstractPhase1)

Set the Phase 1 information of a control chart.
"""
function set_phase1!(CH::C, phase1::PH1) where C <: AbstractChart where PH1 <: AbstractPhase1
    CH.phase1 = phase1
    return phase1
end
export set_phase1!

"""
    new_data(CH::AbstractChart)

Simulate a new observation for the control chart.
"""
new_data(CH::AbstractChart) = new_data(get_phase1(CH))
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
    update_chart!(CH::AbstractChart, x)
    
Update the control chart using a new observation `x`.
"""
function update_chart!(CH::AbstractChart, x)
    CH.t += 1
    update_statistic!(get_statistic(CH), x)
end

function update_chart!(CH::MultipleControlChart, x)
    CH.t += 1
    update_statistic!.(get_statistic(CH), x)
end

function update_chart!(CH::AbstractChart{STAT, LIM, NOM, PH1}, x) where {STAT, LIM <: DynamicLimit, NOM, PH1}
    CH.t += 1
    update_limit!(CH)
    update_statistic!(get_statistic(CH), x)
end

export update_chart!

update_chart(CH::AbstractChart, x) = update_chart!(shallow_copy_sim(CH), x)
export update_chart

function update_limit!(CH::AbstractChart{S,L,N,P1}) where {S, L <: BootstrapLimit, N, P1}
    for b in 1:length(get_limit(CH).sim)
        get_limit(CH).sim[b] = update_statistic!(deepcopy(get_statistic(CH)), new_data(CH))
    end
    update_value!(get_limit(CH), get_nominal(CH))
    resample_sims!(get_limit(CH))
end
export update_limit!#FIXME: tests

#FIXME: organize functions more neatly

# get_limit_value(CH::AbstractChart{STAT,LIM,NOM,PH1}) where {STAT, LIM <: OneSidedCurvedLimit, NOM, PH1} = get_curved_value(get_limit(CH), get_t(CH) + 1, get_statistic(CH))


include("simulate.jl")