using SPM
using Parameters

abstract type AbstractChart end

"""
    get_limit(CC::AbstractChart)

Get the control limit of a control chart.
"""
get_limit(CC::AbstractChart) = CC.limit
export get_limit

"""
    get_limit_value(CC::AbstractChart)

Get the control limit value of a control chart.
"""
get_limit_value(CC::AbstractChart) = get_limit_value(get_limit(CC))
export get_limit_value

"""
    get_statistic(CC::AbstractChart)

Get the statistic of a control chart.
"""
get_statistic(CC::AbstractChart) = CC.stat
export get_statistic

"""
    get_value(CC::AbstractChart)
    
Get the current value of the control chart statistic.
"""
get_value(CC::AbstractChart) = get_value(get_statistic(CC))
export get_value

"""
    get_nominal(CC::AbstractChart)

Get the nominal properties of a control chart.
"""
get_nominal(CC::AbstractChart) = CC.nominal
export get_nominal

"""
    get_phase1(CC::AbstractChart)

Get the Phase 1 information of a control chart.
"""
get_phase1(CC::AbstractChart) = CC.phase1
export get_phase1

"""
    get_parameters(CC::AbstractChart)
    set_parameters(CC::AbstractChart, par)
    
Get and set the parameters of the control chart statistic.
"""
get_param(CC::AbstractChart) = @NamedTuple{}
set_param(CC::AbstractChart, par) = set_param(get_statistic(CC), par)
export get_param

"""
    get_value(CC::AbstractChart)
    
Get the maximum run length of the control chart.
"""
get_maxrl(CC::AbstractChart) = get_maxrl(get_statistic(CC))
export get_maxrl

"""
    update_chart!(CC::AbstractChart, x)
    
Update the control chart using a new observation `x`.
"""
update_chart!(CC::AbstractChart, x) = update_statistic!(get_statistic(CC), x)
export update_chart!

"""
    is_IC(CC::AbstractChart)
    is_OC(CC::AbstractChart)
    
Check whether the control chart is in control or out of control.
"""
is_IC(CC::AbstractChart) = is_IC(get_limit(CC), get_value(CC))
is_OC(CC::AbstractChart) = !is_IC(CC)
export is_IC
export is_OC


#################################################################
#               Generic control chart interface                 #
#################################################################
@with_kw mutable struct ControlChart{STAT <: AbstractStatistic, LIM <: AbstractLimit, NOM <: AbstractNominal, PH1 <: AbstractPhase1} <: AbstractChart
    stat::STAT
    limit::LIM
    nominal::NOM
    phase1::PH1
end

"""
    function set_statistic!(CC::AbstractChart, statistic::AbstractStatistic)

Set the statistic of a control chart.
"""
function set_statistic!(CC::C, statistic::STAT) where C <: AbstractChart STAT <: AbstractStatistic
    CC.stat = statistic
    return statistic
end
export set_statistic!


"""
    function set_limit!(CC::AbstractChart, limit::AbstractLimit)

Set the control limit of a control chart.
"""
function set_limit!(CC::C, limit::LIM) where C <: AbstractChart LIM <: AbstractLimit
    CC.limit = limit
    return limit
end
function set_limit!(CC::C, limit::Vector{Float64}) where C <: AbstractChart 
    set_limit!(get_limit(CC), limit)
    return get_limit(CC)
end
export set_limit! 


"""
    function set_parameter!(CC::C, par)

Set the parameters of the control chart statistic.
"""
function set_parameter!(CC::C, par) where C <: AbstractChart
    set_parameter!(get_statistic(CC), par)
    return par
end
export set_parameter!

"""
    set_phase1!(CC::AbstractChart, phase1::AbstractPhase1)

Set the Phase 1 information of a control chart.
"""
function set_phase1!(CC::C, phase1::PH1) where C <: AbstractChart PH1 <: AbstractPhase1
    CC.phase1 = phase1
    return phase1
end
export set_phase1!

"""
    set_nominal!(CC::AbstractChart, nominal::AbstractNominal)

Set the nominal properties of a control chart.
"""
function set_nominal!(CC::C, nominal::N) where C <: AbstractChart N <: AbstractNominal
    CC.nominal = nominal
    return nominal
end
export set_nominal!