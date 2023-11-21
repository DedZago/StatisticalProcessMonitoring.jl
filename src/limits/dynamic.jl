using StatsBase

abstract type DynamicLimit{T} <: AbstractLimit{T} end
update_value!(L::DynamicLimit) = L.t = L.t + 1

function is_IC(L::AbstractLimit, t, stat::AbstractStatistic)
    val = get_value(stat)
    lim = get_value(L, t)
    return compare_values(lim, val, L)
end


"""
    OneSidedCurvedLimit(h::Float64, upw::Bool)
    OneSidedCurvedLimit(h::Vector{T}, upw::Vector{Bool})

Curved one-sided limit, such that the run length ``RL`` of a control chart is the first time ``t`` in which the statistic ``C_t`` crosses the limit.

* if `upw == true`, ``RL = \\inf\\{t : C_t > h\\cdot f(t)\\}``
* if `upw == false`, ``RL = \\inf\\{t : C_t < -h\\cdot f(t)\\}``

Note that by definition, `h > 0`.
"""
@with_kw mutable struct OneSidedCurvedLimit{T, F <: Function} <: DynamicLimit{T}
    h::T
    fun::F
    upw::Bool = true
    t::Int = 0

    @assert h > 0.0
end
export OneSidedCurvedLimit

OneSidedCurvedLimit(h, f, upw) = OneSidedCurvedLimit(h, f, upw, 0)
get_value(L::OneSidedCurvedLimit, t) = (2.0*float(L.upw) - 1.0) * get_h(L) * L.fun(t)
get_value(L::OneSidedCurvedLimit) = get_value(L, L.t)

function compare_values(lim_val, stat_val, L::OneSidedCurvedLimit)
    for i in 1:length(lim_val)
        if L.upw[i]
            stat_val[i] <= lim_val[i] || return false
        else
            stat_val[i] >= lim_val[i] || return false
        end
    end
    return true
end
 



#TODO: Fix documentation
"""
    TwoSidedCurvedLimit(h::Float64)
    TwoSidedCurvedLimit(h::Vector{T})

Curved one-sided limit, such that the run length ``RL`` of a control chart is the first time ``t`` in which the statistic ``C_t`` crosses the limit.

``RL = \\inf\\{t > 0 : |C_t| > h\\cdot f(t)\\}``.

Note that by definition, `h > 0`.
"""
@with_kw mutable struct TwoSidedCurvedLimit{T, F <: Function} <: DynamicLimit{T}
    h::T
    fun::F
    t::Int = 0

    @assert h > 0.0
end
export TwoSidedCurvedLimit

TwoSidedCurvedLimit(h, f) = TwoSidedCurvedLimit(h, f, 0)
get_value(L::TwoSidedCurvedLimit, t) = [-get_h(L) * L.fun(t), get_h(L) * L.fun(t)]
get_value(L::TwoSidedCurvedLimit) = get_value(L, L.t)

function compare_values(lim_val, stat_val, L::TwoSidedCurvedLimit)
    if (stat_val < lim_val[1]) || (stat_val > lim_val[2])
        return false
    end
    return true
end


abstract type BootstrapLimit{T} <: DynamicLimit{T} end

Base.show(io::IO, L::BootstrapLimit) = print(io, "$(typeof(L))\n  value: $(get_value(L))\n  Bootstrap samples: $(length(L.sim))\n")

"""
    OneSidedBootstrapLimit{T} <: BootstrapLimit{T}

A one-sided bootstrap limit with constant false-alarm rate.

# Fields
- `sim::Vector{T}`: The vector of simulated statistics.
- `h::T`: The current value of the control limit.
- `upw::Bool`: Whether the control limit is an upper or lower control limit.

# Constructors
- `OneSidedBootstrapLimit(S::AbstractStatistic, upw, B::Int)`: Create a new `OneSidedBootstrapLimit` object. The argument `S` is an `AbstractStatistic` object. The argument `upw` determines whether the bootstrap is one-sided and upper-tailed or lower-tailed. The argument `B` is an integer indicating the number of bootstrap replications.
"""
mutable struct OneSidedBootstrapLimit{T} <: BootstrapLimit{T}
    sim::Vector{T}
    h::T
    upw::Bool

    function OneSidedBootstrapLimit(S::AbstractStatistic, upw, B::Int)
        T = typeof(get_value(S))
        new{T}([deepcopy(get_value(S)) for _ in 1:B], deepcopy(get_value(S)), upw)
    end
end
export OneSidedBootstrapLimit


"""
    TwoSidedBootstrapLimit{T} <: BootstrapLimit{T}

A two-sided bootstrap limit with constant false-alarm rate.

# Fields
- `sim::Vector{T}`: The vector of simulated statistics.
- `h::Vector{T}`: The current value of the control limits.

# Constructors
- `TwoSidedBootstrapLimit(S::AbstractStatistic, B::Int)`: Create a new `TwoSidedBootstrapLimit` object. The argument `S` is an `AbstractStatistic` object. The argument `B` is an integer indicating the number of bootstrap replications.
"""
mutable struct TwoSidedBootstrapLimit{T} <: BootstrapLimit{T}
    sim::Vector{T}
    h::Vector{T}

    function TwoSidedBootstrapLimit(S::AbstractStatistic, B::Int)
        T = typeof(get_value(S))
        new{T}([deepcopy(get_value(S)) for _ in 1:B], [deepcopy(get_value(S)) for _ in 1:2])
    end
end
export TwoSidedBootstrapLimit


"""
    update_value!(L::BootstrapLimit, NM::ARL)
    update_value!(L::BootstrapLimit, NM::QRL)

Update the value of the `BootstrapLimit` object using the nominal value obtained from the `ARL` or `QRL` object.

# Arguments
- `L::BootstrapLimit`: The `BootstrapLimit` object to be updated.
- `NM`: The `ARL` or `QRL` object used to obtain the value for the update.
"""
update_value!(L::BootstrapLimit, NM::ARL) = update_value!(L, 1.0/get_value(NM))

function update_value!(L::BootstrapLimit, NM::QRL)
    @assert NM.qtl == 0.5 "Dynamic limits for `QRL` objects are only implemented for the median."
    update_value!(L, 1.0 - 2.0^(-1.0/get_value(NM)))
end
export update_value!


function update_value!(L::OneSidedBootstrapLimit, alpha::Float64)
    if L.upw
        return set_h!(L, quantile(L.sim, 1.0 - alpha))
    else
        return set_h!(L, quantile(L.sim, alpha))
    end
end


function update_value!(L::TwoSidedBootstrapLimit, alpha::Float64)
    return set_h!(L, quantile(L.sim, [alpha/2.0, 1.0 - alpha/2.0]))
end

function is_IC(L::TwoSidedBootstrapLimit, stat::AbstractStatistic)
    val = get_value(stat)
    lim = get_value(L)
    @assert all(typeof(val) .== typeof.(lim))
    return compare_values(lim, val, L)
end

function compare_values(lim_val, stat_val, L::OneSidedBootstrapLimit)
    for i in 1:length(lim_val)
        if L.upw[i]
            stat_val[i] <= lim_val[i] || return false
        else
            stat_val[i] >= lim_val[i] || return false
        end
    end
    return true
end
 
function compare_values(lim_val, stat_val, L::TwoSidedBootstrapLimit)
    if (stat_val < lim_val[1]) || (stat_val > lim_val[2])
        return false
    end
    return true
end


function resample_sims!(L::BootstrapLimit)
    idx = [compare_values(get_value(L), L.sim[b], L) for b in 1:length(L.sim)]
    L.sim .= StatsBase.sample(view(L.sim, idx), length(L.sim))
end
export resample_sims!