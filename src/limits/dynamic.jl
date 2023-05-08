using StatsBase

abstract type DynamicLimit <: AbstractLimit end
abstract type BootstrapLimit <: DynamicLimit end

Base.show(io::IO, L::BootstrapLimit) = print(io, "$(typeof(L))\n  value: $(get_value(L))\n  Bootstrap samples: $(length(L.sim))\n")

mutable struct OneSidedBootstrapLimit{T} <: BootstrapLimit
    sim::Vector{T}
    h::T
    upw::Bool

    function OneSidedBootstrapLimit(S::AbstractStatistic, upw, B::Int)
        T = typeof(get_value(S))
        new{T}([deepcopy(get_value(S)) for _ in 1:B], deepcopy(get_value(S)), upw)
    end
end
export OneSidedBootstrapLimit

mutable struct TwoSidedBootstrapLimit{T} <: BootstrapLimit
    sim::Vector{T}
    h::Vector{T}

    function TwoSidedBootstrapLimit(S::AbstractStatistic, B::Int)
        T = typeof(get_value(S))
        new{T}([deepcopy(get_value(S)) for _ in 1:B], [deepcopy(get_value(S)) for _ in 1:2])
    end
end
export TwoSidedBootstrapLimit


"""
    get_value(L::BootstrapLimit, NM::ARL)
    get_value(L::BootstrapLimit, NM::QRL)
"""
update_value!(L::BootstrapLimit, NM::ARL) = update_value!(L, 1.0/get_value(NM)) #FIXME: test

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
    L.sim[:] = StatsBase.sample(view(L.sim, idx), length(L.sim))
end
export resample_sims!