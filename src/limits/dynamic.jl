abstract type DynamicLimit <: AbstractLimit end
abstract type BootstrapLimit <: DynamicLimit end
# abstract type OneSidedBootstrapLimit <: DynamicLimit end
# abstract type TwoSidedDynamicLimit <: DynamicLimit end

Base.show(io::IO, L::BootstrapLimit) = print(io, "$(typeof(L))\n  value: $(get_value(L))\n  sims: $(length(L.sim))\n")

mutable struct OneSidedBootstrapLimit{T} <: BootstrapLimit
    sim::Vector{T}
    value::T
    upw::Bool

    function OneSidedBootstrapLimit(S::AbstractStatistic, upw, B::Int)
        T = typeof(get_value(S))
        new{T}([deepcopy(get_value(S)) for _ in 1:B], deepcopy(get_value(S)), upw)
    end
end#FIXME: test
export OneSidedBootstrapLimit

mutable struct TwoSidedBootstrapLimit{T} <: BootstrapLimit
    sim::Vector{T}
    value::Vector{T}

    function TwoSidedBootstrapLimit(S::AbstractStatistic, upw, B::Int)
        T = typeof(get_value(S))
        new{T}([deepcopy(get_value(S)) for _ in 1:B], [deepcopy(get_value(S)) for _ in 1:2])
    end
end#FIXME: test
export TwoSidedBootstrapLimit


"""
    get_value(L::BootstrapLimit, NM::ARL)
    get_value(L::BootstrapLimit, NM::QRL)
"""

update_value!(L::BootstrapLimit, NM::ARL) = update_value!(L, 1.0/get_value(NM)) #FIXME: test

function update_value!(L::BootstrapLimit, NM::QRL)
    @assert NM.qtl == 0.5 "Dynamic limits for `QRL` objects are only implemented for the median."
    update_value!(L, 1.0 - 2.0^(-1.0/get_value(NM)))
end#FIXME: test


function update_value!(L::OneSidedBootstrapLimit, alpha::Float64)
    if L.upw
        return set_value!(L, quantile(L.sim, 1.0 - alpha))
    else
        return set_value!(L, quantile(L.sim, alpha))
    end
end#FIXME: test


function update_value!(L::TwoSidedBootstrapLimit, alpha::Float64)
    return set_value!(L, quantile(L.sim, [alpha/2.0, 1.0 - alpha/2.0]))
end#FIXME: test

function resample_sims!(L::OneSidedBootstrapLimit)
    error("Not implemented yet.")
end

function resample_sims!(L::TwoSidedBootstrapLimit)
    error("Not implemented yet.")
end


#TODO: Define update_limit! method to perform bootstrap





#TODO: Define update_limit! method to perform bootstrap