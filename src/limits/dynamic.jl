abstract type DynamicLimit <: AbstractLimit end
abstract type OneSidedDynamicLimit <: DynamicLimit end
abstract type TwoSidedDynamicLimit <: DynamicLimit end

"""
    get_value(L::DynamicLimit, NM::ARL)
    get_value(L::DynamicLimit, NM::QRL)
"""

get_value(L::DynamicLimit, NM::ARL) = get_value(L, 1.0/get_value(NM)) #FIXME: test

function get_value(L::DynamicLimit, NM::QRL)
    @assert NM.qtl == 0.5 "Dynamic limits for `QRL` objects are only implemented for the median."
    get_value(L, 1.0 - 2.0^(-1.0/get_value(NM)))
end#FIXME: test


function get_value(L::OneSidedDynamicLimit, alpha::Float64)
    if L.upw
        return quantile(L.sim, 1.0 - alpha)       
    else
        return quantile(L.sim, alpha)       
    end
end#FIXME: test


function get_value(L::TwoSidedDynamicLimit, alpha::Float64)
    return quantile(L.sim, [alpha/2.0, 1.0 - alpha/2.0])       
end#FIXME: test


@with_kw mutable struct OneSidedBootstrapLimit{T} <: OneSidedDynamicLimit
    sim::Vector{T}
    upw::Bool
end#FIXME: test

#TODO: Define update_limit! method to perform bootstrap

OneSidedBootstrapLimit(S::AbstractStatistic, upw, B::Int) = OneSidedBootstrapLimit([deepcopy(get_value(S)) for _ in 1:B], upw)
export OneSidedBootstrapLimit


@with_kw mutable struct TwoSidedBootstrapLimit{T} <: TwoSidedDynamicLimit
    sim::Vector{T}
end#FIXME: test

TwoSidedBootstrapLimit(S::AbstractStatistic, B::Int) = TwoSidedBootstrapLimit([deepcopy(get_value(S)) for _ in 1:B])
export TwoSidedBootstrapLimit

#TODO: Define update_limit! method to perform bootstrap