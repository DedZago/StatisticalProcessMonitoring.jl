abstract type AbstractLimit{T} end
get_h(L::AbstractLimit) = L.h

# function get_h(L::Vector{T}) where T <: AbstractLimit
#     ret = Vector{Float64}(undef, length(L))
#     for i in 1:length(L)
#         ret[i] = get_h(L[i])
#     end
#     return ret
# end

function get_h(L::Tuple)
    ret = Vector{Float64}(undef, length(L))
    for i in 1:length(L)
        ret[i] = get_h(L[i])
    end
    return ret
end
export get_h

get_value(L::AbstractLimit) = get_h(L)
get_value(L::Vector{AbstractLimit{T}}) where T = get_h(L)
get_value(L::NTuple{N, AbstractLimit}) where N = get_h(L)
export get_value
set_h!(L::AbstractLimit{T}, h::T) where T = L.h = h
set_h!(L::AbstractLimit{T}, h::Vector{T}) where T = L.h = h
function set_h!(L::Tuple, h)
    for i in eachindex(L)
        set_h!(L[i], h)
    end
end
export set_h!

function is_IC(L::AbstractLimit, stat::AbstractStatistic)
    val = get_value(stat)
    lim = get_value(L)
    return compare_values(lim, val, L)
end

is_IC(L::Tuple, stat::Tuple) = all(is_IC_vec(L, stat))
export is_IC

is_OC(L::AbstractLimit, stat::AbstractStatistic) = !is_IC(L, stat)
export is_OC

# function is_IC_vec(L::Vector{LIM}, stat::Vector{STAT}) where {LIM <: AbstractLimit, STAT <: AbstractStatistic} 
#     output = Vector{Bool}(undef, length(L))
#     for i in eachindex(output)
#         output[i] = is_IC(L[i], stat[i])
#     end
#     return output
# end

function is_IC_vec(L::Tuple, stat::Tuple)
    output = Vector{Bool}(undef, length(L))
    for i in eachindex(output)
        output[i] = is_IC(L[i], stat[i])
    end
    return output
end
export is_IC_vec

export get_curved_value

include("fixed.jl")
include("dynamic.jl")