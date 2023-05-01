@with_kw mutable struct OneSidedLimit{T} <: AbstractLimit
    value::Vector{T}
    upw::Vector{Bool} = ones(length(value))
    @assert length(value) == length(upw)
end


function is_IC(L::OneSidedLimit, stat::AbstractStatistic)
    #TODO: run tests to check if it is correct
    val = get_value(stat)
    lim = get_limit_value(L)
    @assert length(val) == length(lim)
    for i in eachindex(val)
        if upw[i]
            val[i] < lim[i] || return false
        else
            val[i] > lim[i] || return false
        end
    end
    return true
end