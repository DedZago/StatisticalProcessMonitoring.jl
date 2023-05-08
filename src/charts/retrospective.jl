#FIXME: construct retrospective interface for plotting and diagnosing.

struct ProcessControl{D<:AbstractVecOrMat, S, A, L}
    x::D                # Observed data
    stat::S             # Value of the statistic
    is_OC::A
    lim::L
end
export ProcessControl

Base.show(io::IO, P::ProcessControl) = println(io, "ProcessControl of length $(length(P.stat))")


"""
    apply_chart(CH::AbstractChart, x::AbstractVecOrMat)
    apply_chart!(CH::AbstractChart, x::AbstractVector)
    apply_chart!(CH::AbstractChart, x::AbstractMatrix)

Apply a control chart to a data vector or data matrix `x`.
"""
function apply_chart!(CH::AbstractChart, x::AbstractVector)
    y = [deepcopy(get_value(CH)) for _ in eachindex(x)]
    lim = [deepcopy(get_limit_value(CH)) for _ in eachindex(x)]
    alarm = [false for _ in eachindex(x)]
    for i in eachindex(x)
        update_chart!(CH, x[i])
        y[i] = get_value(CH)
        lim[i] = get_limit_value(CH)
        if is_OC(CH)
            alarm[i] = true
        end
    end
    return ProcessControl(x, y, alarm, lim)
end


function apply_chart!(CH::AbstractChart, x::AbstractMatrix)
    n, _ = size(x)
    y = [deepcopy(get_value(CH)) for _ in 1:n]
    lim = [deepcopy(get_limit_value(CH)) for _ in 1:n]
    alarm = [false for _ in eachindex(x)]
    for i in 1:n
        update_chart!(CH, view(x, i, :))
        y[i] = get_value(CH)
        lim[i] = get_limit_value(CH)
        if is_OC(CH)
            alarm[i] = true
        end
    end
    return ProcessControl(x, y, alarm, lim)
end
export apply_chart!


function apply_chart(CH::AbstractChart, x::AbstractVecOrMat)
    CH_ = deepcopy(CH)
    return apply_chart!(CH_, x)
end
export apply_chart
