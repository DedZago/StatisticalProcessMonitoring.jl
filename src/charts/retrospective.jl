#FIXME: construct retrospective interface for plotting and diagnosing.

"""
    apply_chart(CH::AbstractChart, x::AbstractVecOrMat)
    apply_chart!(CH::AbstractChart, x::AbstractVector)
    apply_chart!(CH::AbstractChart, x::AbstractMatrix)

Apply a control chart to a data vector or data matrix `x`.
"""
function apply_chart!(CH::AbstractChart, x::AbstractVector)
    y = [deepcopy(get_value(CH)) for _ in eachindex(x)]
    alarm = Vector{Bool}(undef, length(y))
    for i in eachindex(x)
        update_chart!(CH, x[i])
        y[i] = get_value(CH)
        alarm[i] = is_OC(CH)
    end
    return y, alarm
end


function apply_chart!(CH::AbstractChart, x::AbstractMatrix)
    n, _ = size(x)
    y = [deepcopy(get_value(CH)) for _ in 1:n]
    alarm = Vector{Bool}(undef, length(y))
    for i in 1:n
        update_chart!(CH, view(x, i, :))
        y[i] = get_value(CH)
        alarm[i] = is_OC(CH)
    end
    return y, alarm
end
export apply_chart!


function apply_chart(CH::AbstractChart, x::AbstractVecOrMat)
    CH_ = deepcopy(CH)
    return apply_chart!(CH_, x)
end
export apply_chart
