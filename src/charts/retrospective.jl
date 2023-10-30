using DataFrames
struct ProcessControl{D, S, A, L}
    x::D                # Observed data
    stat::S             # Value of the statistic
    is_OC::A
    lim::L
end
export ProcessControl

Base.show(io::IO, P::ProcessControl) = print(io, "ProcessControl of length $(length(P.stat))\n x: $(typeof(P.x))\n stat: $(typeof(P.stat))\n is_OC: $(typeof(P.is_OC))\n lim: $(typeof(P.lim))")


"""
    apply_chart(CH::AbstractChart, x::AbstractVecOrMat)
    apply_chart!(CH::AbstractChart, x::AbstractVector)
    apply_chart!(CH::AbstractChart, x::AbstractMatrix)

Apply a control chart to a data vector or data matrix `x`.
"""
function apply_chart!(CH::AbstractChart, x::AbstractVector)
    y = Vector{typeof(get_value(CH))}(undef, length(x))
    lim = Vector{typeof(get_limit(CH))}(undef, length(x))
    alarm = [false for _ in eachindex(x)]
    for i in eachindex(x)
        update_chart!(CH, x[i])
        y[i] = deepcopy(get_value(CH))
        lim[i] = deepcopy(get_limit(CH))
        if is_OC(CH)
            alarm[i] = true
        end
    end
    return ProcessControl(deepcopy(x), y, alarm, lim)
end


function apply_chart!(CH::AbstractChart, x::AbstractMatrix)
    n, _ = size(x)
    y = Vector{typeof(get_value(CH))}(undef, n)
    lim = Vector{typeof(get_limit(CH))}(undef, n)
    alarm = [false for _ in 1:n]
    for i in 1:n
        update_chart!(CH, view(x, i, :))
        y[i] = deepcopy(get_value(CH))
        lim[i] = deepcopy(get_limit(CH))
        if is_OC(CH)
            alarm[i] = true
        end
    end
    return ProcessControl(deepcopy(x), y, alarm, lim)
end

function apply_chart!(CH::AbstractChart, x::DataFrame)
    n, _ = size(x)
    y = Vector{typeof(get_value(CH))}(undef, n)
    lim = Vector{typeof(get_limit(CH))}(undef, n)
    alarm = [false for _ in 1:n]
    for i in 1:n
        update_chart!(CH, DataFrame(view(x, i, :)))
        y[i] = deepcopy(get_value(CH))
        lim[i] = deepcopy(get_limit(CH))
        if is_OC(CH)
            alarm[i] = true
        end
    end
    return ProcessControl(deepcopy(x), y, alarm, lim)
end
export apply_chart!




function apply_chart(CH::AbstractChart, x::Union{AbstractVecOrMat, DataFrame})
    CH_ = deepcopy(CH)
    return apply_chart!(CH_, x)
end
export apply_chart

function plot_series(PCTL; kw_ind::Dict = Dict(), kw_main::Dict = Dict())
    error("Package Plots is not loaded.")
end
export plot_series