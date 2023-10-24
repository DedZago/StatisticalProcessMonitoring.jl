struct ProcessControl{D<:AbstractVecOrMat, S, A, L}
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
export apply_chart!


function apply_chart(CH::AbstractChart, x::AbstractVecOrMat)
    CH_ = deepcopy(CH)
    return apply_chart!(CH_, x)
end
export apply_chart

#FIXME: construct interface for plotting and diagnosing.

# function plot_series(PCTL::ProcessControl; kw_ind::Dict = Dict(), kw_main::Dict = Dict())
#     vals = hcat(PCTL.stat...)
#     lims = hcat(PCTL.lim...)
#     plt_i = []
#     val_row = eachrow(vals)
#     lim_row = eachrow(lims)
#     for i in 1:length(val_row)
#         push!(plt_i, plot(val_row[i]; title="Chart " * string(i), label=""))
#         limits = hcat(get_value.(lim_row[i])...)
#         for h in axes(limits, 1)
#             plot!(plt_i[i], limits[h,:], label = "", linestyle=:dash)
#         end
#     end
#     plt = plot(plt_i...; kw_main...)
#     plt
# end
# export plot_series