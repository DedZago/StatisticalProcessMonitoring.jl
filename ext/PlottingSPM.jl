module PlottingSPM # Should be same name as the file (just like a normal package)

using SPM, Plots

#TODO: write plot_series for single and multiple control charts separately
function SPM.plot_series(proc::ProcessControl; kw...)
    lims = hcat(vec(proc.lim)...)
    lim_row = eachrow(lims)
    plt = plot(proc.stat)
    for i in 1:length(lim_row)
        limits = hcat(get_value.(lim_row[i])...)
        for h in axes(limits, 1)
            plot!(plt, limits[h,:], label = "", linestyle=:dash, colour="black")
        end
    end
    plt
end

function SPM.plot_series(proc::ProcessControl{X,S,I,Vector{T}}; kw...) where {X,S,I,T <: NTuple}
    vals = hcat(proc.stat...)
    lims = hcat(collect.(proc.lim)...)
    plt_i = []
    val_row = eachrow(vals)
    lim_row = eachrow(lims)
    for i in 1:length(val_row)
        push!(plt_i, plot(val_row[i]; title="Chart " * string(i), label=""))
        limits = hcat(get_value.(lim_row[i])...)
        for h in axes(limits, 1)
            plot!(plt_i[i], limits[h,:], label = "", linestyle=:dash, colour="black")
        end
    end
    plt = plot(plt_i...; kw...)
    plt
end


end # module
