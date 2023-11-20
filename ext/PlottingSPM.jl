module PlottingSPM # Should be same name as the file (just like a normal package)

using SPM, Plots

function SPM._get_row_limit_value(lim_row_i::AbstractArray, stat_vec)
    if eltype(lim_row_i) <: OneSidedCurvedLimit || eltype(lim_row_i) <: TwoSidedCurvedLimit 
        t = 1:length(lim_row_i)
        return hcat(get_value.(lim_row_i, t, stat_vec)...)
    end
    return hcat(get_value.(lim_row_i)...)
end

function SPM.plot_series(proc::ProcessControl; kw...)
    lims = hcat(vec(proc.lim)...)
    lim_row = eachrow(lims)
    plt = plot(get_value.(proc.stat); kw...)
    for i in 1:length(lim_row)
        limits = _get_row_limit_value(lim_row[i], proc.stat)
        for h in axes(limits, 1)
            plot!(plt, limits[h,:], label = "", linestyle=:dash, colour="black")
        end
    end
    plt
end


function SPM.plot_series(proc::ProcessControl{X,S,I,Vector{T}}; kw...) where {X,S,I,T <: NTuple}
    y = ones(3) 
    keywords = Dict(kw)
    smaller_titlefontsize = 14
    if haskey(keywords, :titlefontsize)
        smaller_titlefontsize = keywords[:titlefontsize]
    end
    smaller_titlefontsize = ceil(3/5 * smaller_titlefontsize)
    stats = hcat(collect.(proc.stat)...)
    vals = hcat(get_value.(proc.stat)...)
    lims = hcat(collect.(proc.lim)...)
    plt_i = []
    stat_row = eachrow(stats)
    val_row = eachrow(vals)
    lim_row = eachrow(lims)
    subtitles = "Chart " .* string.(1:length(val_row))
    if haskey(keywords, :subtitles)
        subtitles = keywords[:subtitles]
        delete!(keywords, :subtitles)
    end
    if length(subtitles) != length(val_row)
        error("Please provide a subtitle for each plot. Provided subtitles: $(length(subtitles)), required: $(length(val_row)).")
    end
    for i in 1:length(val_row)
        push!(plt_i, plot(val_row[i]; title=subtitles[i], label="", titlefontsize=smaller_titlefontsize))

        limits = _get_row_limit_value(lim_row[i], stat_row[i])
        for h in axes(limits, 1)
            plot!(plt_i[i], limits[h,:], label = "", linestyle=:dash, colour="black")
        end
    end
    if haskey(keywords, :title)
        global_title = keywords[:title]
        delete!(keywords, :title)
        plt = plot(plt_i...; keywords...)
        # create a transparent scatter plot with an 'annotation' that will become title
        title = Plots.scatter(y, marker=0,markeralpha=0, annotations=(2, y[2], Plots.text(global_title)),axis=false, grid=false, leg=false,size=(200,100))
        return plot(title, plt, layout = grid(2,1,heights=[0.01,0.9]))
    else
        plt = plot(plt_i...; keywords...)
        return plt
    end
end


end # module
