using StatisticalProcessMonitoring
using Test

ti = time()

@testset "StatisticalProcessMonitoring.jl" begin
    include("nominal.jl")
    include("phase2.jl")
    include("limits.jl")
    include("statistics.jl")
    include("charts.jl")
    include("settings.jl")
    include("algorithms.jl")
    include("optimization.jl")
    include("retrospective.jl")
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
