using SPM
using Test

ti = time()

@testset "SPM.jl" begin
    include("limits.jl")
    include("statistics.jl")
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
