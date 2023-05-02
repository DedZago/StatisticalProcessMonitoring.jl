module TestPhase1
using SPM
using Test

@testset "Phase1Data" begin
    x = [1.0, 2.0, 3.0]
    xmat = [1.0 2.0 3.0; 2.0 3.0 4.0]
    @testset "newdata" begin
        PH1 = Phase1Data(x)
        @test isa(PH1, Phase1Data{Vector{Float64}})
        y = new_data(PH1)
        @test y in x
        PH1 = Phase1Data(xmat)
        y = new_data(PH1)
        @test y in [z for z in eachrow(xmat)]
        @test isa(PH1, Phase1Data{Matrix{Float64}})
    end
end
end