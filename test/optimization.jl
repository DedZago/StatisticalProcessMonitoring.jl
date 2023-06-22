module TestOptimization
using Test
using LinearAlgebra

# Unit test

@testset "optimize_design! function" begin
    # Define a simulation to optimize
    Random.seed!(123)
    x = randn(100)
    NM = ARL(200)
    PH1 = Phase1Data(Bootstrap(), x)
    LIM = TwoSidedFixedLimit(1.0)
    λ0 = 0.2
    STAT = EWMA(λ = λ0)
    CH = ControlChart(STAT, LIM, NM, PH1)

    # Test happy path
    @testset "happy path" begin
        # Optimize with NLOpt
        optNL = optimize_design(CH, (x)->run_sim_oc(x, shift=0.5), settings=OptSettings(maxiter_opt = 4, minpar_opt = [0.001], maxpar_opt = [0.999]))
        @test optNL != get_design(CH)
        # Optimize with grid search
        optGR = optimize_design(CH, (x)->run_sim_oc(x, shift=0.5), settings=OptSettings(method_opt = :Grid, m_grid=3, maxiter_opt = 1, minpar_opt = [0.001], maxpar_opt = [0.999]))
        @test optGR != get_design(CH)
        # # Optimize with SPSA
        # optSPSA = optimize_design(CH, (x)->run_sim_oc(x, shift=0.5), settings=OptSettings(maxiter_opt = 2, minpar_opt = [0.001], maxpar_opt = [0.999]))
    end

end
end#module