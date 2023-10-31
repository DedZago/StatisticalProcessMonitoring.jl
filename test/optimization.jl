module TestOptimization
using Test
using LinearAlgebra
using Random
using SPM

# Unit test

@testset "optimize_design! function" begin
    # Define a simulation to optimize
    Random.seed!(123)
    x = randn(100)
    NM = ARL(200)
    PH1 = Phase2(Bootstrap(), x)
    LIM = TwoSidedFixedLimit(1.0)
    λ0 = 0.2
    STAT = EWMA(λ = λ0)
    CH = ControlChart(STAT, LIM, NM, PH1)

    # Test happy path
    settings = OptSettings(verbose=false, minpar = [0.001], maxpar = [0.99], maxiter=1, m_grid=2, nsims=3)
    rlsim_oc = x -> run_sim_oc(x, shift = 0.5)
    @testset "optimize_design" begin
        optimize_design(CH, rlsim_oc, settings, verbose=false)
        optimize_design(CH, rlsim_oc, settings, optimizer=:Grid)
    end

end
end#module