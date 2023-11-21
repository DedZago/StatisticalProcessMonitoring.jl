
module TestAlgorithms
using SPM
using Test
using Distributions

@testset "settings" begin
    @testset "constructor" begin
        OptSettings()
    end
end

@testset "control limit opt" begin
    x = randn(500)
    NM = ARL(200)
    PH1 = Phase2(Bootstrap(), x)
    STAT = EWMA(λ = 0.2)
    f(t) = sqrt(STAT.λ/(2.0 - STAT.λ) * (1.0 - (1.0 - STAT.λ)^(2.0*t)))
    LIM = TwoSidedCurvedLimit(2.0, f)
    CH = ControlChart(STAT, LIM, NM, PH1)
    @testset "call" begin
        saCL(CH, Nmin=1, maxiter=1, Nfixed=1)
        bisectionCL(CH, 1.0, maxiter=1, nsims = 10)
        combinedCL(CH, maxiter=1, maxiter_sa=1, nsims=1)
        bootstrapCL(CH, nsims=1, maxiter=1)
        saCL(CH, Nmin=1, maxiter=1, Nfixed=1, parallel=true)
        bootstrapCL(CH, nsims=1, maxiter=1, parallel=true)
        optimize_limit(CH, :SA, Nmin=1, maxiter=1, Nfixed=1)
        optimize_limit(CH, :Bisection, hmax=1.0, maxiter=1, nsims = 10)
        optimize_limit(CH, :Combined, maxiter=1, maxiter_sa=1, nsims=1)
        optimize_limit(CH, :Bootstrap, nsims=1, maxiter=1)
    end

    x = randn(500)
    p = 0.5
    NM = QRL(200, p)
    PH1 = Phase2(Bootstrap(), x)
    STAT = CUSUM(k = 0.5)
    LIM = OneSidedFixedLimit(2.0, true)
    CH = ControlChart(STAT, LIM, NM, PH1)
    @testset "call" begin
        saCL(CH, Nmin = 1, maxiter = 1, Nfixed = 1)
        bisectionCL(CH, 1.0, maxiter=1, nsims = 10)
        combinedCL(CH, maxiter=1, maxiter_sa=1, nsims=1)
    end

    x = randn(500)
    p = 0.5
    NM = QRL(200, p)
    PH1 = Phase2(Bootstrap(), x)
    STAT = CUSUM(k = 0.5)
    LIM = OneSidedFixedLimit(2.0, true)
    CH = ControlChart(STAT, LIM, NM, PH1)
    @testset "call" begin
        saCL(CH, Nmin=1, maxiter=1, Nfixed=1)
        bisectionCL(CH, 1.0, maxiter=1, nsims = 10)
        combinedCL(CH, maxiter=1, maxiter_sa=1, nsims=1)
    end

    x = randn(500)
    NM = ARL(200)
    STAT1 = EWMA(λ = 0.1)
    STAT2 = CUSUM(k = 1.0)
    g(t) = sqrt(STAT1.λ/(2.0 - STAT1.λ) * (1.0 - (1.0 - STAT1.λ)^(2.0*t)))
    LIM1 = TwoSidedCurvedLimit(2.0, g)
    LIM2 = OneSidedFixedLimit(2.3, true)
    PH1 = Phase2(Bootstrap(), x)
    CH = ControlChart((STAT1, STAT2), (LIM1, LIM2), NM, PH1)
    @testset "call" begin
        saCL(CH, Nmin=1, maxiter=1, Nfixed=1)
    end
end


@testset "multiple charts" begin
    x = randn(500)
    NM = ARL(200)
    PH1 = Phase2Distribution(Normal(0,1))
    STAT1 = EWMA(λ = 0.2)
    STAT2 = Shewhart()
    LIM1 = TwoSidedFixedLimit(2.0)
    LIM2 = TwoSidedFixedLimit(2.0)
    CH = ControlChart((STAT1, STAT2), (LIM1, LIM2), NM, PH1)
    @testset "call" begin
        saCL(CH, Nmin=1, maxiter=1, Nfixed=1)
        bootstrapCL(CH, nsims=1, maxiter=1)
        optimize_limit(CH, :SA, Nmin=1, maxiter=1, Nfixed=1)
        optimize_limit(CH, :Bootstrap, nsims=1, maxiter=1)
    end

    x = randn(1000)
    NM = QRL(200, 0.5)
    PH1 = Phase2(BlockBootstrap(10, x), x)
    CH = ControlChart((STAT1, STAT2), (LIM1, LIM2), NM, PH1)
    @testset "call" begin
        saCL(CH, Nmin=1, maxiter=1, Nfixed=1)
        bootstrapCL(CH, nsims=1, maxiter=1)
        optimize_limit(CH, :SA, Nmin=1, maxiter=1, Nfixed=1)
        optimize_limit(CH, :Bootstrap, nsims=1, maxiter=1)
    end
end
end#module