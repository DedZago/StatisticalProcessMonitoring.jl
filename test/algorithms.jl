
module TestAlgorithms
using SPM
using Test

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
    f(t, STAT) = sqrt(STAT.λ/(2.0 - STAT.λ) * (1.0 - (1.0 - STAT.λ)^(2.0*t)))
    LIM = TwoSidedCurvedLimit(2.0, f)
    CH = ControlChart(STAT, LIM, NM, PH1)
    @testset "call" begin
        saCL(CH, settings=OptSettings(Nmin_sa=1, maxiter_sa=1, Nfixed_sa=1))
        bisectionCL(CH, settings=OptSettings(hmax_bi=1.0, maxiter_bi=1, nsims_bi = 10))
        combinedCL(CH, settings=OptSettings(maxiter_bi=1, maxiter_sa=1, nsims_bi=1))
    end

    x = randn(500)
    p = 0.5
    NM = QRL(200, p)
    PH1 = Phase2(Bootstrap(), x)
    STAT = CUSUM(k = 0.5)
    LIM = OneSidedFixedLimit(2.0, true)
    CH = ControlChart(STAT, LIM, NM, PH1)
    @testset "call" begin
        saCL(CH, settings=OptSettings(Nmin_sa=1, maxiter_sa=1, Nfixed_sa=1))
        bisectionCL(CH, settings=OptSettings(hmax_bi=1.0, maxiter_bi=1, nsims_bi = 10))
        combinedCL(CH, settings=OptSettings(maxiter_bi=1, maxiter_sa=1, nsims_bi=1))
    end

    x = randn(500)
    p = 0.5
    NM = QRL(200, p)
    PH1 = Phase2(Bootstrap(), x)
    STAT = CUSUM(k = 0.5)
    LIM = OneSidedFixedLimit(2.0, true)
    CH = ControlChart(STAT, LIM, NM, PH1)
    @testset "call" begin
        saCL(CH, settings=OptSettings(Nmin_sa=1, maxiter_sa=1, Nfixed_sa=1))
        bisectionCL(CH, settings=OptSettings(hmax_bi=1.0, maxiter_bi=1, nsims_bi = 10))
        combinedCL(CH, settings=OptSettings(maxiter_bi=1, maxiter_sa=1, nsims_bi=1))
    end

    x = randn(500)
    NM = ARL(200)
    STAT1 = EWMA(λ = 0.1)
    STAT2 = CUSUM(k = 1.0)
    LIM1 = TwoSidedCurvedLimit(2.0, f)
    LIM2 = OneSidedFixedLimit(2.3, true)
    PH1 = Phase2(Bootstrap(), x)
    CH = ControlChart([STAT1, STAT2], [LIM1, LIM2], NM, PH1)
    @testset "call" begin
        saCL(CH, settings=OptSettings(Nmin_sa=1, maxiter_sa=1, Nfixed_sa=1))
    end
end
end#module