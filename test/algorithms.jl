
module TestAlgorithms
using SPM
using Test

@testset "saCL" begin
    x = randn(500)
    NM = ARL(200)
    PH1 = Phase1Data(x)
    STAT = EWMA(λ = 0.2)
    f(t, STAT) = sqrt(STAT.λ/(2.0 - STAT.λ) * (1.0 - (1.0 - STAT.λ)^(2.0*t)))
    LIM = TwoSidedCurvedLimit(2.0, f, STAT)
    CH = ControlChart(STAT, LIM, NM, PH1)
    @testset "call" begin
        saCL(CH, verbose=false, Nmin=1, maxiter=1, Nfixed=1)
    end

    x = randn(500)
    p = 0.5
    NM = QRL(200, p)
    PH1 = Phase1Data(x)
    STAT = CUSUM(k = 0.5)
    LIM = OneSidedFixedLimit(2.0, true)
    CH = ControlChart(STAT, LIM, NM, PH1)
    @testset "call" begin
        saCL(CH, verbose=false, Nmin=1, maxiter=1, Nfixed=1)
    end

    x = randn(500)
    p = 0.5
    NM = QRL(200, p)
    PH1 = Phase1Data(x)
    STAT = CUSUM(k = 0.5)
    LIM = OneSidedFixedLimit(2.0, true)
    CH = ControlChart(STAT, LIM, NM, PH1)
    @testset "call" begin
        saCL(CH, verbose=false, Nmin=1, maxiter=1, Nfixed=1)
    end

    x = randn(500)
    NM = ARL(200)
    STAT1 = EWMA(λ = 0.1)
    STAT2 = CUSUM(k = 1.0)
    LIM1 = TwoSidedCurvedLimit(2.0, f, STAT)
    LIM2 = OneSidedFixedLimit(2.3, true)
    PH1 = Phase1Data(x)
    CH = ControlChart([STAT1, STAT2], [LIM1, LIM2], NM, PH1)
    @testset "call" begin
        saCL(CH, verbose=false, Nmin=1, maxiter=1, Nfixed=1)
    end
end
end#module