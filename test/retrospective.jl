module TestRetrospective
using StatisticalProcessMonitoring
using Test
using Random

@testset "Retrospective" begin
    Random.seed!(123)
    x = randn(100)
    NM = ARL(200)
    PH1 = Phase2(Bootstrap(), x)
    LIM = TwoSidedFixedLimit(1.0)
    STAT = EWMA(λ = 0.2)
    CH = ControlChart(STAT, LIM, NM, PH1)

    y = randn(20)
    proc = apply_chart(CH, y)

    x = randn(100, 3)
    NM = ARL(200)
    PH1 = Phase2(Bootstrap(), x)
    LIM = OneSidedFixedLimit(1.0, true)
    STAT = LLCUSUM(0.2, x)
    CH = ControlChart(STAT, LIM, NM, PH1)

    y = randn(20, 3)
    proc = apply_chart(CH, y)
end
end#module
