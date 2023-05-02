module TestCharts
using SPM
using Test

@testset "Charts" begin
    x = randn(100)
    NM = ARL(200)
    PH1 = Phase1Data(x)
    LIM = TwoSidedLimit(1.0)
    STAT = EWMA(λ = 0.2)
    @testset "constructor" begin
        CH = ControlChart(STAT, LIM, NM, PH1)
        @test is_IC(CH)
        @test !is_OC(CH)
        @test get_param(CH) == (λ = 0.2,)
        @test get_phase1(CH) == PH1
        @test get_limit_value(CH) == [1.0]
        @test get_value(CH) == 0.0
        @test get_statistic(CH) == STAT
        @test get_nominal(CH) == NM
        @test get_nominal_value(CH) == 200.0
        update_chart!(CH, 1.0)
        @test get_value(CH) == 0.2
        update_chart!(CH, 100.0)
        @test is_OC(CH)
        @test !is_IC(CH)
        set_nominal!(CH, NM)
        set_phase1!(CH, PH1)
        set_statistic!(CH, STAT)
        set_limit!(CH, LIM)
        set_param!(CH, 0.3)
        @test get_maxrl(CH) == Inf
    end
end
end