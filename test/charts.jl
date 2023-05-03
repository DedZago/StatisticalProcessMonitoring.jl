module TestCharts
using SPM
using Test

@testset "Charts" begin
    x = randn(100)
    NM = ARL(200)
    PH1 = Phase1Data(x)
    LIM = TwoSidedFixedLimit(1.0)
    STAT = EWMA(λ = 0.2)
    @testset "constructor" begin
        CH = ControlChart(STAT, LIM, NM, PH1)
        CH = ControlChart(STAT, LIM, NM, PH1, 0)
        @test is_IC(CH)
        @test !is_OC(CH)
        @test get_t(CH) == 0
        @test get_param(CH) == (λ = 0.2,)
        @test get_phase1(CH) == PH1
        @test get_limit_value(CH) == 1.0
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
        set_limit!(CH, 0.1)
        @test get_limit_value(CH) == 0.1
        set_param!(CH, 0.3)
        @test get_param(CH)[1] == 0.3
        @test get_maxrl(CH) == Inf
        @test new_data(CH) in x
    end

    x = randn(100)
    NM = ARL(200)
    PH1 = Phase1Data(x)
    λ = 0.2
    STAT = EWMA(λ = λ)
    f(t, STAT) = sqrt(STAT.λ/(2.0 - STAT.λ) * (1.0 - (1.0 - STAT.λ)^(2.0*t)))
    LIM = OneSidedCurvedLimit(1.0, true, f, STAT)

    @testset "Curved chart" begin
        CH = ControlChart(STAT, LIM, NM, PH1)
        CH = ControlChart(STAT, LIM, NM, PH1, 0)
        @test get_t(CH) == 0
        @test get_param(CH) == (λ = 0.2,)
        @test get_phase1(CH) == PH1
        @test get_limit_value(CH) == 0.0
        @test get_value(CH) == 0.0
        @test get_statistic(CH) == STAT
        @test get_nominal(CH) == NM
        @test get_nominal_value(CH) == 200.0
        xnew = 0.3
        update_chart!(CH, xnew)
        @test is_IC(CH)
        @test !is_OC(CH)
        @test get_value(CH) == λ * xnew
        update_chart!(CH, 100.0)
        @test is_OC(CH)
        @test !is_IC(CH)
        set_nominal!(CH, NM)
        set_phase1!(CH, PH1)
        set_statistic!(CH, STAT)
        set_limit!(CH, LIM)
        set_limit!(CH, 0.1)
        @test get_limit_value(CH) != 0.0
        set_param!(CH, 0.3)
        @test get_param(CH)[1] == 0.3
        @test get_maxrl(CH) == Inf
        @test new_data(CH) in x
        
    end
end
end