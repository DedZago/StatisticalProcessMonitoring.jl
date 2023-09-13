module TestLimits
using SPM
using Test
using StatsBase
using Distributions

@testset "Fixed" begin
    h = 1.0; upw = true
    @testset "One-sided fixed limit constructors" begin
        L = OneSidedFixedLimit(h, true)
        @test get_h(L) == h
        @test get_h(L) == get_value(L)
        hup = 5.0
        set_h!(L, hup)
        @test get_h(L) == hup
        L = OneSidedFixedLimit(h=h, upw=false)
        @test get_h(L) == h
        @test get_h(L) != get_value(L)
        @test get_value(L) == -h
        hup = 5.0
        set_h!(L, hup)
        @test get_h(L) == hup
        @test get_h(L) != get_value(L)
        @test get_value(L) == -hup
        @test_throws AssertionError OneSidedFixedLimit(-0.1, true)
        @test_throws AssertionError OneSidedFixedLimit(-0.1, false)
    end

    @testset "Two-sided fixed limit constructors" begin
        L = TwoSidedFixedLimit(h)
        @test get_h(L) == h
        @test [-get_h(L), get_h(L)] == get_value(L)
        hup = 5.0
        set_h!(L, hup)
        @test get_h(L) == hup
        @test [-get_h(L), get_h(L)] == get_value(L)
        @test get_value(L) == [-hup, hup]
    end

    @testset "isOC statistic" begin
        STAT = EWMA(λ=0.2, value=1.0)
        L1 = OneSidedFixedLimit(1.5, true)
        @test is_IC(L1, STAT)
        @test !is_OC(L1, STAT)
        L2 = OneSidedFixedLimit(0.5, true)
        @test is_OC(L2, STAT)
        @test !is_IC(L2, STAT)
        L3 = TwoSidedFixedLimit(0.5)
        @test is_OC(L3, STAT)
        @test !is_IC(L3, STAT)
        L4 = TwoSidedFixedLimit(3.0)
        @test is_IC(L4, STAT)
        @test !is_OC(L4, STAT)
        STAT = EWMA(λ=0.2, value=-1.0)
        L5 = TwoSidedFixedLimit(0.8)
        @test !is_IC(L5, STAT)
        @test is_OC(L5, STAT)
    end

    @testset "One-sided curved" begin
        f(t, STAT) = sqrt(STAT.λ/(2.0 - STAT.λ) * (1.0 - (1.0 - STAT.λ)^(2.0*t)))
        λ = 0.2
        STAT = EWMA(λ = λ, value = 0.0)
        h = 0.5
        L = OneSidedCurvedLimit(h, true, f)
        @test get_value(L) == h
        @test get_value(L, 0, STAT) == h * f(0, STAT) 
        @test get_value(L, 1, STAT) == h * f(1, STAT) 
        @test get_value(L, 10^5, STAT) == h * sqrt(λ/(2-λ)) 
        L = OneSidedCurvedLimit(h, false, f)
        @test get_value(L) == -h
        @test get_value(L, 0, STAT) == -h * f(0, STAT) 
        @test get_value(L, 1, STAT) == -h * f(1, STAT) 
        @test get_value(L, 10^5, STAT) == -h * sqrt(λ/(2-λ)) 
        @test_throws AssertionError OneSidedCurvedLimit(-0.5, true, f)
    end
    @testset "Two-sided curved" begin
        f(t, STAT) = sqrt(STAT.λ/(2.0 - STAT.λ) * (1.0 - (1.0 - STAT.λ)^(2.0*t)))
        λ = 0.2
        STAT = EWMA(λ = λ, value = 0.0)
        h = 0.5
        L = TwoSidedCurvedLimit(h, f)
        @test get_h(L) == h
        @test get_value(L, 0, STAT) == h * f(0, STAT) * [-1, 1]
        @test get_value(L, 1, STAT) == h * f(1, STAT) * [-1, 1]
        @test get_value(L, 10^5, STAT) == h * sqrt(λ/(2-λ)) * [-1, 1]
        @test_throws MethodError TwoSidedCurvedLimit(-0.5, true, f)
    end

    @testset "Dynamic limits" begin
        STAT = EWMA(λ = 1.0)
        B = 1000
        L = OneSidedBootstrapLimit(STAT, true, B)
        @test length(L.sim) == B
        @test get_value(L) == 0.0
        NM = ARL(200)
        update_value!(L, NM)
        @test get_value(L) == 0.0
        x = randn(B)
        L.sim = x
        update_value!(L, NM)
        @test get_value(L) == quantile(x, 1.0 - 1.0/get_value(NM))
        resample_sims!(L)
        @test all([L.sim[b] in x for b in 1:B])
        L = OneSidedBootstrapLimit(STAT, false, B)
        is_IC(L, STAT)

        L = OneSidedBootstrapLimit(STAT, false, B)
        NM = ARL(100)
        update_value!(L, NM)
        @test get_value(L) == 0.0
        x = randn(B)
        L.sim = x
        update_value!(L, NM)
        @test get_value(L) == quantile(x, 1.0/get_value(NM))
        @test all([L.sim[b] in x for b in 1:B])

        L = TwoSidedBootstrapLimit(STAT, B)
        NM = ARL(100)
        update_value!(L, NM)
        @test get_value(L) == [0.0, 0.0]
        x = randn(B)
        L.sim = x
        update_value!(L, NM)
        alpha = 1.0/get_value(NM)
        @test get_value(L) == quantile(x, [alpha/2, 1.0 - alpha/2])
        @test all([L.sim[b] in x for b in 1:B])
        is_IC(L, STAT)

        PH2 = Phase2Distribution(Normal(0,1))
        CH = ControlChart(STAT, L, NM, PH2)
        update_chart!(CH, rand(Normal(0,1)))
    end
end
end