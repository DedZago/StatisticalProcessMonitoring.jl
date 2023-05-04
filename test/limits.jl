module TestLimits
using SPM
using Test

@testset "Fixed" begin
    h = 1.0; upw = true
    @testset "One-sided fixed limit constructors" begin
        L = OneSidedFixedLimit(h, true)
        @test get_value(L) == h
        hup = 5.0
        set_value!(L, hup)
        @test get_value(L) == hup
        L = OneSidedFixedLimit(value=h, upw=false)
        OneSidedFixedLimit(-0.1, true)
    end
    @testset "Two-sided fixed limit constructors" begin
        L = TwoSidedFixedLimit(h)
        @test_throws AssertionError TwoSidedFixedLimit(-0.1)
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
        L = OneSidedCurvedLimit(h, true, f, STAT)
        @test get_value(L) == h
        @test get_curved_value(L, 0, STAT) == h * f(0, STAT) 
        @test get_curved_value(L, 1, STAT) == h * f(1, STAT) 
        @test get_curved_value(L, 10^5, STAT) == h * sqrt(λ/(2-λ)) 
        L = OneSidedCurvedLimit(-h, false, f, STAT)
        @test get_value(L) == -h
        @test get_curved_value(L, 0, STAT) == -h * f(0, STAT) 
        @test get_curved_value(L, 1, STAT) == -h * f(1, STAT) 
        @test get_curved_value(L, 10^5, STAT) == -h * sqrt(λ/(2-λ)) 
    end
    @testset "Two-sided curved" begin
        f(t, STAT) = sqrt(STAT.λ/(2.0 - STAT.λ) * (1.0 - (1.0 - STAT.λ)^(2.0*t)))
        λ = 0.2
        STAT = EWMA(λ = λ, value = 0.0)
        h = 0.5
        L = TwoSidedCurvedLimit(h, f, STAT)
        @test get_value(L) == h
        @test get_curved_value(L, 0, STAT) == h * f(0, STAT) 
        @test get_curved_value(L, 1, STAT) == h * f(1, STAT) 
        @test get_curved_value(L, 10^5, STAT) == h * sqrt(λ/(2-λ)) 
    end
end
end