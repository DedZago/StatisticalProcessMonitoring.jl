module TestLimits
using SPM
using Test

@testset "Fixed" begin
    h = 1.0; upw = true
    @testset "One-sided limit constructors" begin
        L = OneSidedFixedLimit(h)
        @test get_value(L) == [h]
        hup = 5.0
        set_value!(L, hup)
        @test get_value(L) == [hup]
        L = OneSidedFixedLimit(h, upw=false)
        L = OneSidedFixedLimit([h], [upw])
        L = OneSidedFixedLimit([h for _ in 1:10], [upw for _ in 1:10])
        L = OneSidedFixedLimit([h for _ in 1:10], [!upw for _ in 1:10])
        @test_throws AssertionError OneSidedFixedLimit([h for _ in 1:5], [!upw for _ in 1:10])
    end
    @testset "Two-sided limit constructors" begin
        L = TwoSidedFixedLimit(h)
        L = TwoSidedFixedLimit([h])
        L = TwoSidedFixedLimit([h for _ in 1:10])
        @test_throws AssertionError TwoSidedFixedLimit([-h for _ in 1:5])
    end
    @testset "isOC statistic" begin
        STAT = EWMA(λ=0.2, value=1.0)
        L1 = OneSidedFixedLimit([1.5], [true])
        @test is_IC(L1, STAT)
        @test !is_OC(L1, STAT)
        L2 = OneSidedFixedLimit([0.5], [true])
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
end
end