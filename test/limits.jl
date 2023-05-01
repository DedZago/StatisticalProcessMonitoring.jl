module TestLimits
using SPM
using Test

@testset "Fixed" begin
    h = 1.0; upw = true
    @testset "One-sided limit constructors" begin
        L = OneSidedLimit(h)
        L = OneSidedLimit(h, upw=false)
        L = OneSidedLimit([h], [upw])
        L = OneSidedLimit([h for _ in 1:10], [upw for _ in 1:10])
        L = OneSidedLimit([h for _ in 1:10], [!upw for _ in 1:10])
        @test_throws AssertionError OneSidedLimit([h for _ in 1:5], [!upw for _ in 1:10])
    end
    @testset "Two-sided limit constructors" begin
        L = TwoSidedLimit(h)
        L = TwoSidedLimit([h])
        L = TwoSidedLimit([h for _ in 1:10])
        @test_throws AssertionError TwoSidedLimit([-h for _ in 1:5])
    end
end
end