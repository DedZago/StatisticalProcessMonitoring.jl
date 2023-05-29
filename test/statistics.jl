module TestStatistics
using SPM
using Test

@testset "EWMA" begin
    λ = 0.1; value = 0.0
    @testset "constructors" begin
        EWMA(λ = 0.2)
        STAT = EWMA(λ, value, 0.0, 1.0)
        @test_throws AssertionError EWMA(λ=-0.1, value=0.0)
        @test_throws AssertionError EWMA(λ=1.2, value=0.0)
        @test_throws AssertionError EWMA(λ=0.1, value=Inf)
        @test_throws AssertionError EWMA(λ=0.1, value=-Inf)
        params = get_design(STAT)
        @test length(params) == 1
        @test params[:λ] == λ
    end
    @testset "update" begin
        x = 1.0
        STAT = EWMA(λ, value, 0.0, 1.0)
        update_statistic!(STAT, x)
        @test get_value(STAT) == λ * x
        lnew = 0.999
        set_design!(STAT, lnew)
        @test get_design(STAT)[1] == lnew
    end
end


@testset "CUSUM" begin
    k = 0.1; value = 0.0
    @testset "constructors" begin
        CUSUM(k = 0.2)
        STAT = CUSUM(k, value, true, 0.0, 1.0)
        @test_throws AssertionError CUSUM(k=-0.1)
        @test_throws AssertionError CUSUM(k=0.0)
        @test_throws AssertionError CUSUM(k=Inf)
        @test_throws AssertionError CUSUM(value=Inf)
        @test_throws AssertionError CUSUM(value=-Inf)
        params = get_design(STAT)
        @test length(params) == 1
        @test params[:k] == k
    end
    @testset "update" begin
        x = 1.0
        STAT = CUSUM(k, value, true, 0.0, 1.0)
        update_statistic!(STAT, x)
        @test get_value(STAT) == x - k
        knew = 0.999
        set_design!(STAT, knew)
        @test get_design(STAT)[1] == knew
    end
end


end#module