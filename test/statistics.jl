module TestStatistics
using SPM
using Test

@testset "EWMA" begin
    λ = 0.1; value = 0.0
    @testset "constructors" begin
        EWMA(λ = 0.2)
        STAT = EWMA(λ, value)
        @test_throws AssertionError EWMA(-0.1, 0.0)
        @test_throws AssertionError EWMA(1.2, 0.0)
        @test_throws AssertionError EWMA(0.0, 0.0)
        @test_throws AssertionError EWMA(0.1, Inf)
        @test_throws AssertionError EWMA(0.1, -Inf)
        params = get_param(STAT)
        @test length(params) == 1
        @test params[:λ] == λ
    end
    @testset "update" begin
        x = 1.0
        STAT = EWMA(λ, value)
        update_statistic!(STAT, x)
        @test get_value(STAT) == λ * x
        lnew = 0.999
        set_param(STAT, lnew)
        @test get_param(STAT)[1] == lnew
    end
end


@testset "CUSUM" begin
    k = 0.1; value = 0.0
    @testset "constructors" begin
        CUSUM(k = 0.2)
        STAT = CUSUM(k, value, true)
        @test_throws AssertionError CUSUM(k=-0.1)
        @test_throws AssertionError CUSUM(k=0.0)
        @test_throws AssertionError CUSUM(k=Inf)
        @test_throws AssertionError CUSUM(value=Inf)
        @test_throws AssertionError CUSUM(value=-Inf)
        params = get_param(STAT)
        @test length(params) == 1
        @test params[:k] == k
    end
    @testset "update" begin
        x = 1.0
        STAT = CUSUM(k, value, true)
        update_statistic!(STAT, x)
        @test get_value(STAT) == x - k
        knew = 0.999
        set_param(STAT, knew)
        @test get_param(STAT)[1] == knew
    end
end


end#module