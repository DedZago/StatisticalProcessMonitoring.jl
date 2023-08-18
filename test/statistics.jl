module TestStatistics
using SPM
using Test
using Random
using StatsBase

@testset "EWMA" begin
    λ = 0.1; value = 0.0
    @testset "constructors" begin
        EWMA(λ = 0.2)
        STAT = EWMA(λ, value)
        @test_throws AssertionError EWMA(λ=-0.1, value=0.0)
        @test_throws AssertionError EWMA(λ=1.2, value=0.0)
        @test_throws AssertionError EWMA(λ=0.1, value=Inf)
        @test_throws AssertionError EWMA(λ=0.1, value=-Inf)
        params = get_design(STAT)
        @test length(params) == 1
        @test params[1] == λ
    end
    @testset "update" begin
        x = 1.0
        STAT = EWMA(λ, value)
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
        STAT = CUSUM(k, value, true)
        @test_throws AssertionError CUSUM(k=-0.1)
        @test_throws AssertionError CUSUM(k=0.0)
        @test_throws AssertionError CUSUM(k=Inf)
        @test_throws AssertionError CUSUM(value=Inf)
        @test_throws AssertionError CUSUM(value=-Inf)
        params = get_design(STAT)
        @test length(params) == 1
        @test params[1] == k
    end
    @testset "update" begin
        x = 1.0
        STAT = CUSUM(k, value, true)
        update_statistic!(STAT, x)
        @test get_value(STAT) == x - k
        knew = 0.999
        set_design!(STAT, knew)
        @test get_design(STAT)[1] == knew
    end
end

@testset "Estimated statistics" begin
    @testset "Location scale" begin
        @testset "EWMA" begin
            Random.seed!(123)
            x = randn(500)
            E = EWMA(λ = 0.2)
            STAT = LocationScaleStatistic(E, x)
            @test get_design(STAT) == get_design(E)
            @test get_statistic(STAT) == E 
            @test get_value(STAT) == get_value(E)
            @test get_maxrl(STAT) == get_maxrl(E)
            xnew = randn()
            znew = (xnew - mean(x)) / std(x)
            @test update_statistic(STAT, xnew) == update_statistic(E, znew)
        end

        @testset "CUSUM" begin
            Random.seed!(123)
            x = randn(500)
            E = CUSUM(k = 0.5)
            STAT = LocationScaleStatistic(E, x)
            @test get_design(STAT) == get_design(E)
            @test get_statistic(STAT) == E 
            @test get_value(STAT) == get_value(E)
            @test get_maxrl(STAT) == get_maxrl(E)
            xnew = randn()
            znew = (xnew - mean(x)) / std(x)
            @test update_statistic(STAT, xnew) == update_statistic(E, znew)
        end
    end
end

end#module