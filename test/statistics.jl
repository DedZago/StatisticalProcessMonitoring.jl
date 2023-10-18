module TestStatistics
using SPM
using Test
using Random
using StatsBase
using Distributions
using GLM

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
        set_design!(STAT, [0.1])
        @test get_design(STAT) == [0.1]
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
        set_design!(STAT, [0.1])
        @test get_design(STAT) == [0.1]
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

@testset "Shewhart" begin
    S = Shewhart()
    @test get_design(S) == Vector{Float64}()
    @test_throws ErrorException set_design!(S, 1)
end

@testset "AEWMA" begin
    λ = 0.1; k = 2.0; value = 0.0
    @testset "constructors" begin
        AEWMA(λ = 0.2, k = 2.0)
        STAT = AEWMA(λ, k, value)
        @test_throws AssertionError AEWMA(λ=-0.1, k=1.0, value=0.0)
        @test_throws AssertionError AEWMA(λ=1.2, value=0.0)
        @test_throws AssertionError AEWMA(λ=0.1, value=Inf)
        @test_throws AssertionError AEWMA(λ=0.1, value=-Inf)
        @test_throws AssertionError AEWMA(λ=0.1, k=-0.1, value=0.0)
        @test_throws AssertionError AEWMA(λ=0.1, k=2.0,value=Inf)
        @test_throws AssertionError AEWMA(λ=0.1, k=3.0,value=-Inf)
        params = get_design(STAT)
        @test length(params) == 2
        @test params[1] == λ
        @test params[2] == k
        set_design!(STAT, [0.1, 5.0])
        @test get_design(STAT) == [0.1, 5.0]
    end
    @testset "update" begin
        x = 1.0
        STAT = AEWMA(λ, k, value)
        update_statistic!(STAT, x)
        @test get_value(STAT) == λ * x
        x = 10.0
        STAT = AEWMA(λ, k, value)
        update_statistic!(STAT, x)
        @test get_value(STAT) == x - (1-λ)*k
        lnew = [0.999, 100.0]
        set_design!(STAT, lnew)
        @test get_design(STAT) == lnew
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

        @testset "Univariate" begin
            STAT = OneSidedEWMA(λ = 0.2)   
            update_statistic!(STAT, randn())
            STAT = OneSidedEWMA(λ = 0.2, upw=false)   
            update_statistic!(STAT, randn())
        end

        @testset "Multivariate" begin
            Random.seed!(123)
            x = randn(100, 3)
            MS = MCUSUM(0.25, x)
            STAT = LocationScaleStatistic(MS, x)
            @test get_design(STAT) == get_design(MS)
            @test get_statistic(STAT) == MS 
            @test get_value(STAT) == get_value(MS)
            @test get_maxrl(STAT) == get_maxrl(MS)
            xnew = randn(3)
            znew = sqrt(inv(cov(x))) * (xnew - mean.(eachcol(x)))
            @test update_statistic(STAT, xnew) ≈ update_statistic(MS, znew)
            update_statistic!(STAT, xnew)
            set_value!(STAT, 0.2)
            @test get_value(STAT) == 0.2
            set_design!(STAT, 0.2)
            @test get_design(STAT) == [0.2]
        end

        @testset "other multivariate" begin
            STAT = MShewhart(randn(100,3))
            update_statistic!(STAT, randn(3))            
            update_statistic(STAT, randn(3))            
            @test get_design(STAT) == Vector{Float64}()
            @test_throws ErrorException set_design!(STAT, 3.0)

            STAT = DiagMEWMA(Λ = [0.2,0.2])
            set_design!(STAT, [0.1,0.1])
            @test get_design(STAT) == [0.1,0.1]
            set_design!(STAT, [0.5])
            @test get_design(STAT) == [0.5,0.5]
            update_statistic!(STAT, [1.0,1.0])
            get_value(STAT)

            STAT = AMCUSUM(0.1, randn(100,3))
            get_value(STAT), set_value!(STAT, 0.1), get_design(STAT), set_design!(STAT, 0.1), set_design!(STAT, [0.5])
            update_statistic!(STAT, randn(3))
        end

        @testset "data categorization" begin
            @testset "backward selection" begin
                Random.seed!(123)
                n = 50
                using DataFrames
                df = DataFrame(y = rand(Poisson(3.0), n), x1 = randn(n), x2 = randn(n))
                blmodel = backward_loglinear(df, :y)
                @test size(blmodel.mm) == (n, 1)
            end
            x = randn(1000,2)
            STAT = LLCUSUM(0.1, x) 
            PH2 = MultinomialBootstrap(STAT)
            update_statistic(STAT, new_data!(PH2, x))

            STAT = WANG2017(0.1, x) 
            PH2 = MultinomialBootstrap(STAT)
            update_statistic(STAT, new_data!(PH2, x))

            STAT = LI2012(0.1, x) 
            PH2 = MultinomialBootstrap(STAT)
            update_statistic(STAT, new_data!(PH2, x))
        end

        @testset "partially observed" begin
            NM = ARL(200)
            p = 5
            q = 2
            x = randn(1,p)
            STAT = RSADA(0.3, 1.0, q, x, sampler = ThompsonSampling())
            LIM = OneSidedFixedLimit(1.0, true)
            DIST = MvNormal(zeros(p), ones(p))
            PH2 = Phase2Distribution(DIST)
            CH = ControlChart(STAT, LIM, NM, PH2)
            update_chart!(CH, rand(DIST))

            STAT = RSADA(0.3, 1.0, q, x, sampler = TopQ())
            CH = ControlChart(STAT, LIM, NM, PH2)
            update_chart!(CH, rand(DIST))
        end

        @testset "covariance matrix" begin
            NM = ARL(200)
            p = 3
            STAT = MEWMC(λ=0.1, p = p)
            LIM = OneSidedFixedLimit(1.0, true)
            DIST = MvNormal(zeros(p), ones(p))
            PH2 = Phase2Distribution(DIST)
            CH = ControlChart(STAT, LIM, NM, PH2)
            update_chart(CH, rand(DIST))
        end

        @testset "Functional data" begin
            using GLM, DataFrames
            struct MyLinearModel
                a
                b
            end
            SPM.predict(mm::MyLinearModel, x::AbstractVector) = mm.a .+ mm.b .* x
            NM = ARL(200)
            n = 500
            nj = 10
            xs = 10 .* rand(n, nj)
            ys = 1.0 .+ 1.0 * xs .+ randn(n, nj)
            g = MyLinearModel(1.0, 1.0)
            dat = FunctionalData(xs, ys)
            STAT = NEWMA(0.2, g, dat)
            z = update_statistic(STAT, dat[1])
            @test 0 <= z < Inf
        end

    end
end

end#module