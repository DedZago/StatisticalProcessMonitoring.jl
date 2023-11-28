module TestNominal
using StatisticalProcessMonitoring
using Test

@testset "NominalProperties" begin
    @testset "constructors" begin
        NM = ARL(200)
        NM = QRL(200, 0.5)
        NM = ARL(200.0)
        NM = QRL(200.0, 0.3)
        @test_throws AssertionError ARL(-1.0)
        @test_throws AssertionError QRL(-1.0, 0.5)
        @test_throws AssertionError QRL(200.0, 0.0)
        @test_throws AssertionError QRL(200.0, -0.1)
        @test_throws AssertionError QRL(200.0, 1.0)
        @test_throws AssertionError QRL(200.0, 1.1)
    end
end
end
