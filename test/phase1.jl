module TestPhase1
using SPM
using Test

@testset "Phase1Data" begin
    x = randn(15)
    xmat = randn(15, 5)
    @testset "Bootstrap" begin
        BB = Bootstrap()
        @test new_data(BB, x) in x
        @test new_data(BB, xmat) in eachrow(xmat)
        @test new_data!(BB, x) in x
        @test new_data!(BB, xmat) in eachrow(xmat)
    end
    @testset "Block bootstrap" begin
        blocksize = 10
        BB = BlockBootstrap(blocksize, x)
        @test get_blocksize(BB) == blocksize
        @test length(get_block(BB)) == blocksize
        @test BB.t == 1
        @test issubset(BB.block, x)
        y = new_data(BB, x)
        @test y in x
        @test y == get_block(BB)[1]
        BB = BlockBootstrap(blocksize, xmat)
        @test get_blocksize(BB) == blocksize
        @test size(get_block(BB)) == (blocksize, size(xmat)[2])
        @test BB.t == 1
        @test new_data(BB, xmat) in eachrow(xmat)
        BB = BlockBootstrap(1, xmat)
        @test get_block(BB)[1, :] in eachrow(xmat)
    end
    @testset "Phase 1" begin
        PH1 = Phase1Data(Bootstrap(), x)
        y = new_data(PH1)
        @test y in x
        PH1 = Phase1Data(Bootstrap(), xmat)
        y = new_data(PH1)
        @test y in eachrow(xmat)

        blocksize = 10
        PH1 = Phase1Data(BlockBootstrap(blocksize, x), x)
        y = new_data(PH1)
        @test y in x
        PH1 = Phase1Data(BlockBootstrap(blocksize, xmat), xmat)
        y = new_data(PH1)
        @test y in eachrow(xmat)
    end
end
end