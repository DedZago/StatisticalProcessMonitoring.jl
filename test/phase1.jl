module TestPhase1
using SPM
using Test
using Random

@testset "Phase1Data" begin
    Random.seed!(123)
    x = randn(1000)
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
        @test get_counter(BB) == 1
        @test issubset(BB.block, x)
        y = new_data(BB, x)
        @test y in x
        @test y == get_block(BB)[1]
        BB = BlockBootstrap(blocksize, xmat)
        @test get_blocksize(BB) == blocksize
        @test size(get_block(BB)) == (blocksize, size(xmat)[2])
        @test get_counter(BB) == 1
        @test new_data(BB, xmat) in eachrow(xmat)
        BB = BlockBootstrap(1, xmat)
        @test get_block(BB)[1, :] in eachrow(xmat)
    end

    @testset "Block bootstrap generation" begin
        Random.seed!(123)
        blocksize = 5
        BB = BlockBootstrap(blocksize, x)
        initblock = deepcopy(x[1:blocksize])
        BB.block = deepcopy(initblock)
        @test get_block(BB) == initblock
        y = new_data!(BB, x)
        @test y == x[1]
        @test get_counter(BB) == 2
        @test get_block(BB) == initblock
        y = new_data!(BB, x)
        @test y == x[2]
        @test get_block(BB) == initblock
        @test get_counter(BB) == 3
        y = new_data!(BB, x)
        @test y == x[3]
        @test get_block(BB) == initblock
        @test get_counter(BB) == 4
        y = new_data!(BB, x)
        @test y == x[4]
        @test get_block(BB) == initblock
        @test get_counter(BB) == 5
        y = new_data!(BB, x)
        @test y == x[5]
        @test get_block(BB) != initblock
        y = new_data!(BB, x)
        @test get_counter(BB) == 2
        @test !(y in initblock)
        @test get_block(BB) != initblock
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