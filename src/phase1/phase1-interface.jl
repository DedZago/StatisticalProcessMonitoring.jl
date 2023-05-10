abstract type AbstractSampling end
export new_data
export new_data!

struct Bootstrap <: AbstractSampling end
export Bootstrap

#FIXME:test
new_data(B::Bootstrap, data::AbstractVector) = data[rand(1:length(data))]
new_data(B::Bootstrap, data::AbstractMatrix) = view(data, rand(1:size(data)[1]), :)
new_data!(B::Bootstrap, data::AbstractVecOrMat) = new_data(B, data) 

#FIXME:test
mutable struct BlockBootstrap{T} <: AbstractSampling
    block::T
    const blocksize::Int
    t::Int

    BlockBootstrap(data::Vector{T}, blocksize::Int) where T = new{Vector{T}}(zeros(T, blocksize), blocksize, 0)
    BlockBootstrap(data::Matrix{T}, blocksize::Int) where T = new{Matrix{T}}(zeros(T, blocksize, size(data)[2]), blocksize, 0)
end
export BlockBootstrap

#FIXME:test
new_data(B::BlockBootstrap, data::AbstractVector) = B.block[B.t]
new_data(B::BlockBootstrap, data::AbstractMatrix) = view(B.block, B.t, :)

#FIXME:implement and test
function update_block!(B::BlockBootstrap, data::AbstractVector)
    B.t = B.t % B.blocksize + 1
    if B.t == 1
        bb = sample(1:length(data))
        for i in 1:B.blocksize
            # wrap around the circle
            B.block[i] = data[(bb + i - 1) % length(data) + 1]
        end
    end
end

function update_block!(B::BlockBootstrap, data::AbstractMatrix)
    B.t = B.t % B.blocksize + 1
    if B.t == 1
        bb = sample(1:length(data))
        for i in 1:B.blocksize
            # wrap around the circle
            B.block[i] = data[(bb + i - 1) % length(data) + 1, :]
        end
    end
end

#FIXME:test
function new_data!(B::BlockBootstrap, data::AbstractVecOrMat)
    update_block!(B, data)
    return new_data(B, data) 
end

#FIXME:test
abstract type AbstractPhase1{S,T} end
# new_data(PH1::AbstractPhase1) = 
new_data!(::AbstractPhase1) = error("Not implemented for abstract interface.")

#FIXME:test
function shallow_copy_sim(PH1::T) where T <: AbstractPhase1
    return T(deepcopy(PH1.samp), PH1.data)
end

#FIXME:test
struct Phase1Data{S,T} <: AbstractPhase1{S,T} 
    samp::S
    data::T
end
export Phase1Data



"""
    new_data(P1::Phase1Data{S,T})
    new_data(P1::Phase1Data{S,AbstractVector})
    new_data(P1::Phase1Data{S,AbstractMatrix})

Generates a new observation based on the observed Phase I (in-control) data.
If it is not overloaded, then it defaults to generating data using a nonparametric bootstrap.
"""
#FIXME:test
new_data(P1::Phase1Data) = new_data(P1.samp, P1.data)
new_data!(P1::Phase1Data) = new_data!(P1.samp, P1.data)
export new_data!


################# TEST TRUE DGP #################
# struct Phase2Distribution{T} <: AbstractPhase1{T}
#     dist::T
# end
# export Phase2Distribution

# new_data(DGP::Phase2Distribution) = rand(DGP.dist)