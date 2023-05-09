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
struct BlockBootstrap{T} <: AbstractSampling
    block::T
    blocksize::Int
    t::Int

    BlockBootstrap(data::Vector{T}, blocksize::Int) where T = new{T}(zeros(T, blocksize), blocksize, 0)
    BlockBootstrap(data::Matrix{T}, blocksize::Int) where T = new{T}(zeros(T, blocksize, size(data)[2]), blocksize, 1)
end
export BlockBootstrap

#FIXME:test
new_data(B::BlockBootstrap, data::AbstractVector) = view(B.block, B.t)
new_data(B::BlockBootstrap, data::AbstractMatrix) = view(B.block, B.t, :)

#FIXME:implement and test
function update_block!(B::BlockBootstrap, data::AbstractVector)
    error("Not implemented yet.")
end
function update_block!(B::BlockBootstrap, data::AbstractMatrix)
    error("Not implemented yet.")
end

#FIXME:test
function new_data!(B::BlockBootstrap, data::Abstract)
    ret = new_data(B, data) 
    update_block!(B, data)
    return ret
end

#FIXME:test
abstract type AbstractPhase1{S,T} end
# new_data(PH1::AbstractPhase1) = 
new_data!(PH1::AbstractPhase1) = error("Not implemented for abstract interface.")

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