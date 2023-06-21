abstract type AbstractSampling end
export new_data
export new_data!

############## IID OBSERVATIONS ##############
struct Bootstrap <: AbstractSampling end
export Bootstrap

new_data(B::Bootstrap, data::AbstractVector) = data[rand(1:length(data))]
new_data(B::Bootstrap, data::AbstractMatrix) = view(data, rand(1:size(data)[1]), :)
new_data!(B::Bootstrap, data::AbstractVecOrMat) = new_data(B, data) 

############## TIME SERIES ##############
include("tsboot.jl")


abstract type AbstractPhase1 end
# new_data(PH1::AbstractPhase1) = 
new_data!(::AbstractPhase1) = error("Not implemented for abstract interface.")

function shallow_copy_sim(PH1::T) where T <: AbstractPhase1
    return T(deepcopy(PH1.samp), PH1.data)
end

struct Phase1Data{S,T} <: AbstractPhase1
    samp::S
    data::T
end
export Phase1Data

get_sampler(PH1::Phase1Data) = PH1.samp
export get_sampler
get_data(PH1::Phase1Data) = PH1.data
export get_data



"""
    new_data(P1::Phase1Data{S,T})
    new_data(P1::Phase1Data{S,AbstractVector})
    new_data(P1::Phase1Data{S,AbstractMatrix})

Generates a new observation based on the observed Phase I (in-control) data.
If it is not overloaded, then it defaults to generating data using a nonparametric bootstrap.
"""
new_data(P1::Phase1Data) = new_data(P1.samp, P1.data)
new_data!(P1::Phase1Data) = new_data!(P1.samp, P1.data)
export new_data!


################# TEST TRUE DGP #################
struct Phase2Distribution{T} <: AbstractPhase1
    dist::T
end
export Phase2Distribution

new_data(DGP::Phase2Distribution) = rand(DGP.dist)