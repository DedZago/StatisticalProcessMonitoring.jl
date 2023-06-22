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

#TODO: see if an abstraction is needed to avoid defining shallow_copy_sim for every type <: AbstractPhase1 that has to be defined.
function shallow_copy_sim(PH1::Phase2Distribution) 
    return Phase2Distribution(deepcopy(PH1.dist))
end

"""
`Phase1Data` is a struct that holds the reference sample data and a sampling method to generate new observations from the reference data. 

# Arguments
- `data`: The data obtained by phase 1.
- `samp = Bootstrap()`: The sampling method to be used to generate new observations. Defaults to `Bootstrap()`.

# Examples
x = randn(500)
PH1 = Phase1Data(data = x)
"""
@with_kw struct Phase1Data{S <: AbstractSampling, T} <: AbstractPhase1
    samp::S = Bootstrap()
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

"""
    Phase2Distribution{T} <: AbstractPhase1

A struct representing Phase II observations, it is used to generate and monitor new data from the true data-generating process. It contains a field `dist` of type `T`, which represents the underlying data-generating process.

# Notes
A method `rand(::T)` is required to generate new data from `dist`.

# Example
    DGP = Phase2Distribution(Normal(0,1))
    new_data(DGP)
"""
struct Phase2Distribution{T} <: AbstractPhase1
    dist::T
end
export Phase2Distribution

new_data(DGP::Phase2Distribution) = rand(DGP.dist)
new_data!(DGP::Phase2Distribution) = rand(DGP.dist)