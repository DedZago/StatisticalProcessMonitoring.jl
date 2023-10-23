
include("bootstrap.jl")

abstract type AbstractPhase2 end
new_data(::AbstractPhase2) = error("Not implemented for abstract interface.")
export new_data
new_data!(::AbstractPhase2) = error("Not implemented for abstract interface.")
export new_data!

function shallow_copy_sim(PH2::T) where T <: AbstractPhase2
    return T(deepcopy(PH2.samp), PH2.data)
end

"""
`Phase2` is a struct that holds the reference sample data and a sampling method to generate new observations from the reference data. 

# Arguments
- `samp = Bootstrap()`: The sampling method to be used to generate new observations. Defaults to `Bootstrap()`.
- `data`: The data obtained after phase 1.

# Examples
x = randn(500)
PH2 = Phase2(data = x)
"""
@with_kw struct Phase2{S <: AbstractSampling, T} <: AbstractPhase2
    samp::S = Bootstrap()
    data::T
end
export Phase2

get_sampler(PH2::Phase2) = PH2.samp
export get_sampler
get_data(PH2::Phase2) = PH2.data
export get_data



"""
    new_data(PH2::Phase2{S,T})
    new_data(PH2::Phase2{S,AbstractVector})
    new_data(PH2::Phase2{S,AbstractMatrix})

Generates a new observation based on the observed Phase II (in-control) data.
"""
new_data(PH2::Phase2) = new_data(PH2.samp, PH2.data)
new_data!(PH2::Phase2) = new_data!(PH2.samp, PH2.data)
export new_data!


################# TEST TRUE DGP #################

"""
    Phase2Distribution{T} <: AbstractPhase2

A struct representing Phase II observations, it is used to generate and monitor new data from the true data-generating process. It contains a field `dist` of type `T`, which represents the underlying data-generating process.

# Notes
A method `rand(::T)` is required to generate new data from `dist`.

# Example
    DGP = Phase2Distribution(Normal(0,1))
    new_data(DGP)
"""
struct Phase2Distribution{T} <: AbstractPhase2
    dist::T
end
export Phase2Distribution

new_data(PH2::Phase2Distribution) = rand(PH2.dist)
new_data!(PH2::Phase2Distribution) = rand(PH2.dist)

#TODO: see if an abstraction is needed to avoid defining shallow_copy_sim for every type <: AbstractPhase2 that has to be defined.
function shallow_copy_sim(PH2::Phase2Distribution) 
    return Phase2Distribution(deepcopy(PH2.dist))
end
