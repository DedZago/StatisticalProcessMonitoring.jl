
include("bootstrap.jl")

abstract type AbstractPhase2 end
new_data(::AbstractPhase2) = error("Not implemented for abstract interface.")
export new_data
new_data!(::AbstractPhase2) = error("Not implemented for abstract interface.")
export new_data!

function shallow_copy_sim(PH2::AbstractPhase2)
    return deepcopy(PH2)
end

"""
`Phase2` is a struct that holds the reference sample data and a sampling method, which is used to generate new observations from the reference data. 

### Arguments
- `samp::AbstractSampling`: The sampling method to be used to generate new observations. Defaults to `Bootstrap()`.
- `data`: The data from which observations need to be resampled.

### Examples
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
function shallow_copy_sim(PH2::Phase2)
    return Phase2(deepcopy(PH2.samp), PH2.data)
end


"""
    new_data(PH2::Phase2{S,T})
    new_data(PH2::Phase2{S,AbstractVector})
    new_data(PH2::Phase2{S,AbstractMatrix})

Generates a new observation using the Phase II object.
"""
new_data(PH2::Phase2) = new_data(get_sampler(PH2), get_data(PH2))
new_data!(PH2::Phase2) = new_data!(get_sampler(PH2), get_data(PH2))
export new_data!


################# TEST TRUE DGP #################

"""
    Phase2Distribution{T} <: AbstractPhase2

A struct that is used to generate and new data from a distribution. It contains a sampleable field `dist` of type `T`, which represents the underlying data-generating process.

### Notes
A method `rand(::T)` is required to generate new data from `dist`.

### Example
    using Distributions
    DGP = Phase2Distribution(Normal(0,1))
    new_data(DGP)
"""
struct Phase2Distribution{T} <: AbstractPhase2
    dist::T
end
export Phase2Distribution

new_data(PH2::Phase2Distribution) = rand(PH2.dist)
new_data!(PH2::Phase2Distribution) = rand(PH2.dist)

function shallow_copy_sim(PH2::Phase2Distribution) 
    return Phase2Distribution(deepcopy(PH2.dist))
end
