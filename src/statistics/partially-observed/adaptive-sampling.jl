using Parameters
using Distributions

"""
    AbstractSampling

Abstract type representing a generic sampling method for partially-observed data.

"""
abstract type AbstractSampling end
export AbstractSampling


"""
    new_layout(sampling::AbstractSampling, local_statistics) -> layout

Abstract function that creates a new layout based on the given sampling method and local statistics.

# Arguments
- `sampler::AbstractSampling`: An instance of an object subtype of `AbstractSampling`.
- `local_statistics`: Array of local statistics used for generating the layout.
- `q::Int`: The number of items to include in the layout.


# Returns
- `layout`: The generated layout to observe at the next iteration.

# Examples
```julia
# Create an instance of ThompsonSampling
ts = ThompsonSampling(2.0)

# Generate a layout using ThompsonSampling
stats = [1, 2, 3, 4]
layout = new_layout(ts, stats, 2)
```
"""
new_layout(sampler::AbstractSampling, local_statistics, q::Int) = error("Not defined for abstract class.")


"""
    ThompsonSampling <: AbstractSampling

Type representing the Thompson Sampling method.

# Fields
- `β::Float64`: A parameter regulating the concentration of the `Dirichlet` distribution for Thompson sampling.

# Examples
```julia
ts = ThompsonSampling(2.0)
```
"""
@with_kw struct ThompsonSampling <: AbstractSampling
    β::Float64 = 1.0
end
export ThompsonSampling

function new_layout(SPL::ThompsonSampling, local_statistics::AbstractVector, q::Int)
    # Use Thompson sampling from the Dirichlet distribution, local statistics should be non-negative 
    obs_new = rand(Dirichlet(1.0 .+ SPL.β .* local_statistics))
    return sortperm(obs_new, rev=true)[1:q]
end

"""
    TopQ <: AbstractSampling

Type representing the TopQ sampling method. New data is generating by taking the `q` local statistics with highest value. 

# References
Mei, Y. (2011). Quickest detection in censoring sensor networks. 2011 IEEE International Symposium on Information Theory Proceedings, 2148-2152. https://doi.org/10.1109/ISIT.2011.6034390
"""
struct TopQ <: AbstractSampling end
export TopQ

function new_layout(SPL::TopQ, local_statistics::AbstractVector, q::Int)
    # Add a small amount of random noise to remove ties in the local statistics
    return sortperm(local_statistics .+ 1e-10*randn(length(local_statistics)), rev=true)[1:q]
end