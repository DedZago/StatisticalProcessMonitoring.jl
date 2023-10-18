using Parameters
using Distributions

"""
    RSADA{F,D}

A struct representing the RSADA control chart, whcih is used for monitoring the mean of partially-observed independent data streams.
The monitoring statistic iteratively samples arms (data streams) and updates a set of local statistics based on the observed values of the data streams.
The RSADA algorithm uses the local statistics to make decisions on which arms to sample at each iteration.

# Fields
- `k`: The allowance constant of the CUSUM Control chart. Default value is 0.1.
- `mu_min`: The minimum value of the mean shift to detect. Default value is 0.2.
- `p::Int`: The number of independent data streams. Default value is 1.
- `q::Int`: The number of data streams that are observable at each iteration. Default value is 1.
- `dist::Distribution`: The distribution of the data streams. Default value is `Normal(0,1)`.
- `value`: The current value of the monitoring statistic. Default value is 0.0.
- `eta`: The vector of augmented variables. Initialized to be a vector of zeros for each data stream.
- `obs::Vector{Int}`: The indices of the data streams to sample. Default value is an array of q random integers between 1 and p.
- `S1::Vector{F}`: The sum of the rewards for each arm.
- `S2::Vector{F}`: The sum of the squared rewards for each arm.

# References
Xian, X., Zhang, C., Bonk, S., & Liu, K. (2019). Online monitoring of big data streams: A rank-based sampling algorithm by data augmentation. Journal of Quality Technology, 53(2), 135-153. https://doi.org/10.1080/00224065.2019.1681924
"""
@with_kw mutable struct RSADA{F,D,S <: AbstractSampling} <: UnivariateStatistic
    k::F = 0.1
    mu_min::F = 0.2
    p::Int = 1
    q::Int = 1
    dist::D = Normal(0,1)
    value::F = 0.0
    eta::Vector{F} = zeros(p)
    obs::Vector{Int} = sample(1:p, q)
    S1::Vector{F} = zeros(p)
    S2::Vector{F} = zeros(p)
    g::Vector{F} = fill(1/p, p)
    sampler::S
end
export RSADA

RSADA(k, mu_min, q, x::AbstractMatrix; dist = Normal(0,1), sampler = ThompsonSampling()) = RSADA(k=k, mu_min=mu_min, p=size(x,2), q=q, dist=dist, sampler = sampler)

get_design(STAT::RSADA) = [STAT.k, STAT.mu_min]

function set_design!(STAT::RSADA, p::AbstractVector)
    STAT.k = p[1]
    STAT.mu_min = p[2]
    return p
end


"""
    augmented_vector(y, obs, mu, dist)

Compute the augmented vector `eta` (Equation 10 of [Xian et Al., 2019]) based on the given parameters.

# Arguments
- `y`: The input vector.
- `obs`: The indices of the observed elements in `y`.
- `mu`: The mean parameter.
- `dist`: The IC distribution of the individual data streams.

# Returns
- `eta`: The computed augmented vector.

# References
Xian, X., Zhang, C., Bonk, S., & Liu, K. (2019). Online monitoring of big data streams: A rank-based sampling algorithm by data augmentation. Journal of Quality Technology, 53(2), 135-153. https://doi.org/10.1080/00224065.2019.1681924
"""
function augmented_vector!(y, obs, mu, eta, dist)
    # According to Xian (2019)
    notObs = setdiff(1:length(y), obs)
    A = sum(pdf.(dist, y[obs] .- mu) ./ pdf.(dist, y[obs]))
    xi, i = findmax(y)
    B = cdf(dist, xi)
    C = cdf(dist, xi - mu)
    d = length(notObs)
    eta .= zeros(length(y))                          # Equation (8)
    eta[i] = (B^d * A + B^(d-1) * C * d) / (A + d)  # Equation (7)
    eta[notObs] .= (1.0 - eta[i]) / d               # Equation (9)
    return eta
end
export augmented_vector

"""
    update_sampling!(STAT::RSADA)

Update the monitored data streams of the RSADA monitoring statistic `STAT` by selecting the top-`q` values of the local monitoring statistics in `STAT.S1`.

# Arguments
- `STAT`: The RSADA object to update.

# Returns
- `STAT.obs`: The selected indices of the top `STAT.q` values of `STAT.S1`.
"""
function update_sampling!(STAT::RSADA)
    STAT.obs .= new_layout(STAT.sampler, STAT.S1, STAT.q)
    return STAT.obs
end


function update_statistic!(STAT::RSADA, x::AbstractVector)
    STAT.eta .= augmented_vector!(x, STAT.obs, STAT.mu_min, STAT.eta, STAT.dist)
    STAT.S1 .= STAT.S1 .+ STAT.eta
    STAT.S2 .= STAT.S2 .+ STAT.g
    C = sum((STAT.S1 - STAT.S2).^2 ./ STAT.S2)
    # Equation (10) of [Xian et Al., 2019]
    if C <= STAT.k
        STAT.S1[:] .= STAT.g[:]
        STAT.S2[:] .= STAT.g[:]
    else
        C = (C-STAT.k)/C
        STAT.S1 .*= C
        STAT.S2 .*= C
    end
    update_sampling!(STAT)
    STAT.value = sum((STAT.S1 .- STAT.S2).^2 ./ STAT.S2)
end