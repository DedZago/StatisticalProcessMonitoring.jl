using Parameters
using Distributions

"""
    TRAS{F,D}

# Fields
A Top-r-based Adaptive Sampling monitoring statistic for monitoring multiple data streams with partial observations.

## Fields
- `k::Float64`: The allowance constant for the CUSUM statistic. Defaults to `0.1`.
- `mu_min::Float64`: The minimum shift to be detected. Defaults to `1.0`.
- `value::Float64`: The current value of the monitoring statistic. Defaults to `0.0`.
- `Δ::Float64`: The compensation coefficient for non-observed variables.
- `r::Int`: The number of largest variables to sum.
- `q::Int`: The total number of observable variables.
- `p::Int`: The total number of variables.
- `W::Vector{Float64}`: The vector of local monitoring statistics. Defaults to a vector of zeros.
- `W1::Vector{Float64}`: The vector used to store the local monitoring statistics to detect increases in the mean. Defaults to a vector of zeros.
- `W2::Vector{Float64}`: The vector used to store the local monitoring statistics to detect decreases in the mean. Defaults to a vector of zeros.
- `obs::Vector{Int}`: The selected indices of the top `q` values of local monitoring statistics. Defaults to a random sample of size `q` from the range `1:p`.
- `sampler::S`: The sampling strategy used to select the next set of `obs` indices. Defaults to `ThompsonSampling()`.

# References
Liu, K., Mei, Y., & Shi, J. (2015). An Adaptive Sampling Strategy for Online High-Dimensional Process Monitoring. Technometrics, 57(3), 305-319. https://doi.org/10.1080/00401706.2014.947005
"""
@with_kw mutable struct TRAS{S <: AbstractSampling} <: AbstractStatistic
    k::Float64 = 0.1
    mu_min::Float64 = 1.0
    value::Float64 = 0.0
    Δ::Float64
    r::Int
    q::Int
    p::Int
    W::Vector{Float64} = zeros(p)
    W1::Vector{Float64} = zeros(p)
    W2::Vector{Float64} = zeros(p)
    obs::Vector{Int} = sample(1:p, q)
    sampler::S = ThompsonSampling()
    @assert 1 <= r < q  "Number of variables to sum ($(r)) must be less than total number of observable variables ($(q))"
    @assert 1 <= q <= p  "Number of observable variables ($(q)) must be less than total number of variables ($(p))"
    @assert k > 0 "Allowance constant must be positive"
    @assert mu_min > 0 "Minimum shift to be detected must be positive"
    @assert Δ > 0 "Compensation coefficient must be positive"
end
export TRAS

get_design(STAT::TRAS) = [STAT.k, STAT.mu_min, STAT.Δ, STAT.r]

function set_design!(STAT::TRAS, p::AbstractVector)
    STAT.k = p[1]
    STAT.mu_min = p[2]
    STAT.Δ = p[3]
    STAT.r = p[4]
    return p
end


"""
    update_sampling!(STAT::TRAS)

Update the monitored data streams of the TRAS monitoring statistic `STAT` using its associated `AbstractSampler` object.

# Arguments
- `STAT`: The TRAS object to update.

# Returns
- `STAT.obs`: The selected indices of the top `STAT.q` values of `STAT.S1`.
"""
function update_sampling!(STAT::TRAS)
    STAT.obs .= new_layout(STAT.sampler, STAT.W, STAT.q)
    return STAT.obs
end


function update_statistic!(STAT::TRAS, x::AbstractVector)
    STAT.W1[STAT.obs] = max.(STAT.W1[STAT.obs] .+ STAT.mu_min * x[STAT.obs] .- STAT.mu_min^2/2, 0.0)
    STAT.W2[STAT.obs] = max.(STAT.W1[STAT.obs] .- STAT.mu_min * x[STAT.obs] .- STAT.mu_min^2/2, 0.0)
    notObs = setdiff(1:STAT.p, STAT.obs)
    STAT.W1[notObs] = STAT.W1[notObs] .+ STAT.Δ
    STAT.W2[notObs] = STAT.W2[notObs] .+ STAT.Δ
    STAT.W[:] = max.(STAT.W1, STAT.W2)

    update_sampling!(STAT)
    # Value of the statistic is the sum of top-r elements
    STAT.value = 0.0
    for i in 1:STAT.r
        STAT.value += STAT.W[i]
    end
    return STAT.value
end
