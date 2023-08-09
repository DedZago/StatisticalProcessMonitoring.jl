using LinearAlgebra

"""
    MShewhart(μ, Σ, value)

Multivariate Shewhart control chart for observations with mean `μ`, variance-covariance `Σ` and initial value `value`.

The update mechanism based on a new observation `x` is given by

``value = (x - μ) Σ^(-1) (x - μ)``.

### References 
* Shewhart, W. A. (1931). Economic Control of Quality Of Manufactured Product. D. Van Nostrand Company.
"""
@with_kw mutable struct MShewhart{VF,VM,V} <: UnivariateStatistic 
    μ::VF
    Σ_m1::VM
    value::V = 0.0
    @assert length(μ) == size(Σ_m1, 1)
    @assert size(Σ_m1, 1) == size(Σ_m1, 2)
    @assert !isinf(value)
end
export MShewhart

#TODO: test MShewhart control chart
MShewhart(x::AbstractMatrix; value = 0.0) = MShewhart(mean.(eachcol(x)), inv(cov(x)), value)

get_design(stat::MShewhart) = Vector{Float64}()
set_design!(stat::MShewhart, ::Float64) = error("Cannot set a design for Shewhart chart.")

update_statistic(stat::MShewhart, x::AbstractVector) = dot(x - stat.μ, stat.Σ_m1, x - stat.μ)
update_statistic!(stat::MShewhart, x::AbstractVector) = stat.value = update_statistic(stat, x)

"""
    MEWMA(Λ, value)

Exponentially weighted moving average with smoothing matrix `Λ` and initial value `value`.

The update mechanism based on a new observation `x` is given by

``Z = (I-Λ)*Z + Λ * x``,

and the chart value is defined as

``value = Z' Λ^(-1) Z``.

### References 
Lowry, C. A., Woodall, W. H., Champ, C. W., & Rigdon, S. E. (1992). A Multivariate Exponentially Weighted Moving Average Control Chart. Technometrics, 34(1), 46-53. https://doi.org/10.1080/00401706.1992.10485232
"""
@with_kw mutable struct DiagMEWMA{L,V} <: UnivariateStatistic 
    Λ::Vector{L}
    value::V = 0.0
    z::Vector{L} = zeros(length(Λ))
    inv_Σz::Matrix{L} = inv(diagm([Λ[k]^2/(2*Λ[k]-Λ[k]^2) for k in 1:length(Λ)]))
    @assert !isinf(value)
end
export DiagMEWMA


get_design(stat::DiagMEWMA) = deepcopy(stat.Λ)

function set_design!(stat::DiagMEWMA, Λ::AbstractVector)
    stat.inv_Σz = inv(diagm(Λ)*diagm(Λ))
    stat.Λ = Λ
end

function update_statistic!(stat::DiagMEWMA, x::AbstractVector)
    for j in eachindex(x)
        stat.z[j] = (1.0 - stat.Λ[j]) * stat.z[j] + stat.Λ[j] * x[j]
    end  
    stat.value = stat.z' * stat.inv_Σz * stat.z
end


"""
    MCUSUM{L,V}(k::L, value::V = 0.0, St::Vector{L})

A mutable struct representing a Multivariate Cumulative Sum (MCUSUM) statistic.

# Arguments
- `k::L`: The value of the allowance parameter.
- `value::V`: The initial value of the statistic. Default is 0.0.
- `St::Vector{L}`: A vector representing the multivariate cumulative sum at the current time `t`.

# Examples
    stat = MCUSUM(0.25, 0.0, [0.0, 0.0])

# References
- Crosier, R. B. (1988). Multivariate Generalizations of Cumulative Sum Quality-Control Schemes. Technometrics, 30(3), 291-303. https://doi.org/10.2307/1270083
"""
@with_kw mutable struct MCUSUM{L,V} <: UnivariateStatistic 
    k::L
    p::Int
    value::V = 0.0
    St::Vector{L} = zeros(p)
    @assert !isinf(value)
    @assert p > 0
    @assert k > 0
    @assert p == length(St)
end
export MCUSUM

MCUSUM(x::AbstractMatrix, k::Real) = MCUSUM(k=k, p=size(x,2))

#TODO: Test MCUSUM control chart 
get_design(stat::MCUSUM) = deepcopy(stat.k)

set_design!(stat::MCUSUM, k::Real) = stat.k = k
set_design!(stat::MCUSUM, k::AbstractVector) = stat.k = first(k)

function update_statistic!(stat::MCUSUM, x::AbstractVector)
    Ct = sqrt(dot(stat.St + x, stat.St + x))
    if Ct <= stat.k
        stat.St .= 0.0
    else
        stat.St = (1.0 - stat.k/Ct) * (stat.St + x)
    end
    stat.value = sqrt(dot(stat.St, stat.St))
    return stat.value
end


"""
    AMCUSUM{C,L,V}(λ::L, p::Int; minshift::L = 0.1, shift::L = 0.0, Et::Vector{V}, t::Int = 0, stat::C = MCUSUM(k=0.1, p=p))

Constructs an instance of the AMCUSUM struct.

# Arguments
- `λ::L`: The value of λ, where 0.0 <= λ <= 1.0.
- `p::Int`: The number of quality variables to monitor.
- `minshift::L`: The minimum shift value to be detected. Default is 0.1.
- `shift::L`: The current shift value. Default is 0.0.
- `Et::Vector{V}`: The vector Et of smoothed deviations from the zero mean. Has to be exactly equal to `zeros(p)`
- `t::Int`: The current value of time. Has to be exactly 0.
- `stat::C`: The underlying classical MCUSUM statistic. Default is MCUSUM(k=0.1, p=p).

# References
- Dai, Y., Luo, Y., Li, Z., & Wang, Z. (2011). A new adaptive CUSUM control chart for detecting the multivariate process mean. Quality and Reliability Engineering International, 27(7), 877-884. https://doi.org/10.1002/qre.1177
"""
@with_kw mutable struct AMCUSUM{C,L,V} <: UnivariateStatistic 
    λ::L
    p::Int
    minshift::L = 0.1
    shift::L = 0.0
    Et::Vector{V} = zeros(p)
    t::Int = 0
    stat::C = MCUSUM(k=0.1, p=p)

    @assert 0.0 <= λ <= 1.0
    @assert minshift > 0.0
    @assert t == 0
    @assert Et == zeros(p)
end
export AMCUSUM

AMCUSUM(x::AbstractMatrix, λ::Real; minshift=0.1) = AMCUSUM(λ=λ, p=size(x,2))

get_value(stat::AMCUSUM) = get_value(stat.stat)
set_value!(stat::AMCUSUM) = set_value!(stat.stat)

set_design!(stat::AMCUSUM, λ::Real) = stat.λ = λ
set_design!(stat::AMCUSUM, λ::AbstractVector) = stat.λ = first(λ)
get_design(stat::AMCUSUM) = [stat.λ]

function update_statistic!(stat::AMCUSUM, x::AbstractVector)
    update_statistic!(stat.stat, x)
    stat.t += 1
    λ = stat.λ
    stat.Et = (1 - λ) * stat.Et + λ * x
    squared_shift_est = 1/(1 - (1 - λ)^stat.t)^2 * (stat.Et'stat.Et - (1-(1-λ)^(2*stat.t))*(λ*length(stat.Et)/(2-λ)))
    squared_shift = max(stat.minshift, (1-λ)*(stat.shift)^2 + λ*squared_shift_est)
    stat.shift = sqrt(squared_shift)
    set_design!(stat.stat, stat.shift/2)
    return get_value(stat)
end
