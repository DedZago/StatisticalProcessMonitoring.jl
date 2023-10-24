using LinearAlgebra

#############################################################################################################
#                                          LOCATION MONITORING                                              #
#############################################################################################################
"""
    MShewhart(value)

Shewhart control chart for monitoring multivariate observations with initial value `value`.

The update mechanism based on a new observation `x` is given by

``value = x'x``.

### References 
* Shewhart, W. A. (1931). Economic Control of Quality Of Manufactured Product. D. Van Nostrand Company.
"""
@with_kw mutable struct MShewhart <: AbstractStatistic 
    value::Float64 = 0.0
    @assert !isinf(value)
end
export MShewhart

get_design(::MShewhart) = Vector{Float64}()
set_design!(::MShewhart, par::Float64) = error("Cannot set a design for Shewhart chart.")
set_design!(::MShewhart, par::AbstractVector) = error("Cannot set a design for Shewhart chart.")

update_statistic(stat::MShewhart, x::AbstractVector) = x'x
update_statistic!(stat::MShewhart, x::AbstractVector) = stat.value = update_statistic(stat, x)

"""
    DiagMEWMA(Λ, value)

Exponentially weighted moving average with diagonal smoothing matrix `Λ` and initial value `value`.

The update mechanism based on a new observation `x` is given by

``Z_t = (I-Λ)Z_{t-1} + Λ x_t``,

and the chart value is defined as

``value_t = Z_t' Λ^{-1} Z_t``.

### Arguments
- `Λ::Vector{Float64}`: Vector of smoothing constants.
- `value::Float64`: Current value of the statistic (default = 0.0).
- `z::Vector{Float64}`: Vector of smoothed observations (Default: `zeros(length(Λ))`).
- `inv_Σz::Matrix{Float64}`: Inverse of the covariance matrix of the control variates.


### References 
Lowry, C. A., Woodall, W. H., Champ, C. W., & Rigdon, S. E. (1992). A Multivariate Exponentially Weighted Moving Average Control Chart. Technometrics, 34(1), 46-53. https://doi.org/10.1080/00401706.1992.10485232
"""
@with_kw mutable struct DiagMEWMA <: AbstractStatistic 
    Λ::Vector{Float64}
    value::Float64 = 0.0
    z::Vector{Float64} = zeros(length(Λ))
    inv_Σz::Matrix{Float64} = inv(diagm([Λ[k]^2/(2*Λ[k]-Λ[k]^2) for k in 1:length(Λ)]))
    @assert !isinf(value)
end
export DiagMEWMA


get_design(stat::DiagMEWMA) = deepcopy(stat.Λ)

function set_design!(stat::DiagMEWMA, Λ::AbstractVector)
    Λ_v = zeros(length(stat.z))
    if length(Λ) == 1 && length(stat.z) > 1
        Λ_v = fill(first(Λ), length(stat.z))
    else
        Λ_v[:] = Λ[:]
    end
    stat.inv_Σz = inv(diagm(Λ_v)*diagm(Λ_v))
    stat.Λ = Λ_v
end

function update_statistic!(stat::DiagMEWMA, x::AbstractVector)
    for j in eachindex(x)
        stat.z[j] = (1.0 - stat.Λ[j]) * stat.z[j] + stat.Λ[j] * x[j]
    end  
    stat.value = stat.z' * stat.inv_Σz * stat.z
end


"""
    MCUSUM(k, p, value = 0.0, St = zeros(p))

A mutable struct representing a Multivariate Cumulative Sum (MCUSUM) statistic.

### Arguments
- `k::Float64`: The value of the allowance parameter.
- `p::Int`: The number of variables to be monitored. 
- `value::Float64`: The initial value of the statistic. Default is 0.0.
- `St::Vector{Float64}`: A vector representing the multivariate cumulative sum at the current time `t`.

### Examples
    stat = MCUSUM(0.25, 2, 0.0, [0.0, 0.0])

### References
- Crosier, R. B. (1988). Multivariate Generalizations of Cumulative Sum Quality-Control Schemes. Technometrics, 30(3), 291-303. https://doi.org/10.2307/1270083
"""
@with_kw mutable struct MCUSUM <: AbstractStatistic 
    k::Float64
    p::Int
    value::Float64 = 0.0
    St::Vector{Float64} = zeros(p)
    @assert !isinf(value)
    @assert p > 0
    @assert k > 0
    @assert p == length(St)
end
export MCUSUM

MCUSUM(k::Real, x::AbstractMatrix) = MCUSUM(k=k, p=size(x,2))

get_design(stat::MCUSUM) = [stat.k]
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
    AMCUSUM(λ, p; minshift = 0.1, shift = 0.0, Et = zeros(p), t = 0, stat = MCUSUM(k=0.1, p=p))

A mutable struct representing an Adaptive Multivariate Cumulative Sum (MCUSUM) statistic.

### Arguments
- `λ::Float64`: The value of λ, where 0.0 <= λ <= 1.0.
- `p::Int`: The number of quality variables to monitor.
- `minshift::Float64`: The minimum shift value to be detected. Default is 0.1.
- `shift::Float64`: The current shift value. Default is 0.0.
- `Et::Vector{Float64}`: The vector Et of smoothed deviations from the zero mean. Has to be exactly equal to `zeros(p)`
- `t::Int`: The current value of time. Has to be exactly 0.
- `stat::C`: The underlying classical MCUSUM statistic. Default is MCUSUM(k=0.1, p=p).

### References
- Dai, Y., Luo, Y., Li, Z., & Wang, Z. (2011). A new adaptive CUSUM control chart for detecting the multivariate process mean. Quality and Reliability Engineering International, 27(7), 877-884. https://doi.org/10.1002/qre.1177
"""
@with_kw mutable struct AMCUSUM{C} <: AbstractStatistic 
    λ::Float64
    p::Int
    minshift::Float64 = 0.1
    shift::Float64 = 0.0
    Et::Vector{Float64} = zeros(p)
    t::Int = 0
    stat::C = MCUSUM(k=0.1, p=p)

    @assert 0.0 <= λ <= 1.0
    @assert minshift > 0.0
    @assert t == 0
    @assert Et == zeros(p)
end
export AMCUSUM

AMCUSUM(λ::Real, x::AbstractMatrix; minshift=0.1) = AMCUSUM(λ=λ, p=size(x,2), minshift=minshift)

get_value(stat::AMCUSUM) = get_value(stat.stat)
set_value!(stat::AMCUSUM, x) = set_value!(stat.stat, x)

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


"""
    MAEWMA(λ, k, value, z::Vector{Float64})

Multivariate Adaptive Exponentially Weighted Moving Average control chart.

The update mechanism based on a new observation `x` is given by

``Z_t = (I-Ω)*Z_{t-1} + Ω * x_t``,

where Ω = ω(e)*I is an adaptive generalization of the classical MEWMA smoothing matrix.
The chart value is defined as

``value_t = Z_t' Z_t``.

### Arguments
- `λ::Float64`: The value of the EWMA smoothing parameter.
- `k::Float64`: The value of the parameter of the Huber score.
- `value::Float64`: The value of the statistic. (default: `0.0`)
- `z::Vector{Float64}`: The vector of smoothed observations.

### References 
Mahmoud, M. A., & Zahran, A. R. (2010). A Multivariate Adaptive Exponentially Weighted Moving Average Control Chart. Communications in Statistics - Theory and Methods, 39(4), 606-625. https://doi.org/10.1080/03610920902755813
"""
@with_kw mutable struct MAEWMA <: AbstractStatistic 
    λ::Float64
    k::Float64
    value::Float64 = 0.0
    z::Vector{Float64}
    @assert !isinf(value)
end
export MAEWMA

function set_design!(stat::MAEWMA, par::AbstractVector)
    stat.λ = par[1]    
    stat.k = par[2]    
    return par
end
get_design(stat::MAEWMA) = [stat.λ, stat.k]


function update_statistic!(stat::MAEWMA, x::AbstractVector)
    e = norm(dot(x - stat.z, x - stat.z)) + eps()
    @show e
    omega = huber(e, stat.λ, stat.k)/e
    @show omega
    stat.z = (I - omega*I)*stat.z + omega*I*x
    @show stat.z
    stat.value = dot(stat.z, stat.z)
    return stat.value
end


#############################################################################################################
#                                        COVARIANCE MONITORING                                              #
#############################################################################################################
"""
    ALT{M}

Generalized Likelihood Ratio statistic for monitoring changes in the variance-covariance matrix introduced by [Alt].

### Fields
- `value::Float64`: The value of the statistic, initialized to 0.0.
- `Ω::M`: The precision matrix of the in-control process.
- `detΣ::Float64`: The determinant of the in-control process variance.

### References
Alt, F. A. (1984). Multivariate quality control. In N. L. Johnson, S. Kotz, & C. R. Read (Eds.), The encyclopedia of statistical sciences (Vol. 6, pp. 110-122). Wiley.
"""
@with_kw mutable struct ALT{M} <: AbstractStatistic
    value::Float64 = 0.0
    Ω::M
    detΣ::Float64 = 1.0/det(Ω)
end
export ALT

ALT(x::AbstractVector) = ALT(value = 0.0, Ω = 1.0/cov(x))
ALT(x::AbstractMatrix) = ALT(value = 0.0, Ω = inv(cov(x)))

get_design(stat::ALT) = []

function update_statistic!(stat::ALT, x::AbstractVector)
    m = length(x)
    @assert m > 1 "Must have at least two observations"
    S_i = var(x)
    stat.value = -(m-1)*(log(S_i) - log(stat.detΣ) - stat.Ω*S_i) 
end

function update_statistic!(stat::ALT, x::AbstractMatrix)
    m,p = size(x)
    @assert m > 1 "Must have at least two observations"
    @assert m >= p "Must have at least as many observations ($(m)) as the number of variables ($(p))"
    S_i = cov(x)
    stat.value = -(m-1)*(p + log(det(S_i)) - log(stat.detΣ) - tr(stat.Ω*S_i)) 
end


update_statistic(stat::ALT, x::Real) = update_statistic!(deepcopy(stat), x)


"""
    MEWMC(λ, value)

Exponentially weighted moving covariance matrix control chart with smoothing constant `λ` and initial value `value`.

The update mechanism based on a new observation `x \\in \\mathbb{R}^p` is given by

``Z_t = (1 - λ)*Z_{t-1} + λ \\cdot xx'``,

and the chart value is defined as

``value_t = \\text{tr}(Z_t) - \\log|Z_t| - p``.

### References 
Hawkins, D. M., & Maboudou-Tchao, E. M. (2008). Multivariate Exponentially Weighted Moving Covariance Matrix. Technometrics, 50(2), 155-166.
"""
@with_kw mutable struct MEWMC{M} <: AbstractStatistic 
    λ::Float64 = 0.1
    p::Int
    value::Float64 = 0.0
    z::M = diagm(ones(p))
    @assert 0.0 < λ <= 1.0
    @assert p > 0
    @assert !isinf(value)
end
export MEWMC


get_design(stat::MEWMC) = [stat.λ]
set_design!(stat::MEWMC, λ::Float64) = stat.λ = λ
set_design!(stat::MEWMC, λ::AbstractVector) = stat.λ = first(λ)


function update_statistic!(stat::MEWMC, x::AbstractVector)
    stat.z = (1 - stat.λ) * stat.z + stat.λ * x*x'
    stat.value = tr(stat.z) - log(abs(det(stat.z))) - stat.p
    return stat.value
end

update_statistic(stat::MEWMC, x::Real) = update_statistic!(deepcopy(stat), x)


"""
    MEWMS(λ, value)

Exponentially weighted moving covariance matrix control chart with smoothing constant `λ`.

The update mechanism based on a new observation `x \\in \\mathbb{R}^p` is given by

``Z_t = (1 - λ)*Z_{t-1} + λ \\cdot xx'``,

and the chart value is defined as

``value_t = \\text{tr}(Z_t).

### References 
Huwang, L., Yeh, A. B., & Wu, C.-W. (2007). Monitoring Multivariate Process Variability for Individual Observations. Journal of Quality Technology, 39(3), 258-278. https://doi.org/10.1080/00224065.2007.11917692
"""
@with_kw mutable struct MEWMS{M} <: AbstractStatistic 
    λ::Float64 = 0.1
    value::Float64 = 0.0
    z::M
    @assert 0.0 < λ < 1.0
    @assert !isinf(value)
end
export MEWMS

MEWMS(λ, x::Real) = MEWMS(λ=λ, value=0.0, z=x^2) 
MEWMS(λ, x::AbstractVector) = MEWMS(λ=λ, value=0.0, z=x'x) 

function MEWMS(λ, x::AbstractMatrix)
    n,p = size(x)
    Σ = Matrix{Float64}(undef, p, p)
    if n > 1
        Σ[:] = cov(x)
    else
        Σ[:] = x'x
    end
    MEWMS(λ=λ, value=0.0, z=Σ) 
end

get_design(stat::MEWMS) = [stat.λ]
set_design!(stat::MEWMS, λ::Float64) = stat.λ = λ
set_design!(stat::MEWMS, λ::AbstractVector) = stat.λ = first(λ)

function update_statistic!(stat::MEWMS, x::Real)
    stat.z = (1 - stat.λ) * stat.z + stat.λ * x^2
    stat.value = tr(stat.z)
    return stat.value
end

function update_statistic!(stat::MEWMS, x::AbstractVector)
    stat.z = (1 - stat.λ) * stat.z + stat.λ * x*x'
    stat.value = tr(stat.z)
    return stat.value
end

update_statistic(stat::MEWMS, x::Real) = update_statistic!(deepcopy(stat), x)
