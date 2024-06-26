"""
    Shewhart(value)

Shewhart control chart.

The update mechanism of ``C_t`` based on a new observation `x` is given by

``C_t = x``.

### References 
Shewhart, W. A. (1931). Economic Control of Quality Of Manufactured Product. D. Van Nostrand Company.
"""
@with_kw mutable struct Shewhart{V} <: AbstractStatistic 
    value::V = 0.0
    @assert !isinf(value)
end
export Shewhart

get_design(::Shewhart) = Vector{Float64}()
set_design!(::Shewhart, x) = error("Cannot set a design for Shewhart chart.")

update_statistic(::Shewhart, x::Real) = x
update_statistic!(stat::Shewhart, x::Real) = stat.value = update_statistic(stat, x)



"""
    EWMA(λ, value)

Exponentially weighted moving average with design parameter `λ` and initial value `value`.

The update mechanism for the statistic ``C_t`` based on a new observation `x` is given by

``C_t = (1-λ)\\cdot C_{t-1} + λ \\cdot x``.

### Arguments
- `λ::Float64`: The smoothing constant. Defaults to `0.1`.
- `value::Float64`: The initial value for the EWMA statistic. Defaults to `0.0`.

### References 
Roberts, S. W. (1959). Control Chart Tests Based on Geometric Moving Averages. Technometrics, 1(3), 239-250. https://doi.org/10.1080/00401706.1959.10489860
"""
@with_kw mutable struct EWMA <: AbstractStatistic 
    λ::Float64 = 0.1
    value::Float64 = 0.0
    @assert 0.0 < λ <= 1.0
    @assert !isinf(value)
end
export EWMA


get_design(stat::EWMA) = [stat.λ]
set_design!(stat::EWMA, λ::Float64) = stat.λ = λ
set_design!(stat::EWMA, λ::AbstractVector) = stat.λ = first(λ)


update_statistic(stat::EWMA, x::Real) = (1.0 - stat.λ) * stat.value + stat.λ * x
update_statistic!(stat::EWMA, x::Real) = stat.value = update_statistic(stat, x)

"""
    OneSidedEWMA(λ, value, upw::Bool)

OneSidedEWMA statistic with design parameter `λ` and initial value `value`.

The update mechanism based on a new observation `x` is given by:
* if `upw == true`, then ``C_t = \\max\\{0, (1-λ)\\cdot C_{t-1} + λ\\cdot x\\}``;
* if `upw == true`, then ``C_t = \\min\\{0, (1-λ)\\cdot C_{t-1} + λ\\cdot x\\}``;

### Arguments
- `λ::Float64`: The smoothing constant. Default is `0.1`.
- `value::Float64`: The current value of the statistic. Default is `0.0`.
- `upw::Bool`: Whether the statistic should monitor increases (`true`) or decrease (`falses`) in the process mean. Default is `true`.
"""
@with_kw mutable struct OneSidedEWMA <: AbstractStatistic 
    λ::Float64 = 0.1
    value::Float64 = 0.0
    upw::Bool = true
    @assert !isinf(value)
    @assert 0.0 < λ <= 1.0
end
export OneSidedEWMA


get_design(stat::OneSidedEWMA) = [stat.λ]
set_design!(stat::OneSidedEWMA, λ::Float64) = stat.λ = λ
set_design!(stat::OneSidedEWMA, λ::AbstractVector) = stat.λ = first(λ)

function update_statistic(stat::OneSidedEWMA, x::Real)
    if stat.upw
        return max(0.0, (1.0 - stat.λ) * stat.value + stat.λ*x)
    else
        return min(0.0, (1.0 - stat.λ) * stat.value + stat.λ*x)
    end
end

update_statistic!(stat::OneSidedEWMA, x::Real) = stat.value = update_statistic(stat, x)



"""
    CUSUM(k, value, upw::Bool)

CUSUM statistic with design parameter `k` and initial value `value`.

The update mechanism based on a new observation `x` is given by:
* if `upw == true`, then ``C_t = \\max\\{0, C_{t-1} + x - k\\}``;
* if `upw == false`, then ``C_t = \\min\\{0, C_{t-1} + x + k\\}``.

### Arguments
- `k::Float64`: The allowance constant of the CUSUM statistic. Defaults to `1.0`.
- `value::Float64`: The current value of the cumulative sum. Defaults to `0.0`.
- `upw::Bool`: A boolean indicating whether the CUSUM statistic is increasing or decreasing. Defaults to `true`.

### References 
Page, E. S. (1954). Continuous Inspection Schemes. Biometrika, 41(1/2), 100. https://doi.org/10.2307/2333009
"""
@with_kw mutable struct CUSUM <: AbstractStatistic 
    k::Float64 = 1.0
    value::Float64 = 0.0
    upw::Bool = true
    @assert !isinf(value)
    @assert (k > 0.0) && !isinf(k)
end
export CUSUM


get_design(stat::CUSUM) = [stat.k]
set_design!(stat::CUSUM, k::Float64) = stat.k = k
set_design!(stat::CUSUM, k::AbstractVector) = stat.k = first(k)


function update_statistic(stat::CUSUM, x::Real)
    if stat.upw
        return max(0.0, stat.value + x - stat.k)
    else
        return min(0.0, stat.value + x + stat.k)
    end
end


update_statistic!(stat::CUSUM, x::Real) = stat.value = update_statistic(stat, x)


"""
    WCUSUM <: AbstractStatistic

A weighted cumulative sum statistic.

# Arguments
- `k::Float64`: The allowance constant of the CUSUM used for calculating the WCUSUM statistic. Default is 1.0.
- `λ::Float64`: The smoothing constant for updating the estimate of the fault signature. Must be a value between 0 and 1. Default is 0.2.
- `value::Float64`: The initial value of the weighted cumulative sum statistic. Default is 0.0.
- `Q::Float64`: The residual value of the weighted cumulative sum. Default is 0.0.
- `upw::Bool`: Whether to monitor increases in the mean (`true`) or decreases (`false`). Default is `true`.

# References
Shu, L., Jiang, W., & Tsui, K.-L. (2008). A Weighted CUSUM Chart for Detecting Patterned Mean Shifts. Journal of Quality Technology, 40(2), 194-213. https://doi.org/10.1080/00224065.2008.11917725

# Examples
```julia
stat = WCUSUM()
get_design(stat)            # returns [1.0, 0.2]
update_statistic(stat, 3.0) # returns 1.2
stat.Q                      # returns 0.6
```
"""
@with_kw mutable struct WCUSUM <: AbstractStatistic 
    k::Float64 = 1.0
    λ::Float64 = 0.2
    value::Float64 = 0.0
    Q::Float64 = 0.0
    upw::Bool = true
    @assert !isinf(value)
    @assert (k > 0.0) && !isinf(k) "k must be finite and positive"
    @assert 0 < λ <= 1
end
export WCUSUM


get_design(stat::WCUSUM) = [stat.k, stat.λ]

function set_design!(stat::WCUSUM, par::AbstractVector)
    @assert length(par) == 2
    stat.k = first(par)
    stat.λ = last(par)
    return par
end

function update_statistic!(stat::WCUSUM, x::Real)
    # Update wcusum residuals
    stat.Q = (1 - stat.λ) * stat.Q + stat.λ*x
    # Calculate fault signature
    φ = abs(stat.Q)
    # Update statistic
    if stat.upw
        return max(0.0, stat.value + (x - stat.k)*φ)
    else
        return min(0.0, stat.value + (x + stat.k)*φ)
    end
end


update_statistic(stat::WCUSUM, x::Real) = update_statistic!(deepcopy(stat), x)

"""
    AEWMA(λ, k, value)

Adaptive exponentially weighted moving average with design parameters `λ`, `k`, and initial value `value`.

The update mechanism for the statistic ``C_t`` based on a new observation `x` is given by

``C_t = (1-\\phi(e))\\cdot C_{t-1} + \\phi(e) \\cdot x``,

where ``\\phi(e)`` is a forecast error function based on the Huber score function.

### Arguments
- `λ::Float64`: The smoothing constant. Default is `0.1`.
- `k::Float64`: The threshold value in the Huber score. Default is `3.0``.
- `value::Float64`: The initial value of the statistic. Default is 0.0.

### References 
Capizzi, G. & Masarotto, G. (2003). An Adaptive Exponentially Weighted Moving Average Control Chart. Technometrics, 45(3), 199-207.
"""
@with_kw mutable struct AEWMA <: AbstractStatistic 
    λ::Float64 = 0.1
    k::Float64 = 3.0
    value::Float64 = 0.0
    @assert 0.0 < λ <= 1.0
    @assert 0.0 < k
    @assert !isinf(value)
end
export AEWMA

get_design(stat::AEWMA) = [stat.λ, stat.k]

function set_design!(stat::AEWMA, par::AbstractVector)
    @assert length(par) == 2
    stat.λ = first(par)
    stat.k = last(par)
    return par
end
 
function huber(e, λ, k)
    if abs(e) <= k
        return λ*e
    elseif e > k
        return e - (1-λ)*k
    else
        return e + (1-λ)*k
    end
end

update_statistic(stat::AEWMA, x::Real) = stat.value + huber(x - stat.value, stat.λ, stat.k)
update_statistic!(stat::AEWMA, x::Real) = stat.value = update_statistic(stat, x)
