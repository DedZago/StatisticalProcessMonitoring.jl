using DataFrames
using Parameters
using LinearAlgebra
import LinearAlgebra: kron

"""
    MOC(x::AbstractMatrix, l::Real; ncuts::AbstractVector = [3 for _ in eachcol(x)], N = 1)

A Multivariate Ordinal Categorical monitoring statistic based on data categorization.

# Fields
- `l::Float64`: The exponentially weighted smoothing constant of the statistic.
- `value::Float64`: The initial value of the statistic (default: 0.0).
- `qtls::Vector{Vector{Float64}}`: A vector of quantiles used for data categorization.
- `f0::Vector{Float64}`: The vector of IC cell probabilities, estimated using a loglinear model.
- `table::Matrix{Float64}`: The matrix of cell combinations.
- `inv_VCOV::Matrix{Float64}`: The matrix containing the inverse covariance to be used in the running statistic.
- `N::Int`: The number of observations at each time point (default: 1)
- `z_k::Vector{Float64}`: The current vector of smoothed values.

# References
Wang, J., Li, J., & Su, Q. (2017). Multivariate Ordinal Categorical Process Control Based on Log-Linear Modeling. Journal of Quality Technology, 49(2), 108-122. https://doi.org/10.1080/00224065.2017.11917983

"""
@with_kw mutable struct MOC <: AbstractStatistic
    l::Float64
    value::Float64 = 0.0
    qtls::Vector{Vector{Float64}}
    f0::Vector{Float64} 
    table::Matrix{Float64}
    inv_VCOV::Matrix{Float64}
    N::Int = 1
    z_k::Vector{Float64} = deepcopy(f0)
end
export MOC

get_design(stat::MOC) = [stat.l]
set_design!(stat::MOC, l::AbstractVector) = stat.l = first(l)
set_design!(stat::MOC, l::Float64) = stat.l = l

function MOC(l::Real, x::AbstractMatrix; ncuts::AbstractVector = [3 for _ in eachcol(x)], N = 1)
    @assert length(ncuts) == size(x,2) "Must provide a number of classes for each variable ($(size(x,2)) total, $(length(ncuts)) provided)"
    df_mat, table, qtls = create_table(x, ncuts)
    df = DataFrame(df_mat, :auto)
    rename!(df, Dict(Symbol("x"*string(size(df,2))) => :y))
    f0 = estimate_ordinal_model_probabilities(df, table)
    @assert isapprox(sum(f0), 1.0) "Sum of probabilities is different from 1 (value is $(sum(f0)))"
    #----- Create GLRT vcov matrix -----#
    Y::Matrix{Float64} = kronecker_matrix(qtls)
    Lambda::Matrix{Float64} = diagm(f0) - f0*f0'
    return MOC(l=l, value=0.0, qtls=qtls, f0=f0, table=table, N = N, inv_VCOV = Y*inv(Y'*Lambda*Y)*Y')
end

"""
    estimate_ordinal_model_probabilities(df, table)

Estimates the probabilities of an ordinal loglinear model based on the observed cell counts.

# Arguments
- `df::DataFrame`: The data frame containing the predictor variables.
- `table::AbstractMatrix`: The matrix of predictor values for which probabilities are to be estimated.

# Returns
- `prob::Vector{Float64}`: A vector of estimated probabilities.
"""
function estimate_ordinal_model_probabilities(df::DataFrame, table; order=2)
    rhs_symbol = setdiff(Symbol.(names(df)), [:y])
    rhs = foldl(*, term.(rhs_symbol))
    terms_order = collect(length.(terms.(rhs)))     # Get the order of the terms in `rhs`
    rhs = rhs[terms_order .<= order]
    mod = glm(term(:y) ~ rhs, df, Poisson())
    return predict(mod) / sum(df.y)
end
export estimate_ordinal_model_probabilities


function update_statistic!(STAT::MOC, x)
    ncells = length(STAT.f0)
    g_n = zeros(ncells)
    gObs = categorize_data(x, STAT.qtls)
    idx = categorical_to_index(gObs, STAT.table)
    g_n[idx] = 1.0
    # @show gObs, g_n, STAT.f0
    
    STAT.z_k = (1.0 - STAT.l) * STAT.z_k + STAT.l * g_n
    # Equation (7) of [Wang2017]
    STAT.value = 1.0/(STAT.N) * (STAT.z_k - STAT.N*STAT.f0)' * STAT.inv_VCOV * (STAT.z_k - STAT.N*STAT.f0)
    return STAT.value
end
