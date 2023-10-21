using DataFrames
using Parameters
using StatsBase
using LinearAlgebra

"""
    LLCUSUM(x::AbstractMatrix, k::Real; ncuts::AbstractVector = [2 for _ in eachcol(x)])

A mutable struct representing the distibution-free CUSUM statistic based on data categorization.

# Fields
- `k::Float64`: The allowance constant of the CUSUM statistic.
- `value::Float64`: The initial value of the statistic (default: 0.0).
- `Sobs::Vector{Float64}`: The vector of cumulative observed categories.
- `Sexp::Vector{Float64}`: The vector of expected observed categories.
- `qtls::Vector{Vector{Float64}}`: A vector of quantiles used for data categorization.
- `f0::Vector{Float64}`: The vector of IC cell probabilities, estimated using a loglinear model.
- `table::Matrix{Float64}`: The matrix of cell combinations.

# References
Qiu, P. (2008). Distribution-free multivariate process control based on log-linear modeling. IIE Transactions, 40(7), 664-677. https://doi.org/10.1080/07408170701744843

"""
@with_kw mutable struct LLCUSUM{M} <: AbstractStatistic
    k::Float64
    value::Float64 = 0.0
    Sobs::Vector{Float64}
    Sexp::Vector{Float64}
    qtls::Vector{Vector{Float64}}
    f0::Vector{Float64} 
    table::M
end
export LLCUSUM

get_design(stat::LLCUSUM) = [stat.k]
set_design!(stat::LLCUSUM, k::AbstractVector) = stat.k = first(k)
set_design!(stat::LLCUSUM, k::Float64) = stat.k = k


function LLCUSUM(k::Float64, x::AbstractMatrix; ncuts::Vector{Int} = [2 for _ in eachcol(x)])
    df_mat, table, qtls = create_table(x, ncuts)
    df = DataFrame(df_mat, :auto)
    rename!(df, Dict(Symbol("x"*string(size(df,2))) => :y))
    f0 = estimate_loglinear_model_probabilities(df, table)
    if !isapprox(sum(f0), 1.0)
        error("Sum of probabilities is different from 1 (value is $(sum(f0)))")
    end
    return LLCUSUM(k=k, value=0.0, Sobs=zeros(length(f0)), Sexp=zeros(length(f0)), qtls=qtls, f0=f0, table=table)
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
function estimate_loglinear_model_probabilities(df, table)
    mod = backward_loglinear(df, :y)
    return predict(mod) / sum(df.y)
end
export estimate_loglinear_model_probabilities

"""
    categorical_to_index(g::AbstractVector, table::AbstractMatrix)

Converts a categorized variable `g` to its corresponding index in a table.

# Arguments
- `g::AbstractVector`: The categorical variable to convert.
- `table::AbstractMatrix`: The table containing all the possible categories.

# Returns
- `Int`: The index of the category in the table.

# Raises
- `ErrorException` if the index is not found

# Example
    x = randn(500, 3)
    df, table, qtls = create_table(x, [2 for _ in eachcol(x)])
    xnew = randn(3)
    g = categorize_data(xnew, qtls)
    idx = categorical_to_index(g, table)
"""
function categorical_to_index(g::AbstractVector, table::AbstractMatrix)
    for i in 1:size(table, 1)
        if g == view(table, i, :)
            return i
        end
    end
    error("Category not found")
end
export categorical_to_index

function update_statistic!(STAT::LLCUSUM, x::AbstractVector)
    ncells = length(STAT.f0)
    g_n = zeros(ncells)
    gObs = categorize_data(x, STAT.qtls)
    idx = categorical_to_index(gObs, STAT.table)
    g_n[idx] = 1.0
    Sobs_nm1 = STAT.Sobs
    Sexp_nm1 = STAT.Sexp
    vv = (Sobs_nm1 - Sexp_nm1) + g_n - STAT.f0
    Lambda = diagm(Sexp_nm1 + STAT.f0)
    C_n = vv' * inv(Lambda) * vv 
    k = STAT.k
    # Equation (4) of [Qiu, 2008]
    if C_n <= k
        STAT.Sobs .= 0.0
        STAT.Sexp .= 0.0
    else
        STAT.Sobs = (STAT.Sobs + g_n) * (C_n - k)/C_n
        STAT.Sexp = (STAT.Sexp + STAT.f0) * (C_n - k)/C_n
    end

    # Add a small random error to make the statistic continuous and
    # attain all values of ARL0
    C_n = C_n + sqrt(0.1)*randn()
    STAT.value = max(0.0, C_n - k)
    return STAT.value
end

update_statistic(STAT::LLCUSUM, x::AbstractVector) = update_statistic!(deepcopy(STAT), x)