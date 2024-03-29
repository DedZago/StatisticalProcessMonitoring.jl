using DataFrames
using Parameters
using LinearAlgebra

"""
    LLD(x::AbstractMatrix, l::Real; ncuts::AbstractVector = [3 for _ in eachcol(x)], N = 1)

A Log-Linear Directional monitoring statistic based on data categorization.

# Fields
- `l::Float64`: The exponentially weighted smoothing constant of the statistic.
- `value::Float64`: The initial value of the statistic (default: 0.0).
- `qtls::Vector{Vector{Float64}}`: A vector of quantiles used for data categorization.
- `f0::Vector{Float64}`: The vector of IC cell probabilities, estimated using a loglinear model.
- `table::Matrix{Float64}`: The matrix of cell combinations.
- `Sigma::Matrix{Float64}`: The covariance matrix of the cell counts.
- `directions::Matrix{Float64}`: The matrix of directions to be used in the hypohesis test (see Equation (6) of [Li, 2012]). This is calculated according to the procedure described in Equations (5) and (6) of [Wang, 2017].
- `N::Int`: The number of observations at each time point (default: 1)
- `z_k::Vector{Float64}`: The current vector of smoothed values.

# References
Li, J., Tsung, F., & Zou, C. (2012). Directional Control Schemes for Multivariate Categorical Processes. Journal of Quality Technology, 44(2), 136â€“154. https://doi.org/10.1080/00224065.2012.11917889

Wang, J., Li, J., & Su, Q. (2017). Multivariate Ordinal Categorical Process Control Based on Log-Linear Modeling. Journal of Quality Technology, 49(2), 108-122. https://doi.org/10.1080/00224065.2017.11917983
"""
@with_kw mutable struct LLD <: AbstractStatistic
    l::Float64
    value::Float64 = 0.0
    qtls::Vector{Vector{Float64}}
    f0::Vector{Float64} 
    table::Matrix{Float64}
    Sigma::Matrix{Float64} = diagm(f0) - f0*f0'
    directions::Matrix{Float64} = kronecker_matrix(qtls)
    N::Int = 1
    z_k::Vector{Float64} = deepcopy(f0)
end
export LLD

get_design(stat::LLD) = [stat.l]
set_design!(stat::LLD, l::AbstractVector) = stat.l = first(l)
set_design!(stat::LLD, l::Float64) = stat.l = l


function LLD(l::Real, x::AbstractMatrix; ncuts::AbstractVector = [3 for _ in eachcol(x)], N = 1)
    @assert length(ncuts) == size(x,2) "Must provide a number of classes for each variable ($(size(x,2)) total, $(length(ncuts)) provided)"
    df_mat, table, qtls = create_table(x, ncuts)
    df = DataFrame(df_mat, :auto)
    rename!(df, Dict(Symbol("x"*string(size(df,2))) => :y))
    f0 = estimate_ordinal_model_probabilities(df, table)
    @assert isapprox(sum(f0), 1.0) "Sum of probabilities is different from 1 (value is $(sum(f0)))"
    #----- Create GLRT vcov matrix -----#
    return LLD(l=l, value=0.0, qtls=qtls, f0=f0, table=table, N = N)
end


function update_statistic!(STAT::LLD, x)
    ncells = length(STAT.f0)
    g_n = zeros(ncells)
    gObs = categorize_data(x, STAT.qtls)
    idx = categorical_to_index(gObs, STAT.table)
    g_n[idx] = 1.0
    STAT.z_k = (1.0 - STAT.l) * STAT.z_k + STAT.l * g_n
    glrt = -Inf
    quadform = 0.0
    for j in 1:length(STAT.qtls)
        x = view(STAT.directions, :, j)
        quadform = 1/(STAT.N) * dot(STAT.z_k - STAT.N*STAT.f0, x*inv(x'STAT.Sigma*x)*x',STAT.z_k - STAT.N*STAT.f0)
        if quadform > glrt
            glrt = quadform
        end
    end
    STAT.value = glrt
    return STAT.value
end
