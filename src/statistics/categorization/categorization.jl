using DataFrames
using StatsBase
using GLM
using Distributions

function create_table(x::AbstractMatrix, q::AbstractVector)
    # Get quantiles needed to divide each column in d elements
    prob_qtls = [quantile_range(0, 1, q[j]) for j in eachindex(q)]
    @assert length(prob_qtls) == size(x, 2) "Length of quantile vector ($(length(prob_qtls))) must be equal to number of columns in x ($(size(x,2)))"
    # Get vector of quantiles for each column of x
    qtls = [quantile(x[:,j], prob_qtls[j]) for j in 1:size(x,2)]
    df = data_to_contingency_table(x, qtls)
    table = get_all_categories(qtls)
    if size(df, 1) != size(table, 1)
        @warn "Some categories ($(size(table,1)) total) are not present in the data ($(size(df,1)) total). Introducing artificial observations to compensate."
        all_categories = [table[i,:] for i in 1:size(table,1)]
        obs_categories = [df[i,1:(end-1)] for i in 1:size(df,1)]
        not_obs = setdiff(all_categories, obs_categories)
        @show not_obs
        for categ in not_obs
            @show categ
            @show df
            df = vcat(df, [categ; [1.0]]')
        end
        df = sortslices(df, dims=1)
    end
    @assert size(df,1) == size(table,1) 
    @assert view(df, :, 1:(size(df, 2)-1)) == table "Categories are not correctly sorted"
    return df, table, qtls
end
export create_table


"""
    quantile_range(lower, upper, m)

Get the set of quantiles needed to divide the data in m classes.
"""
function quantile_range(lower, upper, m)
    return [lower + (upper-lower) * (j / m) for j in 1:(m-1)]
end
export quantile_range


"""
    categorize_data(x::AbstractVector, qtls::Vector{Vector{Float64}})
    categorize_data(x::AbstractMatrix, qtls::Vector{Vector{Float64}})

Categorize continuous data using the medians of the continuous variables.
"""
function categorize_data(x::AbstractVector, qtls::Vector{Vector{Float64}})
    p = length(x)
    output = zeros(p)
    @assert length(qtls) == length(x) "Length of quantiles ($(length(qtls))) is different from number of variables ($(length(x)))"
    qq = 0.0
    for j in 1:p
        for l in 1:length(qtls[j])
            qq = qtls[j][l]
            output[j] = ifelse(x[j] > qq, l, output[j])
        end
    end
    return output
end

function categorize_data(x::AbstractMatrix, qtls::Vector{Vector{Float64}})
    n,p = size(x)
    output = zeros(n,p)
    for i in 1:n
        output[i, :] = categorize_data(x[i,:], qtls)
    end
    return output
end
export categorize_data

function data_to_contingency_table(x::AbstractMatrix, qtls::Vector{Vector{Float64}})
    dd = length(qtls) + 1
    xCat = categorize_data(x, qtls)
    y = countmap([xCat[i,:] for i in 1:size(xCat,1)])
    k, v = collect(keys(y)), collect(values(y))
    _,p = size(xCat)
    ret = zeros(length(k), p+1)
    for i in 1:size(ret,1)
       ret[i,:] = [k[i]; v[i]] 
    end
    return sortslices(ret, dims=1)
end
export data_to_contingency_table


function get_all_categories(qtls)
    dims = [0:length(qtls[j]) for j in 1:length(qtls)]
    zz = Iterators.product(dims...)
    return zz |> Iterators.flatten |> collect |> x->reshape(x, (length(qtls),:))' |> float |> x->sortslices(x, dims=1)
end
export get_all_categories


"""
    kronecker_matrix(qtls::AbstractVector)

Compute the matrix given by the Kronecker vectors according to Equation (6) of Wang et Al. (2017).
This matrix is used to compute the approximate GLRT statistic to test the null hypothesis that the main effects and the second order interactions of a log-linear model are stable.

# References
Wang, J., Li, J., & Su, Q. (2017). Multivariate Ordinal Categorical Process Control Based on Log-Linear Modeling. Journal of Quality Technology, 49(2), 108-122. https://doi.org/10.1080/00224065.2017.11917983
"""
function kronecker_matrix(qtls::AbstractVector)
    dim = prod(length(q) + 1 for q in qtls)
    tmp = Matrix{Float64}(undef, dim, 0)
    p = length(qtls)
    scoreVecs = [collect(0:(length(qtls[j]))) for j in 1:p]
    onesVecs = [ones(length(qtls[j]) + 1) for j in 1:p]
    kronVecs = Vector{Vector{Float64}}(undef, p)
    for j in 1:p
        for l in 1:p
            if l == j
                kronVecs[l] = scoreVecs[l]
            else
                kronVecs[l] = onesVecs[l]
            end
        end
        tmp = hcat(tmp, kron(kronVecs...))
    end
    for j in 1:(p-1)
        for l in (j+1):p
            for h in 1:p
                kronVecs[h] = onesVecs[h]
            end
            kronVecs[j] = scoreVecs[j]
            kronVecs[l] = scoreVecs[l]
            tmp = hcat(tmp, kron(kronVecs...))
        end
    end
    return Matrix(tmp)
end
export kronecker_matrix


function compose(lhs::Symbol, rhs::AbstractVector{Symbol})
    Formula(lhs, Expr(:call, :+, [1;rhs]...))
end

function stepwise_loglinear(df, lhs::Symbol, rhs::AbstractVector{Symbol},
              forward::Bool, use_aic::Bool)
    options = forward ? setdiff(names(df), [lhs; rhs]) : rhs
    fun = use_aic ? aic : bic
    isempty(options) && return (rhs, false)
    best_fun = fun(lm(compose(lhs, rhs), df))
    improved = false
    best_rhs = rhs
    for opt in options
        this_rhs = forward ? [rhs; opt] : setdiff(rhs, [opt])
        this_fun = fun(lm(compose(lhs, this_rhs), df))
        if this_fun < best_fun
            best_fun = this_fun
            best_rhs = this_rhs
            improved = true
        end
    end
    (best_rhs, improved)
end
export stepwise_loglinear

function stepwise(df, lhs::Symbol, forward::Bool, use_aic::Bool)
    rhs = forward ? Symbol[] : setdiff(names(df), [lhs])
    while true
        rhs, improved = step(df, lhs, rhs, forward, use_aic)
        improved || return lm(compose(lhs, sort(rhs)), df)
    end
end







############## MULTINOMIAL SAMPLER ##############
@with_kw struct MultinomialBootstrap{V,M,VV,D} <: SPM.AbstractSampling
    f0::V
    table::M
    qtls::VV
    dist::D = Multinomial(1, f0)
end
export MultinomialBootstrap

MultinomialBootstrap(stat::SPM.AbstractStatistic) = MultinomialBootstrap(deepcopy(stat.f0), deepcopy(stat.table), deepcopy(stat.qtls), Multinomial(1, deepcopy(stat.f0)))

SPM.new_data(B::MultinomialBootstrap, data::AbstractVector) = error("Not implemented yet.")

function SPM.new_data(B::MultinomialBootstrap, data::AbstractMatrix)
    idx = convert(Vector{Bool}, Distributions.rand(B.dist))
    # @show idx
    newx = B.table[idx, :]
    # @show newx
    jitt = sqrt(eps())
    for j in 1:length(newx)
        if newx[j] == 0.0
            newx[j] = B.qtls[j][1] - jitt
        elseif newx[j] == (float(length(B.qtls[j]) - 1))
            newx[j] = B.qtls[j][end] + jitt
        else
            newx[j] = B.qtls[j][Int(newx[j])] + jitt
        end
    end
    return vec(newx)
end

SPM.new_data!(B::MultinomialBootstrap, data::AbstractVecOrMat) = SPM.new_data(B, data) 
