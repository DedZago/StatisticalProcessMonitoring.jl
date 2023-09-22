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

"""
    compose(lhs::Symbol, rhs::Tuple)
    compose(lhs::Symbol, rhs::Term)
    compose(lhs::Symbol, rhs::ConstantTerm)

Compose a response variable `lhs` with a tuple or vector of predictors `rhs` terms.

# Arguments
- `lhs::Symbol`: The left-hand side of the equation.
- `rhs::Tuple`: The right-hand side of the equation.

# Returns
- The composed equation.

If `rhs` is empty, the function returns `lhs` composed with an intercept.
Otherwise, the function returns `lhs` composed with the terms in `rhs`.
"""
function compose(lhs::Symbol, rhs::Tuple)
    if length(rhs) == 0
        return term(lhs) ~ ConstantTerm(1)
    else
        return term(lhs) ~ rhs
    end
end

function compose(lhs::Symbol, rhs::Term)
    return term(lhs) ~ rhs
end

function compose(lhs::Symbol, rhs::ConstantTerm)
    return term(lhs) ~ rhs
end
export compose


function get_term_vector(term::InteractionTerm)
    return term.terms
end

function get_term_vector(term::Term)
    return [term]
end

"""
    get_removable_terms(rhs)

Given an iterable collection of terms `rhs`, this function returns the indices of the terms in `rhs` that can be removed in a backward elimination step.

# Arguments
- `rhs`: An iterable collection of terms in the equation.

# Returns
- An array of indices of the terms in `rhs` that have the maximum order.
"""
function get_removable_terms(rhs)
    terms_order = collect(length.(terms.(rhs)))     # Get the order of the terms in `rhs`
    max_order = maximum(terms_order)                # Maximum order can be removed
    removable_terms = (1:length(rhs))[terms_order .== max_order]        # Get maximum order elements
    max_order == 1 && return removable_terms
    # Get terms that do not appear in higher-order terms
    for j in reverse(1:(max_order-1))
        for k in (1:length(rhs))[terms_order .== j]
            candidate_removable = rhs[k]
            remv = true
            for h in (j+1):max_order
                for tt in rhs[terms_order .== h]
                    if issubset(Set(get_term_vector(candidate_removable)), Set(tt.terms))
                        remv = false
                    end
                end
            end
            if remv
                push!(removable_terms, k)
            end
        end
    end
    return removable_terms
end
export get_removable_terms


function step_backward(df, lhs::Symbol, rhs, use_aic::Bool)
    options = rhs
    fun = use_aic ? aic : bic
    isa(rhs, ConstantTerm) && return(rhs, false)
    isempty(options) && return (rhs, false)
    best_fun = fun(glm(compose(lhs, rhs), df, Poisson()))
    improved = false
    best_rhs = rhs
    @show best_rhs
    removable = get_removable_terms(rhs)
    for i in removable
        opt = options[i]
        this_rhs = foldl(+, setdiff(rhs, [opt]))
        formula_current = compose(lhs, this_rhs)
        @show this_rhs
        this_fun = fun(glm(formula_current, df, Poisson()))
        if this_fun < best_fun
            best_fun = this_fun
            best_rhs = this_rhs
            improved = true
        end
    end
    (best_rhs, improved)
end

function step_backward(df, lhs::Symbol, rhs::Term, use_aic::Bool)
    fun = use_aic ? aic : bic
    best_fun = fun(glm(compose(lhs, rhs), df, Poisson()))
    best_rhs = rhs
    improved = false
    this_rhs = ConstantTerm(1)
    this_fun = fun(glm(term(lhs) ~ this_rhs, df, Poisson()))
    if this_fun < best_fun
        best_rhs = this_rhs
        improved = true
    end
    (best_rhs, improved)
end

function backward_loglinear(df, lhs::Symbol; use_aic::Bool = false)
    rhs_symbol = setdiff(Symbol.(names(df)), [lhs])
    rhs = foldl(*, term.(rhs_symbol))
    @show rhs
    while true
        rhs, improved = step_backward(df, lhs, rhs, use_aic)
        improved || return glm(compose(lhs, rhs), df, Poisson())
    end
end
export backward_loglinear







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
