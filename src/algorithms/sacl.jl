using Statistics

"""
    calculate_limit_gradient(CH::AbstractChart, rl::Real)
    calculate_limit_gradient(nominal::ARL, rl)
    calculate_limit_gradient(nominal::QRL, rl)

Calculate the gradient for the optimization of the control limit.

If the control chart `nominal` attribute is of type `ARL`, then the gradient is calculated according to
Equation (9) of Capizzi and Masarotto (2016).

If the control chart `nominal` attribute is of type `QRL`, then the gradient is calculated using the recursion on page 280 of Capizzi and Masarotto (2009)

### References
Capizzi, G., & Masarotto, G. (2016). "Efficient Control Chart Calibration by Simulated Stochastic Approximation". IIE Transactions 48 (1). https://doi.org/10.1080/0740817X.2015.1055392.

Capizzi, G. & Masarotto, G. (2009) Bootstrap-based design of residual control charts, IIE Transactions, 41:4, 275-286, DOI: https://doi.org/10.1080/07408170802120059 

"""
function calculate_limit_gradient(CH::AbstractChart, rl)
   return calculate_limit_gradient(get_nominal(CH), rl)
end
export calculate_limit_gradient

function calculate_limit_gradient(nominal::ARL, rl)
    #FIXME: test
    return (rl - get_value(nominal)) / get_value(nominal)
end

function calculate_limit_gradient(nominal::ARL, rl::AbstractVector)
    #FIXME: test
    return (minimum(rl) - get_value(nominal) .+ rl .- Statistics.mean(rl)) ./ get_value(nominal)
end

function calculate_limit_gradient(nominal::QRL, rl)
    #FIXME: test
    return -(float(minimum(rl) <= get_value(nominal)) - nominal.qtl)
end

function calculate_limit_gradient(nominal::QRL, rl::AbstractVector)
    #? Penalize the quantile so each chart has the same median run length
    return -(float(minimum(rl) <= get_value(nominal)) - nominal.qtl) .+ (rl .- mean(rl)) ./ get_value(nominal)
end

calculate_gain(D::Float64, Amin, Amax, deltaSA) = 1.0 / (max(1.0/Amax, min(1.0/Amin, D/(2.0*deltaSA))))
calculate_gain(D::Vector{Float64}, Amin, Amax, deltaSA) = @. 1.0 / (max(1.0/Amax, min(1.0/Amin, D/(2.0*deltaSA))))

update_gain(D::Float64, scorePlus, scoreMinus, i) = D + (scorePlus - scoreMinus) / i
update_gain(D::Vector{Float64}, scorePlus, scoreMinus, i) = D .+ (scorePlus - scoreMinus) / i

calculate_limit(h::Float64, D, score, i, q, eps) = max(eps, h - D * score / (i^q))
calculate_limit(h::Vector{Float64}, D, score, i, q, eps) = @. max(eps, h - D * score / (i^q))

update_score(s2::Float64, score::Real, Ndenom) = s2 + (score * score - s2) / Ndenom
update_score(s2::Vector{Float64}, score::AbstractVector, Ndenom) = s2 .+ (score .* score .- s2) ./ Ndenom

"""
    saCL!(CH::ControlChart; kw...)

Computes the control limit to satisfy the nominal properties of a control chart, using the stochastic approximation algorithm described in Capizzi and Masarotto (2016).

### Inputs
* `CH` - A control chart.
* `kw...` - Keyword arguments that control the behaviour of the algorithm. For more information about the specifics of each keyword argument, see Capizzi and Masarotto (2016).

### Keyword arguments:
* `rlsim` - A function that generates new data with signature `rlsim(CH, maxiter, deltaSA)`. If left unspecified, defaults to `new_data_sa`. See the help for `new_data_sa` for more information.
* `Nfixed` - The number of iterations for the gain estimation stage.
* `Afixed` - The fixed gain during the gain estimation stage.
* `Amin` - The minimum allowed value of gain.
* `Amax` - The maximum allowed value of gain.
* `deltaSA` - The shift in control limit used during the gain estimation stage.
* `q` - The power that controls the denominator in the Robbins-Monro algorithm.
* `gamma` - The precision parameter for the stopping criterion of the algorithm.
* `Nmin` - The minimum number of iterations required for the algorithm to end.
* `z` - The quantile of the `Normal(0,1)` that controls the probability of the stopping criterion being satisfied.
* `Cmrl` - The inflation factor for the maximum number of iterations the run length may run for.
* `maxiter` - Maximum number of iterations before the algorithm is forcibly ended.
* `verbose` - Whether to print information to the user about the state of the optimization.

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Capizzi, G., & Masarotto, G. (2016). "Efficient Control Chart Calibration by Simulated Stochastic Approximation". IIE Transactions 48 (1). https://doi.org/10.1080/0740817X.2015.1055392.
"""
function saCL!(CH::ControlChart; rlsim::Function = run_sim_sa, Nfixed::Int=500, Afixed::Real=0.1, Amin::Real=0.1, Amax::Real=100.0, deltaSA::Real=0.1, q::Real=0.55, gamma::Real=0.02, Nmin::Int=1000, z::Real = 3.0, Cmrl::Real=10.0, maxiter::Real = 1e05, verbose::Bool=true)

    #TODO: add checks on the rlsim function

    @assert Nfixed > 0 "Nfixed must be positive"
    @assert Amin > 0 "Amin must be positive"
    @assert Amax > 0 "Amax must be positive"
    @assert deltaSA > 0 "deltaSA must be positive"
    @assert 0 < q < 1 "q must be a number between 0 and 1"
    @assert 0 < gamma < 1 "gamma must be a number between 0 and 1"
    @assert Nmin > 0 "Nmin must be positive"
    @assert z > 0 "z must be positive"
    @assert Cmrl > 0 "Cmrl must be positive"
    @assert maxiter > 0 "maxiter must be positive"

    v = (z/gamma)^2
    eps = 1e-06
    h = deepcopy(get_value(get_limit(CH)))
    score = zero(h)
    s2    = zero(h)
    D     = zero(h)

    if verbose println("Running SA ...") end

    # Stage 1 - adaptive gain estimation
    D          = zero(h)
    rl         = zero(h)
    rlPlus     = zero(h)
    rlMinus    = zero(h)
    scorePlus  = zero(h)
    scoreMinus = zero(h)
    Arl0 = get_nominal_value(CH)
    if verbose println("Running adaptive gain ...") end
    for i in 1:Nfixed
        set_limit!(CH, h)
        rl, rlPlus, rlMinus = rlsim(CH, Cmrl * Arl0 * sqrt(i + Nfixed), deltaSA)
        # @show rl, rlPlus, rlMinus, h
        score = calculate_limit_gradient(CH, rl)
        h = calculate_limit(h, Afixed, score, i, q, eps)
        scorePlus = calculate_limit_gradient(CH, rlPlus)
        scoreMinus = calculate_limit_gradient(CH, rlMinus)
        # @show scorePlus, scoreMinus, D
        # println("")
        D = update_gain(D, scorePlus, scoreMinus, i)
    end

    D = calculate_gain(D, Amin, Amax, deltaSA)

    if verbose println("Estimated gain D = $(D)") end


    # Stage 2 - Stochastic approximations
    hm = zero(h)
    i = 0
    conv = "Maximum number of iterations reached"
    if verbose println("Running optimization ...") end
    while i < (maxiter + Nmin)
        if verbose && (i % floor(maxiter / 50) == 0)
            println("i: $(i)/$(Int(trunc(maxiter)))\th: $(round.(h, digits=5))\thm: $(round.(hm, digits=5))")
        end
        i += 1
        set_limit!(CH, h)
        rl, _, _ = rlsim(CH, Cmrl * Arl0 * sqrt(i + Nfixed), 0.0)
        score = calculate_limit_gradient(CH, rl)
        h = calculate_limit(h, D, score, i, q, eps)
        if i > Nmin
            Ndenom = (i - Nmin)
            hm = hm + (h - hm) / Ndenom 
            s2 = update_score(s2, score, Ndenom)
        end
        if (i > Nmin) && (i > v * maximum(s2))
            if verbose println("i: $(i)/$(Int(trunc(maxiter)))\tConvergence!\n") end
            conv = "Convergence"
            break
        end
    end
    set_limit!(CH, hm)
    return (h=hm, iter=i, status = conv)
end
export saCL!

"""
    saCL(CH::ControlChart; kw...)

Applies the stochastic approximation algorithm of Capizzi and Masarotto (2016) without modifying the control chart object `CH`.
See the documentation of `saCL!` for more information about the algorithm and keyword arguments.

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Capizzi, G., & Masarotto, G. (2016). "Efficient Control Chart Calibration by Simulated Stochastic Approximation". IIE Transactions 48 (1). https://doi.org/10.1080/0740817X.2015.1055392.
"""
function saCL(CH::ControlChart; rlsim::Function = run_sim_sa, Nfixed::Int=500, Afixed::Real=0.1, Amin::Real=0.1, Amax::Real=100.0, deltaSA::Real=0.1, q::Real=0.55, gamma::Real=0.02, Nmin::Int=1000, z::Real = 3.0, Cmrl::Real=10.0, maxiter::Real = 4e05, verbose::Bool=true)
    CH_ = shallow_copy_sim(CH)
    return saCL!(CH_, rlsim = rlsim, Nfixed = Nfixed, Afixed = Afixed, Amin = Amin, Amax = Amax, deltaSA = deltaSA, q = q, gamma = gamma, Nmin = Nmin, z = z, Cmrl = Cmrl, maxiter = maxiter, verbose = verbose)
end
export saCL