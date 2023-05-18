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
    return -(float(rl <= get_value(nominal)) - nominal.qtl)
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
    saCL!(CH::ControlChart[; rlsim::Function, settings::OptSettings])

Computes the control limit to satisfy the nominal properties of a control chart, using the stochastic approximation algorithm described in Capizzi and Masarotto (2016).

### Inputs
* `CH` - A control chart.
* `rlsim` - A function that generates new data with signature `rlsim(CH; maxiter, deltaSA)`. If left unspecified, defaults to `run_sim_sa`. See the help for `run_sim_sa` for more information about the signature of the function.
* `settings` - An `OptSettings` objects which contains variables that control the behaviour of the algorithm. See the `Accepted settings` section below for information about the settings that control the behaviour of the algorithm. For more information about the specifics of each keyword argument, see Capizzi and Masarotto (2016).

### Settings
The following settings control the behaviour of the algorithm: 
* `Nfixed_sa` - The number of iterations for the gain estimation stage.
* `Afixed_sa` - The fixed gain during the gain estimation stage.
* `Amin_sa` - The minimum allowed value of gain.
* `Amax_sa` - The maximum allowed value of gain.
* `deltaSA_sa` - The shift in control limit used during the gain estimation stage.
* `q_sa` - The power that controls the denominator in the Robbins-Monro algorithm.
* `gamma_sa` - The precision parameter for the stopping criterion of the algorithm.
* `Nmin_sa` - The minimum number of iterations required for the algorithm to end.
* `z_sa` - The quantile of the `Normal(0,1)` that controls the probability of the stopping criterion being satisfied.
* `Cmrl_sa` - The inflation factor for the maximum number of iterations the run length may run for.
* `maxiter_sa` - Maximum number of iterations before the algorithm is forcibly ended.
* `verbose_sa` - Whether to print information to the user about the state of the optimization.

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Capizzi, G., & Masarotto, G. (2016). "Efficient Control Chart Calibration by Simulated Stochastic Approximation". IIE Transactions 48 (1). https://doi.org/10.1080/0740817X.2015.1055392.
"""
function saCL!(CH::ControlChart; settings::OptSettings = OptSettings())
    @unpack rlsim, hmin_sa, Nfixed_sa, Afixed_sa, Amin_sa, Amax_sa, delta_sa, q_sa, gamma_sa, Nmin_sa, z_sa, Cmrl_sa, maxiter_sa, verbose_sa = settings

    tmp = rlsim(CH; maxiter=1, deltaSA=0.0)
    @assert haskey(tmp, :rl) "rlsim function must have key :rl"
    @assert haskey(tmp, :rlPlus) "rlsim function must have key :rlPlus"
    @assert haskey(tmp, :rlMinus) "rlsim function must have key :rlMinus"

    @assert Nfixed_sa > 0 "Nfixed_sa must be positive"
    @assert Amin_sa > 0 "Amin_sa must be positive"
    @assert Amax_sa > 0 "Amax_sa must be positive"
    @assert delta_sa > 0 "deltaSA_sa must be positive"
    @assert 0 < q_sa < 1 "q_sa must be a number between 0 and 1"
    @assert 0 < gamma_sa < 1 "gamma_sa must be a number between 0 and 1"
    @assert Nmin_sa > 0 "Nmin_sa must be positive"
    @assert z_sa > 0 "z_sa must be positive"
    @assert Cmrl_sa > 0 "Cmrl_sa must be positive"
    @assert maxiter_sa > 0 "maxiter_sa must be positive"

    v = (z_sa/gamma_sa)^2
    h = deepcopy(get_h(get_limit(CH)))
    score = zero(h)
    s2    = zero(h)
    D     = zero(h)

    if verbose_sa println("Running SA ...") end

    # Stage 1 - adaptive gain estimation
    D          = zero(h)
    rl         = zero(h)
    rlPlus     = zero(h)
    rlMinus    = zero(h)
    scorePlus  = zero(h)
    scoreMinus = zero(h)
    Arl0 = get_nominal_value(CH)
    if verbose_sa println("Running adaptive gain ...") end
    for i in 1:Nfixed_sa
        set_limit!(CH, h)
        rl, rlPlus, rlMinus = rlsim(CH, maxiter=Cmrl_sa * Arl0 * sqrt(i + Nfixed_sa), deltaSA=delta_sa)
        # @show rl, rlPlus, rlMinus, h
        score = calculate_limit_gradient(CH, rl)
        h = calculate_limit(h, Afixed_sa, score, i, q_sa, hmin_sa)
        scorePlus = calculate_limit_gradient(CH, rlPlus)
        scoreMinus = calculate_limit_gradient(CH, rlMinus)
        # @show scorePlus, scoreMinus, D
        # println("")
        D = update_gain(D, scorePlus, scoreMinus, i)
    end

    D = calculate_gain(D, Amin_sa, Amax_sa, delta_sa)

    if verbose_sa println("Estimated gain D = $(D)") end


    # Stage 2 - Stochastic approximations
    hm = zero(h)
    i = 0
    conv = "Maximum number of iterations reached"
    if verbose_sa println("Running optimization ...") end
    while i < (maxiter_sa + Nmin_sa)
        if verbose_sa && (i % floor(maxiter_sa / 50) == 0)
            println("i: $(i)/$(Int(trunc(maxiter_sa + Nmin_sa)))\th: $(round.(h, digits=5))\thm: $(round.(hm, digits=5))")
        end
        i += 1
        set_limit!(CH, h)
        rl, _, _ = rlsim(CH, maxiter=Cmrl_sa * Arl0 * sqrt(i + Nfixed_sa), deltaSA=0.0)
        score = calculate_limit_gradient(CH, rl)
        h = calculate_limit(h, D, score, i, q_sa, hmin_sa)
        if i > Nmin_sa
            Ndenom = (i - Nmin_sa)
            hm = hm + (h - hm) / Ndenom 
            s2 = update_score(s2, score, Ndenom)
        end
        if (i > Nmin_sa) && (i > v * maximum(s2))
            if verbose_sa println("i: $(i)/$(Int(trunc(maxiter_sa)))\tConvergence!\n") end
            conv = "Convergence"
            break
        end
    end
    set_limit!(CH, hm)
    return (h=hm, iter=i, status = conv)
end
export saCL!

"""
    saCL(CH::ControlChart[; rlsim::Function, settings::OptSettings])

Applies the stochastic approximation algorithm of Capizzi and Masarotto (2016) without modifying the control chart object `CH`.
See the documentation of `saCL!` for more information about the algorithm and the keyword arguments.

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Capizzi, G., & Masarotto, G. (2016). "Efficient Control Chart Calibration by Simulated Stochastic Approximation". IIE Transactions 48 (1). https://doi.org/10.1080/0740817X.2015.1055392.
"""
function saCL(CH::ControlChart; settings::OptSettings = OptSettings())
    CH_ = shallow_copy_sim(CH)
    return saCL!(CH_, settings=settings)
end
export saCL