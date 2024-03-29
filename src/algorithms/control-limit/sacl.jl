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
    return (rl - get_value(nominal)) / get_value(nominal)
end

function calculate_limit_gradient(nominal::ARL, rl::AbstractVector)
    return (minimum(rl) - get_value(nominal) .+ rl .- Statistics.mean(rl)) ./ get_value(nominal)
end

function calculate_limit_gradient(nominal::QRL, rl)
    return -(float(rl <= get_value(nominal)) - nominal.qtl)
end

function calculate_limit_gradient(nominal::QRL, rl::AbstractVector)
    #? Penalize the quantile so each chart has the same median run length
    return -(float(minimum(rl) <= get_value(nominal)) - nominal.qtl) .+ (rl .- mean(rl)) ./ get_value(nominal)
end

calculate_gain(D::Float64, Amin, Amax, delta) = 1.0 / (max(1.0/Amax, min(1.0/Amin, D/(2.0*delta))))
calculate_gain(D::Vector{Float64}, Amin, Amax, delta) = @. 1.0 / (max(1.0/Amax, min(1.0/Amin, D/(2.0*delta))))

update_gain(D::Float64, scorePlus, scoreMinus, i) = D + (scorePlus - scoreMinus) / i
update_gain(D::Vector{Float64}, scorePlus, scoreMinus, i) = D .+ (scorePlus - scoreMinus) / i

calculate_limit(h::Float64, D, score, i, q, eps) = max(eps, h - D * score / (i^q))
calculate_limit(h::Vector{Float64}, D, score, i, q, eps) = @. max(eps, h - D * score / (i^q))

update_score(s2::Float64, score::Real, Ndenom) = s2 + (score * score - s2) / Ndenom
update_score(s2::Vector{Float64}, score::AbstractVector, Ndenom) = s2 .+ (score .* score .- s2) ./ Ndenom

"""
    saCL!(CH::ControlChart[; rlsim::Function, kw...])

Computes the control limit to satisfy the nominal properties of a control chart, using the stochastic approximation algorithm described in Capizzi and Masarotto (2016).

### Arguments
* `CH` - A control chart.

### Keyword arguments
* `rlsim` - A function that generates new data with signature `rlsim(CH; maxiter, delta)`. If left unspecified, defaults to `run_sim_sa`. See the help for `run_sim_sa` for more information about the requirements of the function.
* `settings` - An `OptSettings` objects which contains variables that control the behaviour of the algorithm. See the `Accepted settings` section below for information about the settings that control the behaviour of the algorithm. For more information about the specifics of each keyword argument, see Capizzi and Masarotto (2016).
* `Nfixed` - The number of iterations for the gain estimation stage (default: 500).
* `Afixed` - The fixed gain during the gain estimation stage (default: 0.1).
* `Amin` - The minimum allowed value of gain (default: 0.1).
* `Amax` - The maximum allowed value of gain (default: 100.0).
* `delta_sa` - The shift in control limit used during the gain estimation stage (default: 0.1).
* `q` - The power that controls the denominator in the Robbins-Monro algorithm (default: 0.55).
* `gamma` - The precision parameter for the stopping criterion of the algorithm (default: 0.02).
* `Nmin` - The minimum number of iterations to avoid early terminations (default: 1000).
* `z` - The quantile of the `Normal(0,1)` that controls the probability of the stopping criterion being satisfied (default: 3.0).
* `Cmrl` - The inflation factor for the maximum number of iterations the run length may run for (default: 10.0).
* `maxiter` - Maximum number of iterations before the algorithm is forcibly ended (default: 50000).
* `verbose` - Whether to print information to the user about the state of the optimization (default: false).
* `parallel::Bool` - Whether the algorithm should be run in parallel, using available threads (default: false). Parallelization is achieved by averaging `Threads.nthreads` independent replications of the algorithm, each with precision parameter `gamma*sqrt(Threads.nthreads)`. See [Capizzi, 2016] for further discussion on parallelizing the SA algorithm.

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Capizzi, G., & Masarotto, G. (2016). "Efficient Control Chart Calibration by Simulated Stochastic Approximation". IIE Transactions 48 (1). https://doi.org/10.1080/0740817X.2015.1055392.
"""
function saCL!(CH::ControlChart; rlsim::Function = run_sim_sa, hmin::Float64 = sqrt(eps()), Nfixed::Int = 500, Afixed::Float64 = 0.1, Amin::Float64 = 0.1, Amax::Float64 = 100.0, delta_sa::Float64 = 0.1, q::Float64 = 0.55, gamma::Float64 = 0.02, Nmin::Int = 1000, z::Float64 = 3.0, Cmrl::Float64 = 10.0, maxiter::Int = 50_000, verbose::Bool = false, parallel::Bool = false)
    if parallel
        nthr = Threads.nthreads()
        output_list = [[deepcopy(get_h(get_limit(CH))), 0.0, "Convergence"] for _ in 1:nthr]
        Threads.@threads for i in 1:nthr
            output_list[i] = collect(saCL_internal!(CH, rlsim = rlsim, hmin = hmin, Nfixed = Nfixed, Afixed = Afixed, Amin = Amin, Amax = Amax, delta_sa = delta_sa, q = q, gamma = gamma*sqrt(nthr), Nmin = Nmin, z = z, Cmrl = Cmrl, maxiter = maxiter, verbose = verbose))
        end
        status = ifelse(all([x[3] for x in output_list] .== "Convergence"), "Convergence", "Maximum number of iterations reached")
        return (h = mean(x[1] for x in output_list), iter = Int(trunc(mean(x[2] for x in output_list))), status = status)
    else
        return saCL_internal!(CH, rlsim = rlsim, hmin = hmin, Nfixed = Nfixed, Afixed = Afixed, Amin = Amin, Amax = Amax, delta_sa = delta_sa, q = q, gamma = gamma, Nmin = Nmin, z = z, Cmrl = Cmrl, maxiter = maxiter, verbose = verbose)
    end
end
export saCL


function saCL_internal!(CH::ControlChart; rlsim::Function = run_sim_sa, hmin::Float64 = sqrt(eps()), Nfixed::Int = 500, Afixed::Float64 = 0.1, Amin::Float64 = 0.1, Amax::Float64 = 100.0, delta_sa::Float64 = 0.1, q::Float64 = 0.55, gamma::Float64 = 0.02, Nmin::Int = 1000, z::Float64 = 3.0, Cmrl::Float64 = 10.0, maxiter::Int = 50_000, verbose::Bool = false)

    tmp = rlsim(CH; maxiter=1, delta=0.0)
    @assert haskey(tmp, :rl) "rlsim function must have key :rl"
    @assert haskey(tmp, :rlPlus) "rlsim function must have key :rlPlus"
    @assert haskey(tmp, :rlMinus) "rlsim function must have key :rlMinus"

    @assert Nfixed > 0 "Nfixed must be positive"
    @assert Amin > 0 "Amin must be positive"
    @assert Amax > 0 "Amax must be positive"
    @assert delta_sa > 0 "delta_sa must be positive"
    @assert 0 < q < 1 "q must be a number between 0 and 1"
    @assert 0 < gamma < 1 "gamma must be a number between 0 and 1"
    @assert Nmin > 0 "Nmin must be positive"
    @assert z > 0 "z must be positive"
    @assert Cmrl > 0 "Cmrl must be positive"
    @assert maxiter > 0 "maxiter must be positive"

    v = (z/gamma)^2
    h = deepcopy(get_h(get_limit(CH)))
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
        rl, rlPlus, rlMinus = rlsim(CH, maxiter=Cmrl * Arl0 * sqrt(i + Nfixed), delta=delta_sa)
        # @show rl, rlPlus, rlMinus, h
        score = calculate_limit_gradient(CH, rl)
        h = calculate_limit(h, Afixed, score, i, q, hmin)
        scorePlus = calculate_limit_gradient(CH, rlPlus)
        scoreMinus = calculate_limit_gradient(CH, rlMinus)
        # @show scorePlus, scoreMinus, D
        # println("")
        D = update_gain(D, scorePlus, scoreMinus, i)
    end

    D = calculate_gain(D, Amin, Amax, delta_sa)

    if verbose println("Estimated gain D = $(D)") end


    # Stage 2 - Stochastic approximations
    hm = zero(h)
    i = 0
    criter = 0.0
    conv = "Maximum number of iterations reached"
    if verbose println("Running optimization ...") end
    while i < maxiter
        if verbose && (i % floor(maxiter / 50) == 0)
            println("i: $(i)/$(Int(trunc(maxiter)))\th: $(round.(h, digits=5))\thm: $(round.(hm, digits=5))\tstop: $(Int(trunc(criter)))")
        end
        i += 1
        set_limit!(CH, h)
        rl, _, _ = rlsim(CH, maxiter=Cmrl * Arl0 * sqrt(i + Nfixed), delta=0.0)
        score = calculate_limit_gradient(CH, rl)
        h = calculate_limit(h, D, score, i, q, hmin)
        hm = hm + (h - hm) / i 
        s2 = update_score(s2, score, i)
        criter = v * maximum(s2)
        if (i > Nmin) && (i > criter)
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
    saCL(CH::ControlChart[; rlsim::Function, kw...])

Applies the stochastic approximation algorithm of Capizzi and Masarotto (2016) without modifying the control chart object `CH`.

### Keyword arguments
See the documentation of `saCL!` for more information about the algorithm and the keyword arguments.

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Capizzi, G., & Masarotto, G. (2016). "Efficient Control Chart Calibration by Simulated Stochastic Approximation". IIE Transactions 48 (1). https://doi.org/10.1080/0740817X.2015.1055392.
"""
function saCL(CH::ControlChart; kw...)
    CH_ = shallow_copy_sim(CH)
    return saCL!(CH_; kw...)
end
export saCL