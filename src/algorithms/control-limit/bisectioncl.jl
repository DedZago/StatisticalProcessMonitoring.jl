using Statistics

"""
    bisectionCL!(CH::ControlChart[; rlsim::Function, settings::OptSettings])

Computes the control limit to satisfy the nominal properties of a control chart, using the bisection algorithm (see for instance Qiu, 2013)

### Inputs
* `CH` - A control chart.
* `hmax` - The maximum value for the control limit.

### Accepted settings:
* `rlsim` - A function that generates a run length for the control chart with signature `rlsim(CH; maxiter)`. If left unspecified, defaults to `run_sim`. See the help for `run_sim` for more information about the signature of the function.
* `nsims` - The number of run lengths used to estimate the target nominal property (default: 10000).
* `hmin` - The minimum value of the control limit, (default: `sqrt(eps())`).
* `maxiter` - The maximum number of bisection iterations (default: 30).
* `maxrl` - The value at which to maxrlate the run length, to avoid excessive computations (default: Inf, i.e. no maxrlation).
* `x_tol` - Absolute tolerance for the algorithm, which is terminated if
    ``h^{(k+1)} - h^{(k)} < x_{\\text{tol}}``
    (default: 1e-06)
* `f_tol` - Absolute tolerance for the algorithm, which is terminated if
    ``\\text{target}(h^{(k+1)}) - \\text{target}(h^{(k)}) < f_{\\text{tol}}``
    (default: 1.0)
* `verbose` - Whether to print information to the user about the state of the optimization (default: false).
* `parallel::Bool` - Whether the algorithm should be run in parallel, using available threads (default: false)

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Qiu, P. (2013). Introduction to Statistical Process Control. CRC Press.

"""
function bisectionCL!(CH::ControlChart, hmax; rlsim::Function = run_sim, nsims::Int = 10000, hmin::Float64 = sqrt(eps()), maxiter::Int = 30, maxrl::Real = Inf, x_tol::Float64 = 1e-06, f_tol::Float64 = 1.0, verbose::Bool = false, parallel::Bool = false)

    @assert hmin > 0 "hmin must be positive"
    @assert hmax > 0 "hmax must be positive"
    @assert maxiter > 0 "maxiter must be positive"
    @assert maxrl > 0 "maxrl must be positive"
    @assert x_tol > 0 "x_tol must be positive"
    @assert f_tol > 0 "f_tol must be positive"
    @assert nsims > 0 "nsims must be positive"


    if verbose println("Running bisection with endpoints [$(hmin), $(hmax)] ...") end

    hold = hmax + x_tol + 1.0                 # Starting value to assess convergence
    nsims_i = Int(nsims)                      # Number of simulated run lengts
    RLs = Vector{Float64}(undef, nsims_i)        # Vector of simulated run lenghts
    target = get_nominal_value(CH)                  # Target nominal ARL/QRL/...
    E_RL = 0.0                                      # Estimated ARL/QRL/...
    h = deepcopy(get_h(get_limit(CH)))              # Initialize control limit value
    conv = "Maximum number of iterations reached"
    i = 0
    while i < maxiter
        i = i+1
        h = (hmin + hmax) / 2
        if verbose print("i: $(i)/$(maxiter),\th: $(h)\t") end
        # Set control limit value
        set_limit!(CH, h)
        # Simulate run lengths 
        if parallel
            Threads.@threads for j in 1:nsims_i
                RLs[j] = first(rlsim(CH, maxiter=maxrl))
            end
        else
            for j in 1:nsims_i
                RLs[j] = first(rlsim(CH, maxiter=maxrl))
            end
        end
        # Calculate nominal measure (ARL/QRL/...) 
        E_RL = measure(RLs, CH, verbose=verbose)
        # Apply bisection algorithm
        if E_RL > target
            hmax = h
        else
            hmin = h
        end
        # Assess convergence in the run length value
        if abs(E_RL - target) < f_tol
            conv = "Convergence (target)"
            break
        end
        # Assess convergence in the control limit value
        if abs(hold - h) < x_tol
            conv = "Convergence (limit)"
            break
        end
        hold = h
    end
    return (h=h, iter=i, status = conv)
end
export bisectionCL!

"""
    bisectionCL(CH::ControlChart; kw...)

Applies the bisection algorithm to find the control limit of a control chart without modifying the control chart object `CH`.
See the documentation of `bisectionCL!` for more information about the algorithm and keyword arguments.

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Qiu, P. (2013). Introduction to Statistical Process Control. CRC Press.

"""
function bisectionCL(CH::ControlChart, hmax; rlsim::Function = run_sim, nsims::Int=1000, hmin::Float64 = sqrt(eps()), maxiter::Int = 30, maxrl::Real = Inf, x_tol::Float64 = 1e-06, f_tol::Float64 = 1.0, verbose::Bool = false)
    CH_ = shallow_copy_sim(CH)
    return bisectionCL!(CH_, hmax, rlsim=rlsim, nsims=nsims, hmin=hmin, maxiter=maxiter, maxrl=maxrl, x_tol=x_tol, f_tol=f_tol, verbose=verbose)
end
export bisectionCL


function measure(RLs, CH::ControlChart{S,L,N,P}; verbose=true) where {S,L,N<:ARL,P}
    ret = mean(RLs)
    if verbose println("E[RL] = $(ret)") end
    return ret
end

function measure(RLs, CH::ControlChart{S,L,N,P}; verbose=true) where {S,L,N<:QRL,P}
    ret = quantile(RLs, get_nominal(CH).qtl)
    if verbose println("q$(get_nominal(CH).qtl)[RL] = $(ret)") end
    return ret
end
export measure


"""
    combinedCL!(CH::ControlChart[; rlsim::Function, settings::OptSettings])

Computes the control limit to satisfy the nominal properties of a control chart, using the bisection algorithm (see for instance Qiu, 2013). The control limit upper bound `hmax` for the bisection algorithm is found using the stochastic approximation algorithm of Capizzi and Masarotto (2016)

### Inputs
* `CH` - A control chart.

### Accepted settings
* `inflate::Real` - An inflation constant for the starting control limit value so that, on average, the first iteration will move the control limit to lower values. This usually saves computational time (default: 1.05).
* `parallel::Bool` - Whether the algorithm should be run in parallel, using available threads (default: false)

#### Bisection algorithm
* `rlsim` - A function that generates a run length for the control chart with signature `rlsim(CH; maxiter)`. If left unspecified, defaults to `run_sim`. See the help for `run_sim` for more information about the signature of the function.
* `nsims` - The number of run lengths used to estimate the target nominal property (default: 10000).
* `hmin` - The minimum value of the control limit, (default: `sqrt(eps())`).
* `maxiter` - The maximum number of bisection iterations (default: 30).
* `maxrl` - The value at which to maxrlate the run length, to avoid excessive computations (default: Inf, i.e. no maxrlation).
* `x_tol` - Absolute tolerance for the algorithm, which is terminated if
    ``h^{(k+1)} - h^{(k)} < x_{\\text{tol}}``
    (default: 1e-06)
* `f_tol` - Absolute tolerance for the algorithm, which is terminated if
    ``\\text{target}(h^{(k+1)}) - \\text{target}(h^{(k)}) < f_{\\text{tol}}``
    (default: 1.0)
#### SA algorithm
* `Nfixed` - The number of iterations for the gain estimation stage (default: 200).
* `Afixed` - The fixed gain during the gain estimation stage (default: 0.1).
* `Amin` - The minimum allowed value of gain (default: 0.1).
* `Amax` - The maximum allowed value of gain (default: 100).
* `deltaSA` - The shift in control limit used during the gain estimation stage (default: 0.1).
* `q` - The power that controls the denominator in the Robbins-Monro algorithm (default: 0.55).
* `gamma` - The precision parameter for the stopping criterion of the algorithm (default: 0.05).
* `Nmin` - The minimum number of iterations required for the algorithm to end (default: 200).
* `z` - The quantile of the `Normal(0,1)` that controls the probability of the stopping criterion being satisfied (default: 3.0).
* `Cmrl` - The inflation factor for the maximum number of iterations the run length may run for (default: 10.0).
* `maxiter_sa` - Maximum number of iterations before the algorithm is forcibly ended (default: 200).
* `verbose` - Whether to print information to the user about the state of the optimization (default: false).

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Qiu, P. (2013). Introduction to Statistical Process Control. CRC Press.
* Capizzi, G., & Masarotto, G. (2016). Efficient control chart calibration by simulated stochastic approximation. IIE Transactions, 48(1), 57-65. https://doi.org/10.1080/0740817X.2015.1055392

"""
function combinedCL!(CH::ControlChart; rlsim::Function = run_sim, nsims::Int = 10000, hmin::Float64 = sqrt(eps()), maxiter::Int = 30, maxrl::Real = Inf, x_tol::Float64 = 1e-06, f_tol::Float64 = 1.0, verbose::Bool = false, inflate = 1.05, rlsim_sa::Function = run_sim_sa, Nfixed::Int = 200, Nmin::Int = 200, maxiter_sa::Int = 200, parallel::Bool = false)
    h, _, _ = saCL(CH, rlsim = rlsim_sa, Nfixed=Nfixed, Nmin=Nmin, maxiter=maxiter_sa, verbose=verbose, parallel=parallel)
    hmax = 2.0 * inflate * h
    bisectionCL!(CH, hmax, rlsim=rlsim, nsims=nsims, hmin=hmin, maxiter=maxiter, maxrl=maxrl, x_tol=x_tol, f_tol=f_tol, verbose=verbose, parallel=parallel)
end
export combinedCL!


"""
    combinedCL(CH::ControlChart; kw...)

Applies the bisection algorithm to find the control limit of a control chart without modifying the control chart object `CH`. The control limit upper bound `hmax` for the bisection algorithm is found using the stochastic approximation algorithm of Capizzi and Masarotto (2016).
See the documentation of `combinedCL!` for more information about the algorithm and keyword arguments.

### Keyword arguments:
* See the documentation of `combinedCL!` for a list of keyword arguments.

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Qiu, P. (2013). Introduction to Statistical Process Control. CRC Press.
* Capizzi, G., & Masarotto, G. (2016). Efficient control chart calibration by simulated stochastic approximation. IIE Transactions, 48(1), 57-65. https://doi.org/10.1080/0740817X.2015.1055392
"""
function combinedCL(CH::ControlChart; kw...)
    CH_ = shallow_copy_sim(CH)
    return combinedCL!(CH_; kw...)
end
export combinedCL
