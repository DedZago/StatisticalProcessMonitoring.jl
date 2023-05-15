using Statistics

"""
    bisectionCL!(CH::ControlChart[; rlsim::Function, settings::OptSettings])

Computes the control limit to satisfy the nominal properties of a control chart, using the bisection algorithm (see for instance Qiu, 2013)

### Inputs
* `CH` - A control chart.
* `rlsim` - A function that generates a run length for the control chart with signature `rlsim(CH; maxiter)`. If left unspecified, defaults to `run_sim`. See the help for `run_sim` for more information about the signature of the function.
* `settings` - An `OptSettings` objects which contains variables that control the behaviour of the algorithm. See the `Accepted settings` section below for information about the settings that control the behaviour of the algorithm. For more information about the specifics of each keyword argument, see for instance Qiu (2013).

### Accepted settings:
* `hmin_bi` - The minimum value of the control limit, defaults to `sqrt(eps())`.
* `hmax_bi` - The maximum value for the control limit.
* `maxiter_bi` - The maximum number of bisection iterations.
* `nsims_bi` - The number of run lengths used to estimate the target nominal property.
* `trunc_bi` - The maximum run length after which it is trunc_biated, to avoid excessive computations.
* `x_tol_bi` - Absolute tolerance for the algorithm, which is ended if
    ``h^{(k+1)} - h^{(k)} < x_{\\text{tol}}``
* `f_tol_bi` - Absolute tolerance for the algorithm, which is ended if
    ``\\text{target}(h^{(k+1)}) - \\text{target}(h^{(k)}) < f_{\\text{tol}}``

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Qiu, P. (2013). Introduction to Statistical Process Control. CRC Press.

"""
function bisectionCL!(CH::ControlChart; rlsim::Function = run_sim, settings::OptSettings = OptSettings())

    #TODO: consider trunc_biation of the control chart run lengths

    @unpack rlsim, hmin_bi, hmax_bi, maxiter_bi, nsims_bi, trunc_bi, x_tol_bi, f_tol_bi, verbose_bi = settings

    @assert hmin_bi > 0 "hmin_bi must be positive"
    @assert hmax_bi > 0 "hmax_bi must be positive"
    @assert maxiter_bi > 0 "maxiter_bi must be positive"
    @assert trunc_bi > 0 "trunc_bi must be positive"
    @assert x_tol_bi > 0 "x_tol_bi must be positive"
    @assert f_tol_bi > 0 "f_tol_bi must be positive"
    @assert nsims_bi > 0 "nsims_bi must be positive"


    if verbose_bi println("Running bisection with endpoints [$(hmin_bi), $(hmax_bi)] ...") end

    hold = hmax_bi + 1
    nsims_bi_i = Int(nsims_bi)
    RLs = Vector{Float64}(undef, nsims_bi_i)
    target = get_nominal_value(CH)
    E_RL = 0.0
    h = deepcopy(get_h(get_limit(CH)))
    conv = "Maximum number of iterations reached"
    i = 0
    while i < maxiter_bi
        i = i+1
        h = (hmin_bi + hmax_bi) / 2
        if verbose_bi print("i: $(i)/$(maxiter_bi),\th: $(h)\t") end
        set_limit!(CH, h)
        for j in 1:nsims_bi_i
            RLs[j] = first(rlsim(CH, maxiter=trunc_bi))
        end
        E_RL = measure(RLs, CH, verbose=verbose_bi)
        if E_RL > target
            hmax_bi = h
        else
            hmin_bi = h
        end
        if abs(E_RL - target) < f_tol_bi
            conv = "Convergence (target)"
            break
        end
        if abs(hold - h) < x_tol_bi
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
function bisectionCL(CH::ControlChart; rlsim::Function = run_sim, settings::OptSettings = OptSettings())
    CH_ = shallow_copy_sim(CH)
    return bisectionCL!(CH_; rlsim = rlsim, settings = settings)
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


"""
    combinedCL!(CH::ControlChart[; rlsim::Function, settings::OptSettings])

Computes the control limit to satisfy the nominal properties of a control chart, using the bisection algorithm (see for instance Qiu, 2013). The control limit upper bound `hmax_bi` for the bisection algorithm is found using the stochastic approximation algorithm of Capizzi and Masarotto (2016)

### Inputs
* `CH` - A control chart.
* `rlsim` - A function that generates a run length for the control chart with signature `rlsim(CH; maxiter)`. If left unspecified, defaults to `run_sim`. See the help for `run_sim` for more information about the signature of the function.
* `settings` - An `OptSettings` objects which contains variables that control the behaviour of the algorithm. See the `Accepted settings` section below for information about the settings that control the behaviour of the algorithm. For more information about the specifics of each keyword argument, see for instance Qiu (2013).

### Accepted settings
#### Bisection algorithm
* `rlsim` - A function that generates new data with signature `rlsim(CH; maxiter_bi)`. If left unspecified, defaults to `run_sim`.
* `hmin_bi` - The minimum value of the control limit, defaults to `sqrt(eps())`.
* `hmax_bi` - The maximum value for the control limit.
* `maxiter_bi` - The maximum number of bisection iterations.
* `nsims_bi` - The number of run lengths used to estimate the target nominal property.
* `trunc_bi` - The maximum run length after which it is trunc_biated, to avoid excessive computations.
* `x_tol_bi` - Absolute tolerance for the algorithm, which is ended if
    ``h^{(k+1)} - h^{(k)} < x_{\\text{tol}}``
* `f_tol_bi` - Absolute tolerance for the algorithm, which is ended if
    ``\\text{target}(h^{(k+1)}) - \\text{target}(h^{(k)}) < f_{\\text{tol}}``
#### SA algorithm
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
* Qiu, P. (2013). Introduction to Statistical Process Control. CRC Press.
* Capizzi, G., & Masarotto, G. (2016). Efficient control chart calibration by simulated stochastic approximation. IIE Transactions, 48(1), 57-65. https://doi.org/10.1080/0740817X.2015.1055392

"""
function combinedCL!(CH::ControlChart; rlsim::Function = run_sim_sa, settings::OptSettings = OptSettings(Nfixed_sa = 200, Nmin_sa = 200, maxiter_sa = 200))
    h, _, _ = saCL(CH, settings=settings)
    bisectionCL!(CH, settings = OptSettings(settings, hmax_bi = settings.inflate_bi * 2.0 * h))
end
export combinedCL!


"""
    combinedCL(CH::ControlChart; kw...)

Applies the bisection algorithm to find the control limit of a control chart without modifying the control chart object `CH`. The control limit upper bound `hmax_bi` for the bisection algorithm is found using the stochastic approximation algorithm of Capizzi and Masarotto (2016).
See the documentation of `combinedCL!` for more information about the algorithm and keyword arguments.

### Keyword arguments:
* See the documentation of `combinedCL!` for a list of keyword arguments.

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Qiu, P. (2013). Introduction to Statistical Process Control. CRC Press.
* Capizzi, G., & Masarotto, G. (2016). Efficient control chart calibration by simulated stochastic approximation. IIE Transactions, 48(1), 57-65. https://doi.org/10.1080/0740817X.2015.1055392
"""
function combinedCL(CH::ControlChart; rlsim::Function = run_sim_sa, settings::OptSettings = OptSettings(Nfixed_sa = 200, Nmin_sa = 200, maxiter_sa = 200))
    CH_ = shallow_copy_sim(CH)
    return combinedCL!(CH_, rlsim=rlsim, settings=settings)
end
export combinedCL
