using Statistics

"""
    bisectionCL!(CH::ControlChart; hmax, kw...)

Computes the control limit to satisfy the nominal properties of a control chart, using the bisection algorithm (see for instance Qiu, 2013)

### Inputs
* `CH` - A control chart.
* `kw...` - Keyword arguments that control the behaviour of the algorithm. 

### Keyword arguments:
* `rlsim` - A function that generates new data with signature `rlsim(CH; maxiter)`. If left unspecified, defaults to `run_sim`.
* `hmin` - The minimum value of the control limit, defaults to `sqrt(eps())`.
* `hmax` - The maximum value for the control limit.
* `maxiter` - The maximum number of bisection iterations.
* `nsims` - The number of run lengths used to estimate the target nominal property.
* `trunc` - The maximum run length after which it is truncated, to avoid excessive computations.
* `x_tol` - Absolute tolerance for the algorithm, which is ended if
    ``h^{(k+1)} - h^{(k)} < x_{\\text{tol}}``
* `f_tol` - Absolute tolerance for the algorithm, which is ended if
    ``\\text{target}(h^{(k+1)}) - \\text{target}(h^{(k)}) < f_{\\text{tol}}``

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Qiu, P. (2013). Introduction to Statistical Process Control. CRC Press.

"""
function bisectionCL!(CH::ControlChart; rlsim::Function = run_sim, hmin::Float64 = sqrt(eps()), hmax::Float64, maxiter::Int = 50, nsims::Real = 10000, trunc::Float64 = 20*get_nominal_value(CH), x_tol::Float64 = 1e-06, f_tol::Float64 = 1.0, verbose::Bool=true)

    #TODO: consider truncation of the control chart run lengths

    @assert hmin > 0 "hmin must be positive"
    @assert hmax > 0 "hmax must be positive"
    @assert maxiter > 0 "maxiter must be positive"
    @assert trunc > 0 "trunc must be positive"
    @assert x_tol > 0 "x_tol must be positive"
    @assert f_tol > 0 "f_tol must be positive"


    if verbose println("Running bisection with endpoints [$(hmin), $(hmax)] ...") end

    hold = hmax + 1
    nsims_i = Int(nsims)
    RLs = Vector{Float64}(undef, nsims_i)
    target = get_nominal_value(CH)
    E_RL = 0.0
    h = deepcopy(get_h(get_limit(CH)))
    conv = "Maximum number of iterations reached"
    i = 0
    while i < maxiter
        i = i+1
        h = (hmin + hmax) / 2
        if verbose print("i: $(i)/$(maxiter),\th: $(h)\t") end
        set_limit!(CH, h)
        for j in 1:nsims_i
            RLs[j] = first(rlsim(CH, maxiter=trunc))
        end
        E_RL = measure(RLs, CH, verbose)
        if E_RL > target
            hmax = h
        else
            hmin = h
        end
        if abs(E_RL - target) < f_tol
            conv = "Convergence (target)"
            break
        end
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
function bisectionCL(CH::ControlChart; kw...)
    CH_ = shallow_copy_sim(CH)
    return bisectionCL!(CH_; kw...)
end
export bisectionCL


function measure(RLs, CH::ControlChart{S,L,N,P}, verbose) where {S,L,N<:ARL,P}
    ret = mean(RLs)
    if verbose println("E[RL] = $(ret)") end
    return ret
end

function measure(RLs, CH::ControlChart{S,L,N,P}, verbose) where {S,L,N<:QRL,P}
    ret = quantile(RLs, get_nominal(CH).qtl)
    if verbose println("q$(get_nominal(CH).qtl)[RL] = $(ret)") end
    return ret
end


"""
    combinedCL!(CH::ControlChart; kw...)

Computes the control limit to satisfy the nominal properties of a control chart, using the bisection algorithm (see for instance Qiu, 2013). The control limit upper bound `hmax` for the bisection algorithm is found using the stochastic approximation algorithm of Capizzi and Masarotto (2016)

### Inputs
* `CH` - A control chart.
* `kw...` - Keyword arguments that control the behaviour of the algorithm. 

### Keyword arguments:
* `rlsim` - A function that generates new data with signature `rlsim(CH; maxiter)`. If left unspecified, defaults to `run_sim`.
* `hmin` - The minimum value of the control limit, defaults to `sqrt(eps())`.
* `hmax` - The maximum value for the control limit.
* `maxiter` - The maximum number of bisection iterations.
* `nsims` - The number of run lengths used to estimate the target nominal property.
* `trunc` - The maximum run length after which it is truncated, to avoid excessive computations.
* `x_tol` - Absolute tolerance for the algorithm, which is ended if
    ``h^{(k+1)} - h^{(k)} < x_{\\text{tol}}``
* `f_tol` - Absolute tolerance for the algorithm, which is ended if
    ``\\text{target}(h^{(k+1)}) - \\text{target}(h^{(k)}) < f_{\\text{tol}}``

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Qiu, P. (2013). Introduction to Statistical Process Control. CRC Press.
* Capizzi, G., & Masarotto, G. (2016). Efficient control chart calibration by simulated stochastic approximation. IIE Transactions, 48(1), 57-65. https://doi.org/10.1080/0740817X.2015.1055392

"""
function combinedCL!(CH::ControlChart; rlsim::Function = run_sim_sa, inflate::Float64 = 1.25, Nfixed::Int = 200, Nmin::Int = 200, maxiter_sa::Int = 200, q = 0.55, hmin::Float64 = sqrt(eps()), maxiter_bi::Int = 50, nsims_bi::Real = 10000, x_tol::Float64 = 1e-06, f_tol::Float64 = 1.0,trunc::Float64 = 20*get_nominal_value(CH), verbose::Bool = true)
    h, _, _ = saCL(CH, rlsim = rlsim, Nfixed = Nfixed, Nmin = Nmin, maxiter = maxiter_sa, q=q, verbose=verbose)
    bisectionCL!(CH, rlsim = rlsim, hmin = hmin, hmax = inflate*h, maxiter = maxiter_bi, nsims = nsims_bi, trunc=trunc, x_tol = x_tol, f_tol = f_tol, verbose=verbose)
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
