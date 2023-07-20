using Statistics
# FIXME: test double bootstrap functions

"""
    doubleBootstrap!(CH::ControlChart[; rlsim::Function, settings::OptSettings])

Computes the control limit to satisfy the nominal properties of a control chart, using the bisection algorithm on bootstrapped paths (see for instance Qiu, 2013).

### Inputs
* `CH` - A control chart.
* `rlsim` - A function that generates a path of the control chart statistic with signature `rlsim(CH; maxiter)`. If left unspecified, defaults to `run_path_sim`. See the help for `run_path_sim` for more information about the signature of the function.
* `settings` - An `OptSettings` objects which contains variables that control the behaviour of the algorithm. See the `Accepted settings` section below for information about the settings that control the behaviour of the algorithm. For more information about the specifics of each keyword argument, see for instance Qiu (2013).

### Accepted settings:
* `maxiter_bi` - The maximum number of bisection iterations.
* `nsims_bi` - The number of run lengths used to estimate the target nominal property.
* `trunc_bi` - The maximum run length after which the run length is truncated, to avoid excessive computations.
* `x_tol_bi` - Absolute tolerance for the algorithm, which is ended if
    ``h^{(k+1)} - h^{(k)} < x_{\\text{tol}}``
* `f_tol_bi` - Absolute tolerance for the algorithm, which is ended if
    ``\\text{target}(h^{(k+1)}) - \\text{target}(h^{(k)}) < f_{\\text{tol}}``

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Qiu, P. (2013). Introduction to Statistical Process Control. CRC Press.

"""
function doubleBootstrap!(CH::ControlChart; settings::OptSettings = OptSettings(ic_solver=:Double))

    #TODO: consider trunc_bi of the control chart run lengths when calculating the control limit?

    @unpack rlsim, maxiter_bi, nsims_bi, trunc_bi, x_tol_bi, f_tol_bi, verbose_bi = settings

    @assert maxiter_bi > 0 "maxiter_bi must be positive"
    @assert nsims_bi > 0 "nsims_bi must be positive"
    @assert trunc_bi > 0 "trunc_bi must be positive"
    @assert x_tol_bi > 0 "x_tol_bi must be positive"
    @assert f_tol_bi > 0 "f_tol_bi must be positive"

    tmp_rlpath = rlsim(CH, maxiter=2)
    @assert isa(tmp_rlpath, Vector) "rlsim must be a vector"
    @assert length(tmp_rlpath) == 2 "rlsim must be of length maxiter"

    if verbose_bi println("Generating $(nsims_bi) run length paths ...") end

    nsims_bi_i = Int(nsims_bi)                              # Number of simulated run lengts
    maxrl = Int(min(get_maxrl(CH), 10*get_nominal_value(CH)))

    rl_paths = Matrix{Float64}(undef, nsims_bi_i, maxrl)    # Generated run length paths
    for i in 1:nsims_bi_i
        rl_paths[i, :] = rlsim(CH, maxiter = maxrl)
    end

    hmin_bi, hmax_bi = extrema(rl_paths)
    if verbose_bi println("Running bisection on simulated paths with endpoints [$(hmin_bi), $(hmax_bi)] ...") end

    hold = hmax_bi + x_tol_bi + 1.0                 # Starting value to assess convergence
    RLs = Vector{Float64}(undef, nsims_bi_i)        # Vector of simulated run lenghts
    target = get_nominal_value(CH)                  # Target nominal ARL/QRL/...
    E_RL = 0.0                                      # Estimated ARL/QRL/...
    h = 0.0                                         # Initialize control limit value
    conv = "Maximum number of iterations reached"
    i = 0
    idx = Vector{Int}(undef, nsims_bi_i)
    rows = 1:nsims_bi_i
    while i < maxiter_bi
        i = i+1
        h = (hmin_bi + hmax_bi) / 2
        if verbose_bi print("i: $(i)/$(maxiter_bi),\th: $(h)\t") end

        idx .= sample(rows, nsims_bi_i)
        # Calculate run length on simulated paths
        for j in 1:nsims_bi_i
            for k in 1:maxrl
                #FIXME: use set_h! and set_value! to make it more general to double-sided limits etc...
                if rl_paths[idx[j],k] > h
                    RLs[j] = k
                    break
                end 
                if k == maxrl
                    RLs[j] = maxrl
                end
            end
        end
        # Calculate nominal measure (ARL/QRL/...) 
        E_RL = measure(RLs, CH, verbose=verbose_bi)
        # Apply bisection algorithm
        if E_RL > target
            hmax_bi = h
        else
            hmin_bi = h
        end
        # Assess convergence in the run length value
        if abs(E_RL - target) < f_tol_bi
            conv = "Convergence (target)"
            break
        end
        # Assess convergence in the control limit value
        if abs(hold - h) < x_tol_bi
            conv = "Convergence (limit)"
            break
        end
        hold = h
    end
    return (h=h, iter=i, status = conv)
end
export doubleBootstrap!

"""
    doubleBootstrap(CH::ControlChart; kw...)

Applies the bisection algorithm on simulated run length paths to find the control limit of a control chart without modifying the control chart object `CH`.
See the documentation of `doubleBootstrap!` for more information about the algorithm and keyword arguments.

### Returns
* A `NamedTuple` containing the estimated control limit `h`, the total number of iterations `iter`, and information `status` about the convergence of the algorithm.

### References
* Qiu, P. (2013). Introduction to Statistical Process Control. CRC Press.

"""
function doubleBootstrap(CH::ControlChart; settings::OptSettings = OptSettings(ic_solver=:Double))
    CH_ = shallow_copy_sim(CH)
    return doubleBootstrap!(CH_; settings = settings)
end
export doubleBootstrap

