using NLopt

"""
    optimize_nlopt(CH::ControlChart, rlconstr::Function, settings::OptSettings)

Optimizes the Control Chart design parameter using the NLOpt library.

### Args
* `CH`: The ControlChart object to be optimized.
* `rlconstr`: The objective function to be minimized.
* `settings`: The settings for the optimization process. Includes:
    - `minpar`: The lower bounds for the parameters.
    - `maxpar`: The upper bounds for the parameters.
    - `x_tol`: The relative tolerance for convergence.
    - `maxiter`: The maximum number of iterations.
    - `optimizer`: The optimization method to be used.

### Returns
* `minx::Vector{Float64}`: The set of optimal parameters for the Control Chart.
"""
function optimize_nlopt(CH::ControlChart, rlconstr::Function, settings::OptSettings; optimizer::Symbol = :LN_BOBYQA)
    @unpack minpar, maxpar, x_tol, maxiter = settings
    @show minpar, maxpar
    par0 = collect(get_design(CH))
    opt = NLopt.Opt(optimizer, length(par0))
    opt.lower_bounds = minpar
    opt.upper_bounds = maxpar
    opt.xtol_rel = x_tol
    opt.maxeval = maxiter
    opt.min_objective = rlconstr
    (minf, minx, ret) = NLopt.optimize(opt, par0)
    return minx
end
export optimize_nlopt