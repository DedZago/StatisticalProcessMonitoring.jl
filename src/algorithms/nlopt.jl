using NLopt

"""
    optimize_nlopt(CH::ControlChart, rlconstr::Function, settings::OptSettings)

Optimizes the Control Chart design parameter using the NLOpt library.

### Args
* `CH`: The ControlChart object to be optimized.
* `rlconstr`: The objective function to be minimized.
* `settings`: The settings for the optimization process. Includes:
    - `minpar_opt`: The lower bounds for the parameters.
    - `maxpar_opt`: The upper bounds for the parameters.
    - `x_tol_opt`: The relative tolerance for convergence.
    - `maxiter_opt`: The maximum number of iterations.
    - `method_opt`: The optimization method to be used.

### Returns
* `minx::Vector{Float64}`: The set of optimal parameters for the Control Chart.
"""
function optimize_nlopt(CH::ControlChart, rlconstr::Function, settings::OptSettings)
    @unpack minpar_opt, maxpar_opt, x_tol_opt, maxiter_opt, method_opt = settings
    par0 = collect(get_parameter(CH))
    opt = NLopt.Opt(method_opt, length(par0))
    opt.lower_bounds = minpar_opt
    opt.upper_bounds = maxpar_opt
    opt.xtol_rel = x_tol_opt
    opt.maxeval = maxiter_opt
    opt.min_objective = rlconstr
    (minf, minx, ret) = NLopt.optimize(opt, par0)
    return minx
end
export optimize_nlopt