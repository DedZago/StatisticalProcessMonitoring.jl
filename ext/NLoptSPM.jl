module NLoptSPM

using StatisticalProcessMonitoring, NLopt

"""
    optimize_nlopt(CH::ControlChart, rlconstr::Function, settings::OptSettings)

Optimizes the Control Chart design parameter using the NLOpt library.

### Args
* `CH::ControlChart`: The ControlChart object whose parameters must be optimized.
* `rlconstr::Function`: The objective function to be minimized.
* `settings::OptSettings`: The settings for the optimization process. Includes:
    - `minpar`: The lower bounds for the parameters.
    - `maxpar`: The upper bounds for the parameters.
    - `x_tol`: The relative tolerance for convergence.
    - `maxiter`: The maximum number of iterations.
    - `optimizer`: The optimization method to be used.

### Returns
* `minx::Vector{Float64}`: The set of optimal parameters for the control chart.
"""
function StatisticalProcessMonitoring.optimize_nlopt(CH::ControlChart, rlconstr::Function, settings::OptSettings; optimizer::Symbol = :LN_BOBYQA)
    minpar = settings.minpar
    maxpar = settings.maxpar
    x_tol = settings.x_tol
    maxiter = settings.maxiter

    # @show minpar, maxpar
    par0 = collect(get_design(CH))
    opt = NLopt.Opt(optimizer, length(par0))
    opt.lower_bounds = minpar
    opt.upper_bounds = maxpar
    opt.xtol_rel = x_tol
    opt.maxeval = maxiter
    opt.min_objective = rlconstr
    # @show par0, opt
    (minf, minx, ret) = NLopt.optimize(opt, par0)
    return minx
end

end # module
