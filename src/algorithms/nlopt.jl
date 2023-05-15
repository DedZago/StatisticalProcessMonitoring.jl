using NLopt

"""
    optimize_nlopt(CH::ControlChart, rlconstr::Function, settings::OptSettings)

#FIXME: docstring
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