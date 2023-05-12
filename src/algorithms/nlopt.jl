using NLopt

function optimize_nlopt(CH, rlconstr::Function, minpar::Vector{Float64}, maxpar::Vector{Float64}, x_tol::Float64, maxiter::Real, method::Symbol = LN_BOBYQA)
    par0 = collect(get_parameter(CH))
    opt = NLopt.Opt(method, length(par0))
    opt.lower_bounds = minpar
    opt.upper_bounds = maxpar
    opt.xtol_rel = x_tol
    opt.maxeval = maxiter
    opt.min_objective = rlconstr
    (minf, minx, ret) = NLopt.optimize(opt, par0)
    return minx
end
export optimize_nlopt