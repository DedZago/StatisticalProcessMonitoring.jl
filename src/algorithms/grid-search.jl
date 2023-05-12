using NLopt

#FIXME: grid search algorithm to optimize the parameters of the control chart
function optimize_parameter!(CH; rlsim::Function = run_sim_sa, minpar::Vector{Float64} = fill(-Inf, length(get_parameter(CH))), maxpar::Vector{Float64} = fill(Inf, length(get_parameter(CH))), x_tol::Float64 = 1e-05, rlsim_oc::Function, nsims::Real = 10000, maxiter::Real = 500,  method = :BOBYQA, ic_solver::Symbol = :SA)
    CH_ = shallow_copy_sim(CH)

    function rlconstr(par::Vector)::Float64
        set_parameter!(CH_, par)
        if ic_solver == :SA
            saCL!(CH_, rlsim=rlsim)
        elseif ic_solver == :Combined
            combinedCL!(CH_, rlsim=rlsim)
        end
        return SPM.measure([rlsim_oc(CH_) for _ in 1:nsims], CH_)
    end

    if method == :BOBYQA
        set_parameter!(CH, optimize_bobyqa(CH, rlconstr, minpar, maxpar, x_tol, maxiter))
    elseif method == :Grid
        set_parameter!(CH, optimize_grid(CH, rlconstr, minpar, maxpar, x_tol, maxiter))
    end
end
export optimize_parameter!

function optimize_grid(CH, rlconstr::Function, minpar::Vector{Float64}, maxpar::Vector{Float64}, x_tol::Float64, maxiter::Real)
    error("Not implemented yet")
end
export optimize_grid

function optimize_bobyqa(CH, rlconstr::Function, minpar::Vector{Float64}, maxpar::Vector{Float64}, x_tol::Float64, maxiter::Real)
    par0 = collect(get_parameter(CH))
    opt = NLopt.Opt(:LN_BOBYQA, length(par0))
    opt.lower_bounds = minpar
    opt.upper_bounds = maxpar
    opt.xtol_rel = x_tol
    opt.maxeval = maxiter
    opt.min_objective = rlconstr
    (minf, minx, ret) = NLopt.optimize(opt, par0)
    return minx
end
export optimize_bobyqa