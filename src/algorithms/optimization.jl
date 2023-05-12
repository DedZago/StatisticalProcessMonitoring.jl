function optimize_parameter!(CH; rlsim::Function = run_sim_sa, minpar::Vector{Float64} = fill(-Inf, length(get_parameter(CH))), maxpar::Vector{Float64} = fill(Inf, length(get_parameter(CH))), x_tol::Float64 = 1e-05, rlsim_oc::Function, nsims::Real = 10000, maxiter::Real = 500,  method::Symbol = :LN_BOBYQA, ic_solver::Symbol = :SA, trace::Int = 0)
    CH_ = shallow_copy_sim(CH)

    function rlconstr(par::Vector, grad::Vector)::Float64
        set_parameter!(CH_, par)
        if ic_solver == :SA
            saCL!(CH_, rlsim=rlsim, verbose=false)
        elseif ic_solver == :Combined
            combinedCL!(CH_, rlsim=rlsim, verbose=false)
        end
        if trace > 0
            print("$(round.(par, digits=6))\t")
        end
        ret = SPM.measure([rlsim_oc(CH_) for _ in 1:nsims], CH_, verbose=trace > 0)
        return ret
    end

    if method == :Grid
        set_parameter!(CH, optimize_grid(CH, rlconstr, minpar, maxpar, x_tol, maxiter))
    elseif method == :SPSA
        set_parameter!(CH, optimize_SPSA(CH, rlconstr, minpar, maxpar, x_tol, maxiter))
    else
        set_parameter!(CH, optimize_nlopt(CH, rlconstr, minpar, maxpar, x_tol, maxiter, method))
    end
    set_h!(get_limit(CH), get_h(get_limit(CH_)))
    return get_parameter(CH)
end
export optimize_parameter!
