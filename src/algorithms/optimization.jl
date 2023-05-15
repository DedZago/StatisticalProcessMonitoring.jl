function optimize_parameter!(CH; settings)
    CH_ = shallow_copy_sim(CH)
    @unpack nsims_opt, trace, method_opt = settings

    function rlconstr(par::Vector, grad::Vector)::Float64
        set_parameter!(CH_, par)
        optimizeLimit!(CH_, settings)
        if trace_opt > 0
            print("$(round.(par, digits=6))\t")
        end
        ret = SPM.measure([rlsim_oc(CH_) for _ in 1:nsims_opt], CH_, verbose=trace_opt > 0)
        return ret
    end

    #FIXME: from here
    if method_opt == :Grid
        set_parameter!(CH, optimize_grid(CH, rlconstr, minpar, maxpar, x_tol, maxiter))
    elseif method_opt == :SPSA
        set_parameter!(CH, optimize_SPSA(CH, rlconstr, minpar, maxpar, x_tol, maxiter))
    else
        set_parameter!(CH, optimize_nlopt(CH, rlconstr, minpar, maxpar, x_tol, maxiter, method_opt))
    end
    set_h!(get_limit(CH), get_h(get_limit(CH_)))
    return get_parameter(CH)
end
export optimize_parameter!
