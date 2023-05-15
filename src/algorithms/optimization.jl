#FIXME: test
function optimize_parameter!(CH, rlsim_oc; settings = OptSettings())
    CH_ = shallow_copy_sim(CH)
    @unpack nsims_opt, trace, method_opt = settings

    function rlconstr(par::Vector, grad::Vector)::Float64
        set_parameter!(CH_, par)
        optimize_limit!(CH_, settings=settings)
        if trace > 0
            print("$(round.(par, digits=6))\t")
        end
        ret = SPM.measure([rlsim_oc(CH_) for _ in 1:nsims_opt], CH_, verbose=trace > 0)
        return ret
    end

    if method_opt == :Grid
        set_parameter!(CH, optimize_grid(CH, rlconstr, settings))
    elseif method_opt == :SPSA
        #TODO: implement SPSA
        set_parameter!(CH, optimize_SPSA(CH, rlconstr, settings))
    else
        set_parameter!(CH, optimize_nlopt(CH, rlconstr, settings))
    end
    set_h!(get_limit(CH), get_h(get_limit(CH_)))
    return get_parameter(CH)
end
export optimize_parameter!
