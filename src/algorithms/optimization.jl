#FIXME: test
"""
    optimize_parameter!(CH, rlsim_oc[; settings = OptSettings()])

Optimizes the parameter of a simulation `CH` with respect to a given objective function `rlsim_oc`. 

#### Arguments
- `CH` : The simulation to optimize.
- `rlsim_oc` : The objective function.
- `settings` (optional, default=OptSettings()) : Optimization settings.

#### Returns
- `get_parameter(CH)` : The optimized parameter.
"""
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

"""
    optimize_parameter(CH, rlsim_oc[; settings = OptSettings()])

Optimize a parameter using the specified CH and rlsim_oc.

### Args:
    `CH`: the CH parameter.
    `rlsim_oc`: the rlsim_oc parameter.
    `settings`: the optimization settings.

### Returns:
    The optimized parameter values.

"""
function optimize_parameter(CH, rlsim_oc; settings = OptSettings())
    CH_ = deepcopy(CH)
    optimize_parameter!(CH_, rlsim_oc, settings=settings)
end
export optimize_parameter
