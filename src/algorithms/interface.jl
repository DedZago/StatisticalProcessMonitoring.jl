"""
    optimize_limit!(CH::ControlChart[; settings = OptSettings()])

Optimizes the control limit of a ControlChart object.

### Args
* CH (ControlChart): The ControlChart object to optimize.
* settings (OptSettings, optional): Optimization settings. Defaults to OptSettings().

### Returns
    The optimized control limit value.

### Raises
    ValueError: If the optimization method specified in settings is unknown.

### Example
    optimize_limit!(my_chart, settings=OptSettings(ic_solver=:SA))
"""
function optimize_limit!(CH::ControlChart, solver::Symbol=:SA; kw...)
    if solver == :SA
        return saCL!(CH; kw...)
    elseif solver == :Bisection
        return bisectionCL!(CH; kw...)
    elseif solver == :Combined
        combinedCL!(CH; kw...)
    else
        error("Unknown optimization method.")
    end
end
export optimize_limit!

"""
    optimize_limit(CH::ControlChart[; settings = OptSettings()])

Optimizes the control limit of a ControlChart object, without modifying the original ControlChart object.

### Args
* CH (ControlChart): The ControlChart object to optimize.
* settings (OptSettings, optional): Optimization settings. Defaults to OptSettings().

### Returns
    The optimized control limit value.

### Raises
    ValueError: If the optimization method specified in settings is unknown.

### Example
    optimize_limit(my_chart, settings=OptSettings(ic_solver=:SA))
"""
function optimize_limit(CH::ControlChart, solver::Symbol = :SA; kw...)
    CH_ = deepcopy(CH)
    optimize_limit!(CH_, solver=solver, kw...)
end
export optimize_limit


"""
    optimize_design!(CH, rlsim_oc[; settings = OptSettings()])

Optimizes the parameter of a simulation `CH` with respect to a given objective function `rlsim_oc`. 

### Arguments
- `CH` : The simulation to optimize.
- `rlsim_oc` : The objective function.
- `settings` (optional, default=OptSettings()) : Optimization settings.

### Returns
- `get_design(CH)` : The optimized parameter.
"""
function optimize_design!(CH, rlsim_oc; settings::OptSettings = OptSettings())
    CH_ = shallow_copy_sim(CH)
    @unpack nsims_opt, trace, method_opt = settings

    function rlconstr(par::Vector, grad::Vector)::Float64
        set_design!(CH_, par)
        optimize_limit!(CH_, settings=settings)
        if trace > 0
            print("$(round.(par, digits=6))\t")
        end
        ret = SPM.measure([rlsim_oc(CH_) for _ in 1:nsims_opt], CH_, verbose=trace > 0)
        return ret
    end

    if method_opt == :Grid
        set_design!(CH, optimize_grid(CH, rlconstr, settings))
    elseif method_opt == :SPSA
        #TODO: implement SPSA
        set_design!(CH, optimize_SPSA(CH, rlconstr, settings))
    else
        set_design!(CH, optimize_nlopt(CH, rlconstr, settings))
    end
    set_h!(get_limit(CH), get_h(get_limit(CH_)))
    return get_design(CH)
end
export optimize_design!

"""
    optimize_design(CH, rlsim_oc[; settings = OptSettings()])

Optimize a parameter using the specified CH and rlsim_oc.

### Args
    `CH`: the CH parameter.
    `rlsim_oc`: the rlsim_oc parameter.
    `settings`: the optimization settings.

### Returns
    The optimized parameter values.
"""
function optimize_design(CH, rlsim_oc; settings::OptSettings = OptSettings())
    CH_ = deepcopy(CH)
    optimize_design!(CH_, rlsim_oc, settings=settings)
end
export optimize_design