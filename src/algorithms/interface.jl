#FIXME: test
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
function optimize_limit!(CH::ControlChart; settings::OptSettings = OptSettings())
    @unpack ic_solver, rlsim = settings
    if ic_solver == :SA
        return saCL!(CH, settings = settings)
    elseif ic_solver == :Bisection
        return bisectionCL!(CH, settings = settings)
    elseif ic_solver == :Combined
        combinedCL!(CH, settings = settings)
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
function optimize_limit(CH::ControlChart; settings::OptSettings = OptSettings())
    CH_ = deepcopy(CH)
    optimize_limit!(CH_, settings=settings)
end
export optimize_limit