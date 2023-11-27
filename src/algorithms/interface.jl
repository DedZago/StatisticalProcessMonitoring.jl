"""
    optimize_limit(CH::ControlChart, solver = :Bootstrap; hmax = 20.0, kw...)

Optimizes the control limit of a ControlChart object.

### Arguments
- `CH::ControlChart`: The ControlChart object to optimize.
- `solver::Symbol`: The solver algorithm to use (default: `:Bootstrap`).

### Keyword Arguments
- `hmax::Float64`: The maximum value of the control limit. Only used for the bisection algorithm (default: 100.0)
- `kw...`: Additional keyword arguments to pass to the algorithm.

### Returns
    The optimized control limit value.

### Raises
    ValueError: If the optimization method specified in settings is unknown.
"""
function optimize_limit!(CH::ControlChart, solver::Symbol=:Bootstrap; hmax::Float64 = 100.0, kw...)
    valid_solvers = [:SA, :Bisection, :Combined, :Bootstrap]
    @assert solver in valid_solvers "Unknown control limit solver. Valid solvers are $(valid_solvers)"

    if solver == :SA
        return saCL!(CH; kw...)
    elseif solver == :Bisection
        return bisectionCL!(CH, hmax; kw...)
    elseif solver == :Combined
        return combinedCL!(CH; kw...)
    elseif solver == :Bootstrap
        return bootstrapCL!(CH; kw...)
    end
end
export optimize_limit!

"""
    optimize_limit(CH::ControlChart, solver::Symbol = :Bootstrap; kw...)

Optimizes the control limit of a ControlChart object, without modifying the original ControlChart object.

### Arguments
- `CH::ControlChart`: The ControlChart object to optimize.
- `solver::Symbol`: The solver algorithm to use (default: `:Bootstrap`).

### Keyword Arguments
- `hmax::Float64`: The maximum value of the control limit. Used only for the bisection algorithm (default: 100.0)
- `kw...`: Additional keyword arguments to pass to the algorithm.

### Returns
    The optimized control limit value.

### Raises
    ValueError: If the optimization method specified in settings is unknown.

### Example
    optimize_limit(my_chart, settings=OptSettings(ic_solver=:SA))
"""
function optimize_limit(CH::ControlChart, solver::Symbol = :Bootstrap; hmax::Float64 = 100.0, kw...)
    CH_ = deepcopy(CH)
    optimize_limit!(CH_, solver; hmax = hmax, kw...)
end
export optimize_limit


"""
    optimize_design!(CH::ControlChart, rlsim_oc::Function, settings::OptSettings=OptSettings(CH); optimizer = :LN_BOBYQA, solver = :Bootstrap, hmax::Float64 = 20.0, kw...)

Optimizes the design of a control chart `CH` using a specified optimization algorithm.

### Arguments
- `CH::ControlChart`: The control chart object to optimize.
- `rlsim_oc::Function`: A function that simulates the out-of-control state of the control chart.
- `settings::OptSettings`: The optimization settings that control the optimization routine (default: `OptSettings(CH)`).

### Keyword Arguments
- `optimizer::Symbol`: The optimization algorithm to use (default: `:LN_BOBYQA`).
- `solver::Symbol`: The root-finding algorithm to use for control limit estimation (default: `:Bootstrap`).
- `hmax::Float64`: The maximum value of the control limit, used only when the solver is set to `:Bisection` (default: 100.0)
- `kw...`: Additional keyword arguments to pass to the solver algorithm.

# Returns
The optimized design parameters of the control chart.
"""
function optimize_design!(CH::ControlChart, rlsim_oc::Function, settings::OptSettings=OptSettings(CH); optimizer::Symbol = :Grid, solver::Symbol = :Bootstrap, hmax::Float64 = 100.0, kw...)
    CH_ = deepcopy(CH)

    @unpack nsims = settings

    valid_optimizers = [:Grid, :LN_BOBYQA, :LN_COBYLA, :LN_SBPLX, :LN_NELDERMEAD, :LN_PRAXIS, :LN_NEWOA]
    @assert optimizer in valid_optimizers "Unknown optimization algorithm. Valid optimizers are: $(valid_optimizers)"

    function rlconstr(par::Vector, grad::Vector)::Float64
        set_design!(CH_, par)
        optimize_limit!(CH_, solver; hmax=hmax, kw...)
        if settings.verbose
            print("$(round.(par, digits=6))\t")
        end
        ret = SPM.measure([rlsim_oc(CH_) for _ in 1:nsims], CH_, verbose=settings.verbose)
        return ret
    end

    if optimizer == :Grid
        set_design!(CH, optimize_grid(CH, rlconstr, settings))
    else
        set_design!(CH, optimize_nlopt(CH, rlconstr, settings, optimizer=optimizer))
    end
    set_h!(get_limit(CH), get_h(get_limit(CH_)))
    return get_design(CH)
end
export optimize_design!

"""
    optimize_design(CH::ControlChart, rlsim_oc::Function, settings::OptSettings=OptSettings(CH); optimizer::Symbol = :LN_BOBYQA, solver::Symbol = :SACL, nsims_opt::Int = 1000, kw...)

Optimizes the design of a control chart using a specified optimization algorithm.

### Arguments
- `CH::ControlChart`: The control chart object to optimize.
- `rlsim_oc::Function`: A function that simulates the out-of-control state of the control chart.
- `settings::OptSettings`: The optimization settings that control the optimization routine (default: `OptSettings(CH)`).

### Keyword Arguments
- `optimizer::Symbol`: The optimization algorithm to use (default: `:LN_BOBYQA`).
- `solver::Symbol`: The root-finding algorithm to use for control limit estimation (default: `:Bootstrap`).
- `hmax::Float64`: The maximum value of the control limit, used only for the bisection algorithm (default: 100.0)
- `kw...`: Additional keyword arguments to pass to the solver algorithm.

# Returns
The optimized design parameters of the control chart.
"""
function optimize_design(CH::ControlChart, rlsim_oc::Function, settings::OptSettings=OptSettings(CH); kw...)
    CH_ = deepcopy(CH)
    optimize_design!(CH_, rlsim_oc, settings; kw...)
end
export optimize_design