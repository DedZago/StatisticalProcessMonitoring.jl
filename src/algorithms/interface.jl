"""
    optimize_limit(CH::ControlChart, solver::Symbol = :SA; kw...)

Optimizes the control limit of a ControlChart object.

### Args
- `CH::ControlChart`: The ControlChart object to optimize.
- `solver::Symbol`: The solver algorithm to use. Defaults to `:SA`.
- `kw...`: Additional keyword arguments to pass to the algorithm.

### Returns
    The optimized control limit value.

### Raises
    ValueError: If the optimization method specified in settings is unknown.

### Example
    optimize_limit!(my_chart, settings=OptSettings(ic_solver=:SA))
"""
function optimize_limit!(CH::ControlChart, solver::Symbol=:SA; kw...)
    valid_solvers = [:SA, :Bisection, :Combined, :Bootstrap]
    @assert solver in valid_solvers "Unknown control limit solver. Valid solvers are $(valid_solvers)"

    #FIXME: better assertion for solver
    if solver == :SA
        return saCL!(CH; kw...)
    elseif solver == :Bisection
        return bisectionCL!(CH; kw...)
    elseif solver == :Combined
        return combinedCL!(CH; kw...)
    elseif solver == :Bootstrap
        error("To be implemented yet.")
    end
end
export optimize_limit!

"""
    optimize_limit(CH::ControlChart, solver::Symbol = :SA; kw...)

Optimizes the control limit of a ControlChart object, without modifying the original ControlChart object.

### Args
- `CH::ControlChart`: The ControlChart object to optimize.
- `solver::Symbol`: The solver algorithm to use. Defaults to `:SA`.
- `kw...`: Additional keyword arguments to pass to the algorithm.

### Returns
    The optimized control limit value.

### Raises
    ValueError: If the optimization method specified in settings is unknown.

### Example
    optimize_limit(my_chart, settings=OptSettings(ic_solver=:SA))
"""
function optimize_limit(CH::ControlChart, solver::Symbol = :SA; kw...)
    CH_ = deepcopy(CH)
    optimize_limit!(CH_, solver; kw...)
end
export optimize_limit


"""
    optimize_design!(CH::ControlChart, rlsim_oc::Function, settings::OptSettings=OptSettings(CH); optimizer::Symbol = :LN_BOBYQA, solver::Symbol = :SACL, nsims_opt::Int = 1000, trace::Int, kw...)

Optimizes the design of a control chart using a specified optimization algorithm.

# Arguments
- `CH::ControlChart`: The control chart object to optimize.
- `rlsim_oc::Function`: A function that simulates the out-of-control state of the control chart.
- `settings::OptSettings`: The optimization settings to use. Defaults to `OptSettings(CH)`.
- `optimizer::Symbol`: The optimization algorithm to use. Defaults to `:LN_BOBYQA`.
- `solver::Symbol`: The solver algorithm to use. Defaults to `:SACL`.
- `trace::Int`: The trace level for printing optimization progress. Set to `0` to disable printing.
- `kw...`: Additional keyword arguments to pass to the solver algorithm.

# Returns
The optimized design of the control chart.
"""
function optimize_design!(CH::ControlChart, rlsim_oc::Function, settings::OptSettings=OptSettings(CH); optimizer::Symbol = :LN_BOBYQA, solver::Symbol = :SA, trace::Int, kw...)
    CH_ = deepcopy(CH)

    @unpack nsims = settings

    valid_optimizers = [:Grid, :SPSA, :LN_BOBYQA, :LN_COBYLA, :LN_SBPLX, :LN_NELDERMEAD, :LN_PRAXIS, :LN_NEWOA]
    @assert optimizer in valid_optimizers "Unknown tuning parameter optimizer. Valid optimizers are $(valid_optimizers)"

    function rlconstr(par::Vector, grad::Vector)::Float64
        set_design!(CH_, par)
        optimize_limit!(CH_, solver, kw...)
        if trace > 0
            print("$(round.(par, digits=6))\t")
        end
        ret = SPM.measure([rlsim_oc(CH_) for _ in 1:nsims], CH_, verbose=trace > 0)
        return ret
    end

    if optimizer == :Grid
        set_design!(CH, optimize_grid(CH, rlconstr, settings))
    elseif optimizer == :SPSA
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
    optimize_design(CH::ControlChart, rlsim_oc::Function, settings::OptSettings=OptSettings(CH); optimizer::Symbol = :LN_BOBYQA, solver::Symbol = :SACL, nsims_opt::Int = 1000, trace::Int, kw...)

Optimizes the design of a control chart using a specified optimization algorithm.

# Arguments
- `CH::ControlChart`: The control chart object to optimize.
- `rlsim_oc::Function`: A function that simulates the out-of-control state of the control chart.
- `settings::OptSettings`: The optimization settings to use. Defaults to `OptSettings(CH)`.
- `optimizer::Symbol`: The optimization algorithm to use. Defaults to `:LN_BOBYQA`.
- `solver::Symbol`: The solver algorithm to use. Defaults to `:SACL`.
- `trace::Int`: The trace level for printing optimization progress. Set to `0` to disable printing.
- `kw...`: Additional keyword arguments to pass to the solver algorithm.

# Returns
The optimized design of the control chart.
"""
function optimize_design(CH::ControlChart, rlsim_oc::Function, settings::OptSettings=OptSettings(CH); optimizer::Symbol = :LN_BOBYQA, solver::Symbol = :SA, trace::Int, kw...)
    CH_ = deepcopy(CH)
    optimize_design!(CH_, rlsim_oc, settings; optimizer=optimizer, solver=solver, trace=trace, kw...)
end
export optimize_design