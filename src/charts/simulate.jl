"""
    run_sim(CH::AbstractChart)

Simulates a run length for the control chart `CH` by sampling new data from the Phase I object.

### Inputs
* `CH` - A control chart.

### Returns
* An `Int`.
"""
function run_sim(CH::AbstractChart)
    CH_ = shallow_copy_sim(CH)
    maxrl = get_maxrl(CH_)
    i = 0
    while is_IC(CH_) && i < maxrl
        i = i + 1
        update_chart!(CH_, new_data(CH_))
    end
    return i
end
export run_sim

"""
    run_sim_sa(CH::AbstractChart, maxiter::Real, deltaSA::Real)

Simulates a run length for the control chart `CH` by sampling new data from the Phase I object, to be used by the stochastic approximation algorithm implemented in the `saCL!` function.

### Inputs
* `CH` - A control chart.
* `maxiter` - The maximum value of the run length.
* `deltaSA` - A value controlling how much the control limit must be shifted for the gain estimation during the first stage.

### Returns
* A `NamedTuple` containing the simulated run length, `rl`, the simulated run length with control limit shifted by `deltaSA`, `rlPlus`, and the simulated run length with control limit shifted by `-deltaSA`, `rlMinus`.
"""
function run_sim_sa(CH::AbstractChart, maxiter::Real, deltaSA::Real)
    @assert deltaSA >= 0.0
    CH_ = shallow_copy_sim(CH)
    maxrl = min(get_maxrl(CH_), maxiter)
    rl = rlPlus = rlMinus = Int(round(maxrl))
    notDoneP = notDoneM = (deltaSA > 0.0)
    notDoneRl = true
    notDone = notDoneRl + notDoneM + notDoneP
    i = 0
    h = deepcopy(get_value(get_limit(CH)))
    while i < maxrl && (notDone > 0)
        i = i+1
        update_chart!(CH_, new_data(CH_))
        if notDoneRl
            set_limit!(CH_, h)
            if is_OC(CH_)
                notDoneRl = false
                rl = i
            end
        end
        if notDoneP
            set_limit!(CH_, h .+ deltaSA)
            if is_OC(CH_)
                notDoneP = false
                rlPlus = i
            end
        end
        if notDoneM
            set_limit!(CH_, h .- deltaSA)
            if is_OC(CH_)
                notDoneM = false
                rlMinus = i
            end
        end
        notDone = notDoneRl + notDoneM + notDoneP
    end
    set_limit!(CH_, h)
    return (rl = rl, rlPlus = rlPlus, rlMinus = rlMinus)
end
export run_sim_sa