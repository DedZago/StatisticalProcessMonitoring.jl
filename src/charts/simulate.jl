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
    while i < maxrl
        i = i + 1
        update_chart!(CH_, new_data(CH_))
        is_IC(CH_) || break
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
    rl = rlPlus = rlMinus = round(maxrl)
    notDoneP = notDoneM = (deltaSA > 0.0)
    notDoneRl = true
    notDone = notDoneRl + notDoneM + notDoneP
    i = 0
    h = deepcopy(get_h(get_limit(CH)))
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
            set_limit!(CH_, h + deltaSA)
            if is_OC(CH_)
                notDoneP = false
                rlPlus = i
            end
        end
        if notDoneM
            set_limit!(CH_, h - deltaSA)
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


function run_sim_sa(CH::MultipleControlChart, maxiter::Real, deltaSA::Real)
    @assert deltaSA >= 0.0
    CH_ = shallow_copy_sim(CH)
    maxrl = min(get_maxrl(CH_), maxiter)
    nstat = length(get_statistic(CH_))
    rl = fill(round(maxrl), nstat)
    rlPlus = fill(round(maxrl), nstat)
    rlMinus = fill(round(maxrl), nstat)
    notDoneP = [deepcopy(deltaSA > 0.0) for _ in 1:nstat]
    notDoneM = [deepcopy(deltaSA > 0.0) for _ in 1:nstat]
    notDoneRl = fill(true, nstat)
    notDone = notDoneRl + notDoneM + notDoneP
    i = 0
    h = deepcopy(get_value(get_limit(CH_)))
    OC_vector = BitVector(undef, nstat)
    while i < maxrl && any(notDone .> 0)
        i = i+1
        update_chart!(CH_, new_data(CH_))
        # @show get_value(CH_)
        for j in 1:nstat
            if notDoneRl[j]
                set_limit!(CH_, h)
                OC_vector[:] = is_OC_vec(CH_)
                if OC_vector[j]
                    notDoneRl[j] = false
                    rl[j] = i
                end
            end
            if notDoneP[j]
                set_limit!(CH_, h .+ deltaSA)
                OC_vector[:] = is_OC_vec(CH_)
                if OC_vector[j]
                    notDoneP[j] = false
                    rlPlus[j] = i
                end
            end
            if notDoneM[j]
                set_limit!(CH_, h .- deltaSA)
                OC_vector[:] = is_OC_vec(CH_)
                if OC_vector[j]
                    notDoneM[j] = false
                    rlMinus[j] = i
                end
            end
            notDone[:] = notDoneRl + notDoneM + notDoneP
        end
    end
    set_limit!(CH_, h)
    return (rl = rl, rlPlus = rlPlus, rlMinus = rlMinus)
end
export run_sim_sa