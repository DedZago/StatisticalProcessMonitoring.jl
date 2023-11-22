"""
    run_sim(CH::AbstractChart)

Simulates a run length for the control chart `CH` by sampling new data from its Phase II object.

### Inputs
* `CH` - A control chart.

### Returns
* An `Int` representing the simulated run length.
"""
function run_sim(CH::AbstractChart; maxiter::Real = Inf)
    CH_ = shallow_copy_sim(CH)
    maxrl = get_maxrl(CH_)
    i = 0.0
    while i < maxiter && i < maxrl
        i = i + 1
        update_chart!(CH_, new_data!(CH_))
        is_IC(CH_) || break
    end
    return i
end
export run_sim

"""
    run_sim(CH::AbstractChart, DGP::AbstractPhase2)

Simulates a run length for the control chart `CH` by sampling new data from the provided data-generating process `DGP`.

### Inputs
* `CH` - A control chart.
* `DGP` - An AbstractPhase2 object.

### Returns
* An `Int` representing the simulated run length.
"""
function run_sim(CH::AbstractChart, DGP::AbstractPhase2)
    CH_ = shallow_copy_sim(CH)
    maxrl = get_maxrl(CH_)
    i = 0.0
    while i < maxrl
        i = i + 1
        update_chart!(CH_, new_data!(DGP))
        is_IC(CH_) || break
    end
    return i
end

"""
    run_sim_sa(CH::AbstractChart, maxiter::Real, delta::Real)

Simulates a run length for the control chart `CH` by sampling new data from the Phase II object, to be used by the stochastic approximation algorithm implemented in the `saCL!` function.

### Inputs
* `CH` - A control chart.
* `maxiter` - The maximum value of the run length.
* `delta` - A value controlling how much the control limit must be shifted for the gain estimation during the first stage.

### Returns
* A `NamedTuple` containing the simulated run length, `rl`, the simulated run length with control limit shifted by `delta`, `rlPlus`, and the simulated run length with control limit shifted by `-delta`, `rlMinus`.
"""
function run_sim_sa(CH::AbstractChart; maxiter::Real = Inf, delta::Real = 0.0)
    @assert delta >= 0.0
    CH_ = shallow_copy_sim(CH)
    maxrl = min(get_maxrl(CH_), maxiter)
    rl = rlPlus = rlMinus = round(maxrl)
    notDoneP = notDoneM = (delta > 0.0)
    notDoneRl = true
    notDone = notDoneRl + notDoneM + notDoneP
    i = 0.0
    h = deepcopy(get_h(get_limit(CH)))
    while i < maxrl && (notDone > 0)
        i = i+1
        update_chart!(CH_, new_data!(CH_))
        if notDoneRl
            set_limit!(CH_, h)
            if is_OC(CH_)
                notDoneRl = false
                rl = i
            end
        end
        if notDoneP
            set_limit!(CH_, h + delta)
            if is_OC(CH_)
                notDoneP = false
                rlPlus = i
            end
        end
        if notDoneM
            set_limit!(CH_, h - delta)
            if is_OC(CH_)
                notDoneM = false
                rlMinus = i
            end
        end
        notDone = notDoneRl + notDoneM + notDoneP
    end
    set_limit!(CH_, h)
    if delta == 0.0
        rlPlus = rl
        rlMinus = rl
    end
    return (rl = rl, rlPlus = rlPlus, rlMinus = rlMinus)
end


function run_sim_sa(CH::MultipleControlChart; maxiter::Real = Inf, delta::Real = 0.0)
    @assert delta >= 0.0
    CH_ = shallow_copy_sim(CH)
    maxrl = min(get_maxrl(CH_), maxiter)
    nstat = length(get_statistic(CH_))
    rl = fill(round(maxrl), nstat)
    rlPlus = fill(round(maxrl), nstat)
    rlMinus = fill(round(maxrl), nstat)
    notDoneP = [deepcopy(delta > 0.0) for _ in 1:nstat]
    notDoneM = [deepcopy(delta > 0.0) for _ in 1:nstat]
    notDoneRl = fill(true, nstat)
    notDone = notDoneRl + notDoneM + notDoneP
    i = 0
    h = deepcopy(get_h(get_limit(CH_)))
    OC_vector = BitVector(undef, nstat)
    while i < maxrl && any(notDone .> 0)
        i = i+1
        update_chart!(CH_, new_data!(CH_))
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
                set_limit!(CH_, h .+ delta)
                OC_vector[:] = is_OC_vec(CH_)
                if OC_vector[j]
                    notDoneP[j] = false
                    rlPlus[j] = i
                end
            end
            if notDoneM[j]
                set_limit!(CH_, h .- delta)
                OC_vector[:] = is_OC_vec(CH_)
                if OC_vector[j]
                    notDoneM[j] = false
                    rlMinus[j] = i
                end
            end
            notDone[:] = notDoneRl + notDoneM + notDoneP
            # @show notDone
        end
    end
    set_limit!(CH_, h)
    return (rl = rl, rlPlus = rlPlus, rlMinus = rlMinus)
end
export run_sim_sa

"""
    run_sim_oc(CH::AbstractChart; shift = 0.0)

Simulates a run length under location shift for the control chart `CH` by sampling new data from its Phase II object.

### Inputs
* `CH::AbstractChart`:  A control chart.
* `shift::Float64` - The magnitude of location shift.

### Returns
* An `Int` representing the simulated run length.
"""
function run_sim_oc(CH::AbstractChart; shift = 0.0, maxiter::Real = Inf)
    CH_ = shallow_copy_sim(CH)
    maxrl = get_maxrl(CH_)
    i = 0.0
    while i < maxiter && i < maxrl
        i = i + 1
        update_chart!(CH_, new_data!(CH_) .+ shift)
        is_IC(CH_) || break
    end
    return i
end
export run_sim_oc


"""
    run_path_sim(CH::AbstractChart; maxiter)
    run_path_sim(CH::MultipleControlChart; maxiter)

Simulates a run length path for the control chart `CH` by sampling new data from its Phase II object.

### Inputs
* `CH::AbstractChart` - A control chart.
* `maxiter::Real` - The maximum value of the run length. Defaults to `min(maxrl(CH), 10*get_nominal_value(CH))`

### Returns
* A vector containing the simulated values of the control chart.
"""
function run_path_sim(CH::AbstractChart; maxiter::Real = min(get_maxrl(CH), 10*get_nominal_value(CH)))
    CH_ = shallow_copy_sim(CH)
    i = 0
    out = zeros(Int(maxiter))
    while i < maxiter
        i = i + 1
        update_chart!(CH_, new_data!(CH_))
        out[i] = get_value(CH_)
    end
    return out
end

function run_path_sim(CH::MultipleControlChart; maxiter::Real = min(get_maxrl(CH), 10*get_nominal_value(CH)))
    CH_ = shallow_copy_sim(CH)
    i = 0
    out = zeros(Int(maxiter), length(get_statistic(CH)))
    while i < maxiter
        i = i + 1
        update_chart!(CH_, new_data!(CH_))
        out[i, :] = get_value(CH_)
    end
    return out
end
export run_path_sim