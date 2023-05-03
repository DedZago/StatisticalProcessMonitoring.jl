"""
	run_sim(CH::AbstractChart)

Simulates a run length for the control chart `CH` by sampling new data from the Phase I object.
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

function run_sim_sa(CH::AbstractChart, maxiter, deltaSA)
    CH_ = shallow_copy_sim(CH)
    maxrl = min(get_maxrl(CH_), maxiter)
    rl = rlPlus = rlMinus = Int(round(maxrl))
    notDoneP = notDoneM = (deltaSA > 0.0)
    notDoneRl = true
    notDone = notDoneRl + notDoneM + notDoneP
    i = 0
    h = deepcopy(get_limit_value(CH))
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