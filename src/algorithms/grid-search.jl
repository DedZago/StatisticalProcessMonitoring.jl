#FIXME: grid search algorithm to optimize the parameters of the control chart
function optimize_parameter!(CH; rlsim::Function = run_sim_sa, rlsim_oc::Function, nsims::Real = 10000, method = :BOBYQA, ic_solver::Symbol = :SA)
    CH_ = shallow_copy_sim(CH)

    function rlconstr(par::Vector)
        set_parameter!(CH_, par)
        if ic_solver == :SA
            saCL!(CH_, rlsim=rlsim)
        elseif ic_solver == :Combined
            combinedCL!(CH_, rlsim=rlsim)
        end
        return SPM.measure([rlsim_oc(CH_) for _ in 1:nsims], CH_)
    end

    if method == :BOBYQA
        set_parameter!(CH, optimize_bobyqa(rlconstr))
    elseif method == :Grid
        set_parameter!(CH, optimize_grid(rlconstr))
    end
end


function optimize_grid(rlconst::Function)
    error("Not implemented yet")
end

function optimize_bobyqa(rlconst::Function)
    error("Not implemented yet")
end