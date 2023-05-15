#FIXME: test
function optimize_limit!(CH::ControlChart; rlsim::Function = run_sim_sa, settings = OptSettings())
    ic_solver = settings.ic_solver
    if ic_solver == :SA
        return saCL!(CH, rlsim=rlsim, settings = settings)
    elseif ic_solver == :Bisection
        return bisectionCL!(CH, rlsim=rlsim, settings = settings)
    elseif ic_solver == :Combined
        combinedCL!(CH, rlsim=rlsim, settings = settings)
    else
        error("Unknown optimization method.")
    end
end
export optimize_limit!

#FIXME: test
function optimize_limit(CH::ControlChart; rlsim::Function = run_sim_sa, settings = OptSettings())
    CH_ = deepcopy(CH)
    optimize_limit!(CH_, rlsim=rlsim, settings=settings)
end
export optimize_limit