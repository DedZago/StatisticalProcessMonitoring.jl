#FIXME: test
function optimize_limit!(CH::ControlChart; settings = OptSettings())
    @unpack ic_solver, rlsim = settings
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
function optimize_limit(CH::ControlChart; settings = OptSettings())
    CH_ = deepcopy(CH)
    optimize_limit!(CH_, settings=settings)
end
export optimize_limit