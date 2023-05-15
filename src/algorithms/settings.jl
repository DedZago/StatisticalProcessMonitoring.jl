using Parameters 

@with_kw struct OptimizationSettings{F, I, S, B} 
    # Global options
    trace::I = 0
    ic_solver::S = :SA

    # saCL options
    hmin_sa::F = sqrt(eps())
    Nfixed_sa::I = 500
    Afixed_sa::F = 0.1
    Amin_sa::F = 0.1
    Amax_sa::F = 100.0  
    delta_sa::F = 0.1
    q_sa::F = 0.55
    gamma_sa::F = 0.02
    Nmin_sa::I = 1000
    z_sa::F = 3.0
    Cmrl_sa::F = 10.0
    maxiter_sa::I = 100_000
    verbose_sa::B = trace > 0

    # Bisection options
    hmin_bi::F = sqrt(eps())
    hmax_bi::F = 50.0
    maxiter_bi::I = 50
    nsims_bi::I = 10_000
    trunc_bi::I = 10_000
    x_tol_bi::F = 1e-06
    f_tol_bi::F = 1.0
    verbose_bi::B = trace > 0
    inflate_bi::F = 1.03

    # Global parameter optimization options
    x_tol_opt::F = 1e-05
    nsims_opt::I = 10_000
    maxiter_opt::I = 100
    method_opt::S = :LN_BOBYQA
    verbose_opt::I = trace
end
export OptimizationSettings