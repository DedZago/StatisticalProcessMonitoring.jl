using Parameters 

@with_kw mutable struct OptSettings{F, I, S, B, F1} 
    # Global options
    trace::I = 0
    ic_solver::S = :SA
    rlsim::F1 = ifelse(ic_solver == :Bisection, run_sim, run_sim_sa)

    # saCL options
    hmin_sa::F = sqrt(eps())
    Nfixed_sa::I = ifelse(ic_solver==:Combined, 100, 500)
    Afixed_sa::F = 0.1
    Amin_sa::F = 0.1
    Amax_sa::F = 100.0  
    delta_sa::F = 0.1
    q_sa::F = 0.55
    gamma_sa::F = 0.02
    Nmin_sa::I = ifelse(ic_solver==:Combined, 100, 1000)
    z_sa::F = 3.0
    Cmrl_sa::F = 10.0
    maxiter_sa::I = ifelse(ic_solver==:Combined, 100, 100_000)
    verbose_sa::B = trace > 1

    # Bisection options
    hmin_bi::F = sqrt(eps())
    hmax_bi::F = 50.0
    maxiter_bi::I = 50
    nsims_bi::I = 10_000
    trunc_bi::I = 10_000
    x_tol_bi::F = 1e-06
    f_tol_bi::F = 1.0
    verbose_bi::B = trace > 1
    inflate_bi::F = 1.05

    # Global parameter optimization options
    x_tol_opt::F = 1e-05
    nsims_opt::I = 10_000
    maxiter_opt::I = 100
    method_opt::S = :LN_BOBYQA
    verbose_opt::I = trace
    minpar_opt::Vector{F} = [-Inf]
    maxpar_opt::Vector{F} = [Inf]

    # Grid settings
    m_grid::I = 10

    @assert trace >= 0
    @assert ic_solver in [:SA, :Bisection, :Combined]
    @assert hmin_sa > 0
    @assert Nfixed_sa > 0
    @assert Afixed_sa > 0
    @assert Amin_sa > 0
    @assert Amax_sa > Amin_sa
    @assert delta_sa > 0
end
export OptSettings