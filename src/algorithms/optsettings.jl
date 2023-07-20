using Parameters 

@with_kw mutable struct OptSettings{F, I, S, B, F1} 
    # Global options
    trace::I = 0
    ic_solver::S = :SA
    rlsim::F1 = ifelse(ic_solver == :Double, run_path_sim, ifelse(ic_solver == :Bisection, run_sim, run_sim_sa))
    method_opt::S = :LN_BOBYQA

    # saCL options
    hmin_sa::F = sqrt(eps())
    Nfixed_sa::I = ifelse(ic_solver==:Combined || method_opt == :SPSA, 100, 500)
    Afixed_sa::F = 0.1
    Amin_sa::F = 0.1
    Amax_sa::F = 100.0  
    delta_sa::F = 0.1
    q_sa::F = 0.55
    gamma_sa::F = 0.02
    Nmin_sa::I = ifelse(ic_solver==:Combined || method_opt == :SPSA, 100, 1000)
    z_sa::F = 3.0
    Cmrl_sa::F = 10.0
    maxiter_sa::I = ifelse(ic_solver==:Combined || method_opt == :SPSA, 100, 100_000)
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

    # CombinedCL
    inflate_bi::F = 1.05
    
    # Double bootstrap bisection

    # Global parameter optimization options
    x_tol_opt::F = 1e-05
    nsims_opt::I = ifelse(method_opt == :SPSA, 100, 10_000)
    maxiter_opt::I = 1000
    verbose_opt::I = trace
    minpar_opt::Vector{F} = [-Inf]
    maxpar_opt::Vector{F} = [Inf]

    # Grid settings
    m_grid::I = 10

    # SPSA settings
    initial_step_size_spsa::F = 0.05
    expected_n_loss_eval_spsa::I = 100
    n_adaptive_gain_spsa::I = 20
    gamma_spsa::F = 0.02

    @assert trace >= 0
    @assert ic_solver in [:SA, :Bisection, :Combined, :Double]
    @assert hmin_sa > 0
    @assert Nfixed_sa > 0
    @assert Afixed_sa > 0
    @assert Amin_sa > 0
    @assert Amax_sa > Amin_sa
    @assert delta_sa > 0

    @assert hmin_bi > 0
    @assert hmax_bi > hmin_bi
    @assert maxiter_bi > 0
    @assert nsims_bi > 0
    @assert trunc_bi > 0
    @assert x_tol_bi > 0
    @assert f_tol_bi > 0

    @assert x_tol_opt > 0
    @assert nsims_opt > 0
    @assert maxiter_opt > 0

    @assert initial_step_size_spsa > 0
    @assert expected_n_loss_eval_spsa > 0
    @assert n_adaptive_gain_spsa > 0
    @assert gamma_spsa > 0
end
export OptSettings