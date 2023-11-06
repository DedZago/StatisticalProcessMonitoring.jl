using Parameters 

@with_kw mutable struct OptSettings{F, I, B} 
    # Global parameter optimization options
    x_tol::F = 1e-03                            # Tolerance for the parameter
    nsims::I = 1000                             # Number of simulations for approximating the OC RL
    maxiter::I = 1000                           # Maximum number of iterations in the algorithm
    verbose::B = true                           # Verbosity of the optimization algorithm
    minpar::Vector{F} = [-Inf]                  # Minimum value(s) of the tuning parameter(s)
    maxpar::Vector{F} = [Inf]                   # Maximum value(s) of the tuning parameter(s)

    # Grid settings
    m_grid::I = 10                              # Number of segments in the grid search algorithm

    # SPSA settings
    initial_step_size_spsa::F = 0.05
    expected_n_loss_eval_spsa::I = 100
    n_adaptive_gain_spsa::I = 20
    gamma_spsa::F = 0.02

    @assert x_tol > 0
    @assert nsims > 0
    @assert maxiter > 0

    @assert initial_step_size_spsa > 0
    @assert expected_n_loss_eval_spsa > 0
    @assert n_adaptive_gain_spsa > 0
    @assert gamma_spsa > 0
end

function OptSettings(CH::ControlChart; kw...)
    println(kw)
    kw_dict = Dict(kw)
    kw_dict = convert(Dict{Symbol, Any}, kw_dict)
    println(kw_dict)
    d = length(get_design(CH))
    if haskey(kw_dict, :minpar)
        @assert length(kw_dict[:minpar]) == d
    else
        kw_dict[:minpar] = [-Inf for _ in 1:d]
    end
    if haskey(kw_dict, :maxpar)
        @assert length(kw_dict[:maxpar]) == d
    else
        kw_dict[:maxpar] = [Inf for _ in 1:d]
    end
    return OptSettings(; kw_dict...)
end
export OptSettings