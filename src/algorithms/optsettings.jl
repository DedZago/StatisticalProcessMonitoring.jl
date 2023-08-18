using Parameters 

@with_kw mutable struct OptSettings{F, I, B} 
    # Global parameter optimization options
    x_tol::F = 1e-05
    nsims::I = 1000
    maxiter::I = 1000
    verbose::B = true
    minpar::Vector{F} = [-Inf]
    maxpar::Vector{F} = [Inf]

    # Grid settings
    m_grid::I = 10

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
        @assert length(minpar) == d
    else
        kw_dict[:minpar] = [-Inf for _ in 1:d]
    end
    if haskey(kw_dict, :maxpar)
        @assert length(minpar) == d
    else
        kw_dict[:maxpar] = [Inf for _ in 1:d]
    end
    return OptSettings(; kw_dict...)
end
export OptSettings