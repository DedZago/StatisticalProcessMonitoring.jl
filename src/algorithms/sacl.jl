mean(x) = sum(x) / length(x)

"""
    calculate_limit_gradient(CH::AbstractChart, rl::Real)

Calculate the gradient for the optimization of the control limit, according to
Equation (9) of Capizzi, G., and Masarotto, G. (2016). "Efficient Control Chart
Calibration by Simulated Stochastic Approximation". IIE Transactions 48 (1).
https://doi.org/10.1080/0740817X.2015.1055392.

"""
function calculate_limit_gradient(CH::AbstractChart, rl::Real)
   return calculate_limit_gradient(get_nominal(CH), rl)
end
export calculate_limit_gradient

function calculate_limit_gradient(nominal::ARL, rl)
    #TODO: extend to multivariate
    return (minimum(rl) - get_value(nominal) .+ rl .- mean(rl)) ./ get_value(nominal)
end

function calculate_limit_gradient(nominal::QRL, rl)
    #TODO: check gradient correctness
    #TODO: extend to multivariate
    return -(float(rl <= get_value(nominal)) - nominal.qtl)
end


"""
    saCL(CH::ControlChart, rlsim::Function, kw...)

Computes the control limit for a control chart such that it satisfies E[RL] = Arl0.

Keyword arguments are:
"""
function saCL(CH::ControlChart; rlsim::Function = run_sim_sa, Nfixed::Int=500, Afixed::Real=0.1, Amin::Real=0.1, Amax::Real=100.0, deltaSA::Real=0.1, q::Real=0.55, gamma::Real=0.02, Nmin::Int=1000, z::Real = 3.0, Cmrl::Real=10.0, maxiter::Real = 4e05, verbose::Bool=true)

    #TODO: Find out why is it not working
    v = (z/gamma)^2
    eps = 1e-06
    h = deepcopy(get_limit_value(CH))
    score = zero(h)
    s2    = zero(h)
    D     = zero(h)

    if verbose println("Running SA ...") end

    # Stage 1 - adaptive gain estimation
    D          = zero(h)
    rl         = zero(h)
    rlPlus     = zero(h)
    rlMinus    = zero(h)
    scorePlus  = zero(h)
    scoreMinus = zero(h)
    Arl0 = get_nominal_value(CH)
    if verbose println("Running adaptive gain ...") end
    for i in 1:Nfixed
        set_limit!(CH, h)
        rl, rlPlus, rlMinus = rlsim(CH, Cmrl * Arl0 * sqrt(i + Nfixed), deltaSA)
        calculate_limit_gradient(CH, rl)
        h = max.(eps, h .- Afixed * score ./ (i^q))
        calculate_limit_gradient(CH, rlPlus)
        calculate_limit_gradient(CH, rlMinus)
        D .+= (scorePlus .- scoreMinus) ./ i
    end
    @. D = 1.0 / (max(1.0/Amax, min(1.0/Amin, D/(2.0*deltaSA))))
    if verbose println("Estimated gain D = $(D)") end


    # Stage 2 - Stochastic approximations
    hm = zero(h)
    i = 0
    if verbose println("Running optimization ...") end
    while i < (maxiter + Nmin)
        if verbose && (i % floor(maxiter / 20) == 0)
            println("i: $(i)/$(Int(trunc(maxiter)) + Nmin)\th: $(h)\thm: $(hm)")
        end
        i += 1
        set_limit!(CH, h)
        rl, _, _ = rlsim(CH, Cmrl * Arl0 * sqrt(i + Nfixed), 0.0)
        score = calculate_limit_gradient(CH, rl)
        h = max.(eps, h .- D * score ./ (i^q))
        if i > Nmin
            Ndenom = (i - Nmin)
            hm = hm .+ (h .- hm) ./ Ndenom 
            s2 = s2 .+ (score .* score .- s2) ./ Ndenom
        end
        if (i > Nmin) && (i > v * maximum(s2))
            if verbose println("Convergence!") end
            break
        end
    end
    
    if verbose println("Done.") end

    return (h=hm, iter=i)
end
export saCL