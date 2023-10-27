"""
    optimize_grid(CH::ControlChart, rlconstr::Function, settings::OptSettings)
Optimizes a control chart by finding the best set of parameters using a grid search.

### Arguments
- `CH::ControlChart`: The control chart object whose parameters must be optimized.
- `rlconstr::Functiom`: The function that evaluates the OC performance of the control chart.
- `settings::OptSettings`: The optimization settings.

### Returns
- par_current (Vector{Float64}): the optimal set of parameters found by the optimization algorithm.

### References
Qiu, P. (2008). Distribution-Free Multivariate Process Control Based on Log-Linear Modeling. IIE Transactions, 40(7), 664-677. https://doi.org/10.1080/07408170701744843
"""
function optimize_grid(CH::ControlChart, rlconstr::Function, settings::OptSettings)
    @unpack minpar, maxpar, maxiter, nsims, m_grid, x_tol = settings

    p = length(minpar)
    grid_vec = Vector{Vector{Float64}}(undef, p)
    step = Vector{Float64}(undef, p)
    par_current = zeros(p)
    par_old = zeros(p)
    par = zeros(p)
    i = 0
    cnt = 0
    best_RL = 0.0
    RL = 0.0
    min_par = deepcopy(minpar)
    max_par = deepcopy(maxpar)
    while i < maxiter
        i = i+1
        par_old = deepcopy(par_current)
        step = [max_par[w] - min_par[w] for w in 1:p]
        for w in 1:p
            grid_vec[w] = [min_par[w] + step[w] * j/m_grid for j in 0:m_grid]
        end
        parameter_grid = Base.product(grid_vec...)
        cnt = 0
        for par_tuple in parameter_grid
            cnt += 1
            par = collect(par_tuple)
            RL = rlconstr(par, [0.0])
            if (cnt == 1) | (RL < best_RL)
                best_RL = RL       
                par_current = par
            end
        end
        if sum(abs.(par_current - par_old)) < x_tol
            @show par_current, par_old
            break
        end
        for w in 1:p
            if par_current[w] == max_par[w]
                par_current[w] -= step[w]/m_grid
            elseif par_current[w] == min_par[w]
                par_current[w] += step[w]/m_grid
            end
            max_par[w] = par_current[w] + step[w]/m_grid
            min_par[w] = par_current[w] - step[w]/m_grid
        end
    end
    return par_current
end
export optimize_grid
