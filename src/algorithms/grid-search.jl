"""
    optimize_grid(CH::ControlChart, rlconstr::Function, settings::OptSettings)

#FIXME: docstring
#FIXME: test
"""
function optimize_grid(CH::ControlChart, rlconstr::Function, settings::OptSettings)
    @unpack minpar_opt, maxpar_opt, maxiter_opt, nsims_opt, m_grid, x_tol_opt = settings

    p = length(minpar_opt)
    grid_vec = Vector{Vector{Float64}}(undef, p)
    step = Vector{Float64}(undef, p)
    par_current = zeros(p)
    par_old = zeros(p)
    par = zeros(p)
    i = 0
    cnt = 0
    best_RL = 0.0
    RL = 0.0
    min_par = deepcopy(minpar_opt)
    max_par = deepcopy(maxpar_opt)
    while i < maxiter_opt
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
        if sum(abs.(par_current - par_old)) < x_tol_opt
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
