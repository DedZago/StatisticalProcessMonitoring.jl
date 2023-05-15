using SPM
using Test

# Unit tests for OptSettings struct
@testset "OptSettings struct" begin
    # Global options
    @testset "Global options" begin
        settings = OptSettings(trace = 2, ic_solver = :Bisection)
        @test settings.trace == 2
        @test settings.ic_solver == :Bisection
        @test settings.rlsim === run_sim
    end

    # saCL options
    @testset "saCL options" begin
        settings = OptSettings(
            hmin_sa = 0.1,
            Nfixed_sa = 200,
            Afixed_sa = 0.2,
            Amin_sa = 0.05,
            Amax_sa = 50.0,
            delta_sa = 0.05,
            q_sa = 0.6,
            gamma_sa = 0.01,
            Nmin_sa = 2000,
            z_sa = 2.0,
            Cmrl_sa = 5.0,
            maxiter_sa = 200_000,
            verbose_sa = true
        )
        @test settings.hmin_sa ≈ 0.1
        @test settings.Nfixed_sa == 200
        @test settings.Afixed_sa ≈ 0.2
        @test settings.Amin_sa ≈ 0.05
        @test settings.Amax_sa ≈ 50.0
        @test settings.delta_sa ≈ 0.05
        @test settings.q_sa ≈ 0.6
        @test settings.gamma_sa ≈ 0.01
        @test settings.Nmin_sa == 2000
        @test settings.z_sa ≈ 2.0
        @test settings.Cmrl_sa ≈ 5.0
        @test settings.maxiter_sa == 200_000
        @test settings.verbose_sa === true
    end

    # Bisection options
    @testset "Bisection options" begin
        settings = OptSettings(
            hmin_bi = 0.01,
            hmax_bi = 100.0,
            maxiter_bi = 100,
            nsims_bi = 20_000,
            trunc_bi = 5_000,
            x_tol_bi = 1e-07,
            f_tol_bi = 2.0,
            verbose_bi = false,
            inflate_bi = 1.1
        )
        @test settings.hmin_bi ≈ 0.01
        @test settings.hmax_bi ≈ 100.0
        @test settings.maxiter_bi == 100
        @test settings.nsims_bi == 20_000
        @test settings.trunc_bi == 5_000
        @test settings.x_tol_bi ≈ 1e-07
        @test settings.f_tol_bi ≈ 2.0
        @test settings.verbose_bi === false
        @test settings.inflate_bi ≈ 1.1
    end

    # Global parameter optimization options
    @testset "Global parameter optimization options" begin
        settings = OptSettings(
            x_tol_opt = 1e-06,
            nsims_opt = 5_000,
            maxiter_opt = 500,
            method_opt = :LN_BOBYQA,
            verbose_opt = 2,
            minpar_opt = [0.0, -Inf],
            maxpar_opt = [10.0, Inf])

        @test settings.x_tol_opt ≈ 1e-06
        @test settings.nsims_opt == 5_000
        @test settings.maxiter_opt == 500
        @test settings.method_opt == :LN_BOBYQA
        @test settings.verbose_opt == 2
        @test settings.minpar_opt ≈ [0.0, -Inf]
        @test settings.maxpar_opt ≈ [10.0, Inf]
    end

    # Grid settings
    @testset "Grid settings" begin
        settings = OptSettings(m_grid = 5)
        @test settings.m_grid == 5
    end

    # SPSA settings
    @testset "SPSA settings" begin
        settings = OptSettings(
            initial_step_size_spsa = 0.1,
            expected_n_loss_eval_spsa = 500,
            n_adaptive_gain_spsa = 10,
            gamma_spsa = 0.05
        )
        @test settings.initial_step_size_spsa ≈ 0.1
        @test settings.expected_n_loss_eval_spsa == 500
        @test settings.n_adaptive_gain_spsa == 10
        @test settings.gamma_spsa ≈ 0.05
    end

    # Constraint validations
    @testset "Constraint validations" begin
        @test_throws AssertionError OptSettings(trace = -1)
        @test_throws AssertionError OptSettings(ic_solver = :InvalidSolver)
        @test_throws AssertionError OptSettings(hmin_sa = 0.0)
        @test_throws AssertionError OptSettings(Nfixed_sa = 0)
        @test_throws AssertionError OptSettings(Afixed_sa = 0.0)
        @test_throws AssertionError OptSettings(Amin_sa = 0.0)
        @test_throws AssertionError OptSettings(Amax_sa = 0.0)
        @test_throws AssertionError OptSettings(delta_sa = 0.0)

        @test_throws AssertionError OptSettings(hmin_bi = 0.0)
        @test_throws AssertionError OptSettings(hmin_bi = 1.0, hmax_bi=0.5)
        @test_throws AssertionError OptSettings(maxiter_bi = 0)
        @test_throws MethodError OptSettings(nsims_bi = 0.0)
        @test_throws AssertionError OptSettings(nsims_bi = 0)
        @test_throws MethodError OptSettings(trunc_bi = 0.0)
        @test_throws AssertionError OptSettings(trunc_bi = 0)
        @test_throws AssertionError OptSettings(x_tol_bi= 0.0)
        @test_throws AssertionError OptSettings(f_tol_bi = 0.0)

        @test_throws AssertionError OptSettings(initial_step_size_spsa = 0.0)
        @test_throws AssertionError OptSettings(expected_n_loss_eval_spsa = 0)
        @test_throws AssertionError OptSettings(n_adaptive_gain_spsa = 0)
        @test_throws AssertionError OptSettings(gamma_spsa = 0.0)
    end
end