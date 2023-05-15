using SPM
using Test

# Test OptSettings struct
@testset "OptSettings struct" begin
    # Happy path
    @testset "Happy path" begin
        opt = OptSettings(trace=1, ic_solver=:SA)
        @test opt.trace == 1
        @test opt.ic_solver == :SA
        @test opt.rlsim == run_sim_sa
        @test opt.hmin_sa == sqrt(eps())
        @test opt.Nfixed_sa == 500
        @test opt.Afixed_sa == 0.1
        @test opt.Amin_sa == 0.1
        @test opt.Amax_sa == 100.0
        @test opt.delta_sa == 0.1
        @test opt.q_sa == 0.55
        @test opt.gamma_sa == 0.02
        @test opt.Nmin_sa == 1000
        @test opt.z_sa == 3.0
        @test opt.Cmrl_sa == 10.0
        @test opt.maxiter_sa == 100_000
        @test opt.hmin_bi == sqrt(eps())
        @test opt.hmax_bi == 50.0
        @test opt.maxiter_bi == 50
        @test opt.nsims_bi == 10_000
        @test opt.trunc_bi == 10_000
        @test opt.x_tol_bi == 1e-06
        @test opt.f_tol_bi == 1.0
        @test opt.inflate_bi == 1.05
        @test opt.x_tol_opt == 1e-05
        @test opt.nsims_opt == 10_000
        @test opt.maxiter_opt == 100
        @test opt.method_opt == :LN_BOBYQA
        @test opt.verbose_opt == 1
        @test opt.minpar_opt == [-Inf]
        @test opt.maxpar_opt == [Inf]
        @test opt.m_grid == 10
    end

    # Edge cases
    @testset "Edge cases" begin
        # trace is negative
        @test_throws AssertionError OptSettings(trace=-1, ic_solver=:SA)

        # ic_solver is not a symbol
        @test_throws MethodError OptSettings(trace=1, ic_solver="SA")

        # ic_solver is not a valid option
        @test_throws AssertionError OptSettings(trace=1, ic_solver=:invalid)

        # hmin_sa is negative
        @test_throws AssertionError OptSettings(trace=1, ic_solver=:SA, hmin_sa=-1.0)

        # Nfixed_sa is negative
        @test_throws AssertionError OptSettings(trace=1, ic_solver=:SA, Nfixed_sa=-1)

        # Afixed_sa is negative
        @test_throws AssertionError OptSettings(trace=1, ic_solver=:SA, Afixed_sa=-1.0)

        # Amin_sa is negative
        @test_throws AssertionError OptSettings(trace=1, ic_solver=:SA, Amin_sa=-1.0)

        # Amax_sa is negative
        @test_throws AssertionError OptSettings(trace=1, ic_solver=:SA, Amax_sa=-1.0)

        # delta_sa is negative
        @test_throws AssertionError OptSettings(trace=1, ic_solver=:SA, delta_sa=-1.0)
    end
end