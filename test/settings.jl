using SPM
using Test

# Test suite for OptSettings
@testset "OptSettings Tests" begin
    # Test default construction
    @testset "Default Construction" begin
        opt_settings = OptSettings{Float64, Int, Bool}()
        @test opt_settings.x_tol == 1e-03
        @test opt_settings.nsims == 1000
        @test opt_settings.maxiter == 1000
        @test opt_settings.verbose == true
        @test opt_settings.minpar == [-Inf]
        @test opt_settings.maxpar == [Inf]
        @test opt_settings.m_grid == 10
        @test opt_settings.initial_step_size_spsa == 0.05
        @test opt_settings.expected_n_loss_eval_spsa == 100
        @test opt_settings.n_adaptive_gain_spsa == 20
        @test opt_settings.gamma_spsa == 0.02
    end

    # Test construction with custom arguments
    @testset "Custom Construction" begin
        opt_settings = OptSettings{Float64, Int, Bool}(
            x_tol = 1e-04,
            nsims = 500,
            maxiter = 500,
            verbose = false,
            minpar = [0.0, -1.0],
            maxpar = [1.0, 2.0],
            m_grid = 5,
            initial_step_size_spsa = 0.1,
            expected_n_loss_eval_spsa = 200,
            n_adaptive_gain_spsa = 30,
            gamma_spsa = 0.1
        )
        @test opt_settings.x_tol == 1e-04
        @test opt_settings.nsims == 500
        @test opt_settings.maxiter == 500
        @test opt_settings.verbose == false
        @test opt_settings.minpar == [0.0, -1.0]
        @test opt_settings.maxpar == [1.0, 2.0]
        @test opt_settings.m_grid == 5
        @test opt_settings.initial_step_size_spsa == 0.1
        @test opt_settings.expected_n_loss_eval_spsa == 200
        @test opt_settings.n_adaptive_gain_spsa == 30
        @test opt_settings.gamma_spsa == 0.1
    end

    # Test construction using OptSettings function
    @testset "Construction using OptSettings function" begin
        p = 3
        stat = DiagMEWMA(Λ = fill(0.2, p))
        lim = OneSidedFixedLimit(0.2, true)
        nom = ARL(200)
        x = randn(200, 3)
        ph2 = Phase2(Bootstrap(), x)
        ch = ControlChart(stat, lim, nom, ph2)
        opt_settings = OptSettings(ch, x_tol = 1e-04, nsims = 500)
        @test opt_settings.x_tol == 1e-04
        @test opt_settings.nsims == 500
        @test opt_settings.maxiter == 1000  # Default value
        @test opt_settings.verbose == true   # Default value
        @test opt_settings.minpar == fill(-Inf, p)
        @test opt_settings.maxpar == fill(Inf, p)
        @test opt_settings.m_grid == 10      # Default value
        @test opt_settings.initial_step_size_spsa == 0.05  # Default value
        @test opt_settings.expected_n_loss_eval_spsa == 100  # Default value
        @test opt_settings.n_adaptive_gain_spsa == 20       # Default value
        @test opt_settings.gamma_spsa == 0.02               # Default value
    end
end