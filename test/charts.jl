module TestCharts
using SPM
using Test

# ASSUMTION: `e1` and `e2` have the same run-time type
# CHECK FIELDS ARE EQUAL BETWEEN TWO STRUCTS
@generated structEqual(e1, e2) = begin
    if fieldcount(e1) == 0
        return :(true)
    end
    mkEq    = fldName -> :(e1.$fldName == e2.$fldName)
    # generate individual equality checks
    eqExprs = map(mkEq, fieldnames(e1))
    # construct &&-expression for chaining all checks
    mkAnd   = (expr, acc) -> Expr(:&&, expr, acc)
    # no need in initial accumulator because eqExprs is not empty
    foldr(mkAnd, eqExprs)
end

@testset "Charts" begin
    x = randn(100)
    NM = ARL(200)
    PH1 = Phase1Data(Bootstrap(), x)
    LIM = TwoSidedFixedLimit(1.0)
    STAT = EWMA(λ = 0.2)
    @testset "constructor" begin
        CH = ControlChart(STAT, LIM, NM, PH1)
        CH = ControlChart(STAT, LIM, NM, PH1, 0)
        @test is_IC(CH)
        @test !is_OC(CH)
        @test get_t(CH) == 0
        @test get_design(CH) == (λ = 0.2,)
        @test get_phase1(CH) == PH1
        @test get_limit_value(CH) == 1.0 * [-1, 1]
        @test get_value(CH) == 0.0
        @test get_statistic(CH) == STAT
        @test get_nominal(CH) == NM
        @test get_nominal_value(CH) == 200.0
        update_chart!(CH, 1.0)
        @test get_value(CH) == 0.2
        update_chart(CH, 1.0)
        @test get_value(CH) == 0.2
        update_chart!(CH, 100.0)
        @test is_OC(CH)
        @test !is_IC(CH)
        set_nominal!(CH, NM)
        set_phase1!(CH, PH1)
        set_statistic!(CH, STAT)
        set_limit!(CH, LIM)
        set_limit!(CH, 0.1)
        @test get_limit_value(CH) == 0.1 * [-1, 1]
        set_design!(CH, 0.3)
        @test get_design(CH)[1] == 0.3
        @test get_maxrl(CH) == Inf
        @test new_data(CH) in x
    end

    @testset "run_sim" begin
        x = randn(100)
        NM = ARL(200)
        PH1 = Phase1Data(Bootstrap(), x)
        LIM = TwoSidedFixedLimit(1.0)
        STAT = EWMA(λ = 0.2)
        CH = ControlChart(STAT, LIM, NM, PH1)
        CH_ = deepcopy(CH)
        run_sim(CH)
        @test run_sim(CH, maxiter = 1) == 1
        @test structEqual(get_statistic(CH), get_statistic(CH_))
        @test structEqual(get_limit(CH), get_limit(CH_))
        @test structEqual(get_nominal(CH), get_nominal(CH_))
        @test structEqual(get_sampler(get_phase1(CH)), get_sampler(get_phase1(CH_)))
        @test structEqual(get_data(get_phase1(CH)), get_data(get_phase1(CH_)))
        rsa = run_sim_sa(CH, maxiter=Inf, delta=0.0)
        @test length(rsa) == 3
        @test allequal(collect(rsa))

        # --- Block bootstrap
        x = randn(100)
        NM = ARL(200)
        PH1 = Phase1Data(BlockBootstrap(5, x), x)
        LIM = TwoSidedFixedLimit(1.0)
        STAT = EWMA(λ = 0.2)
        CH = ControlChart(STAT, LIM, NM, PH1)
        CH_ = deepcopy(CH)
        run_sim(CH)
        @test run_sim(CH, maxiter = 1) == 1
        @test structEqual(get_statistic(CH), get_statistic(CH_))
        @test structEqual(get_limit(CH), get_limit(CH_))
        @test structEqual(get_nominal(CH), get_nominal(CH_))
        @test structEqual(get_sampler(get_phase1(CH)), get_sampler(get_phase1(CH_)))
        @test structEqual(get_data(get_phase1(CH)), get_data(get_phase1(CH_)))
        rsa = run_sim_sa(CH, maxiter=Inf, delta=0.0)
        @test length(rsa) == 3
        @test allequal(collect(rsa))
    end

    x = randn(100)
    NM = ARL(200)
    PH1 = Phase1Data(Bootstrap(), x)
    λ = 0.2
    STAT = EWMA(λ = λ)
    f(t, STAT) = sqrt(STAT.λ/(2.0 - STAT.λ) * (1.0 - (1.0 - STAT.λ)^(2.0*t)))
    LIM = OneSidedCurvedLimit(1.0, true, f)

    @testset "Curved chart" begin
        CH = ControlChart(STAT, LIM, NM, PH1)
        CH = ControlChart(STAT, LIM, NM, PH1, 0)
        @test get_t(CH) == 0
        @test get_design(CH) == (λ = 0.2,)
        @test get_phase1(CH) == PH1
        @test get_limit_value(CH) == 0.0
        @test get_limit_value(CH) != get_value(get_limit(CH))
        @test get_value(CH) == 0.0
        @test get_statistic(CH) == STAT
        @test get_nominal(CH) == NM
        @test get_nominal_value(CH) == 200.0
        xnew = 0.3
        update_chart!(CH, xnew)
        @test is_IC(CH)
        @test !is_OC(CH)
        @test get_value(CH) == λ * xnew
        update_chart!(CH, 100.0)
        @test is_OC(CH)
        @test !is_IC(CH)
        set_nominal!(CH, NM)
        set_phase1!(CH, PH1)
        set_statistic!(CH, STAT)
        set_limit!(CH, LIM)
        set_limit!(CH, 0.1)
        @test get_limit_value(CH) != 0.0
        @test get_limit_value(CH) != get_value(get_limit(CH))
        set_design!(CH, 0.3)
        @test get_design(CH)[1] == 0.3
        @test get_maxrl(CH) == Inf
        @test new_data(CH) in x
        
    end

    @testset "Multiple charts" begin
        NM = QRL(200, 0.2)
        PH1 = Phase1Data(Bootstrap(), x)
        h = 1.0
        LIM = OneSidedFixedLimit(h, true)
        λ1 = 0.2
        λ2 = 0.5
        STAT1 = EWMA(λ = λ1)
        STAT2 = EWMA(λ = λ2)
        CH = ControlChart([STAT1, STAT2], [deepcopy(LIM), deepcopy(LIM)], NM, PH1)
        @test get_t(CH) == 0
        @test typeof(get_design(CH)) <: Vector
        @test length(get_design(CH)) == 2
        @test get_design(CH)[1] == (λ = λ1,)
        @test get_design(CH)[2] == (λ = λ2,)
        @test get_phase1(CH) == PH1
        @test get_limit_value(CH) == fill(h, 2)
        @test isa(get_limit_value(CH), Vector)
        @test get_value(CH) == zeros(2)
        @test typeof(get_statistic(CH)) <: Vector
        @test length(get_statistic(CH)) == 2
        @test get_nominal(CH) == NM
        @test get_nominal_value(CH) == 200.0
        xnew = 0.3
        update_chart!(CH, xnew)
        val = get_value(CH)
        @test val[1] < val[2]
        @test is_IC(CH)
        @test !is_OC(CH)
        @test get_value(CH) == [λ1*xnew, λ2*xnew]
        update_chart!(CH, 100.0)
        @test is_OC(CH)
        @test !is_IC(CH)
        set_nominal!(CH, NM)
        set_phase1!(CH, PH1)
        newh = 0.1
        set_limit!(CH, newh)
        @test get_limit_value(CH) == fill(newh, 2)
        newh = [0.1, 0.5]
        set_limit!(CH, newh)
        @test get_limit_value(CH) == newh
        set_design!(CH, [0.5, 0.7])
        @test get_design(CH)[1][1] == 0.5
        @test get_design(CH)[2][1] == 0.7
        @test get_maxrl(CH) == Inf
        @test new_data(CH) in x
    end
    @testset "shallow copy" begin
        @testset "ControlChart" begin
            NM = ARL(200)
            PH1 = Phase1Data(Bootstrap(), x)
            h = 1.0
            LIM = TwoSidedFixedLimit(h)
            k = 0.5
            val = 0.0
            STAT = CUSUM(k = k, value=val)
            CH = ControlChart(STAT, LIM, NM, PH1)
            CH_ = shallow_copy_sim(CH)
            knew = 0.2
            set_design!(CH_, knew)
            xnew = 1.0
            update_chart!(CH_, xnew)
            set_limit!(CH_, h*0.5)
            set_nominal!(CH_, ARL(300))
            @test get_value(CH) == val
            @test get_design(CH)[1] == k
            @test get_value(CH) != get_value(CH_)
            @test get_design(CH) != get_design(CH_)
            @test get_nominal(CH) != get_nominal(CH_)
            @test get_nominal(CH) != get_nominal(CH_)
            update_chart!(CH_, 10^5)
            @test is_IC(CH)
            @test !is_OC(CH)
            @test is_OC(CH_)
            @test !is_IC(CH_)
        end

        @testset "MultipleChart" begin
            NM = ARL(200)
            PH1 = Phase1Data(Bootstrap(), x)
            h = 1.0
            LIM = TwoSidedFixedLimit(h)
            k = 0.5
            val = 0.0
            STAT1 = CUSUM(k = k, value=val)
            λ = 0.2
            STAT2 = EWMA(λ = λ, value=val)
            CH = ControlChart([STAT1,STAT2], [deepcopy(LIM), deepcopy(LIM)], NM, PH1)
            CH_ = shallow_copy_sim(CH)
            knew = 0.2
            λnew = 0.5
            set_design!(CH_, [knew, λnew])
            xnew = 1.0
            update_chart!(CH_, xnew)
            set_limit!(CH_, h*0.5)
            set_nominal!(CH_, ARL(300))
            @test get_value(CH) == [val, val]
            @test get_design(CH)[1][1] == k
            @test get_design(CH)[2][1] == λ
            @test get_value(CH) != get_value(CH_)
            @test get_design(CH) != get_design(CH_)
            @test get_nominal(CH) != get_nominal(CH_)
            @test get_nominal(CH) != get_nominal(CH_)
            update_chart!(CH_, 10^5)
            @test is_IC(CH)
            @test !is_OC(CH)
            @test is_OC(CH_)
            @test !is_IC(CH_)
        end
    end
end
end