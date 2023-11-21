abstract type NominalProperties end

get_value(NM::NominalProperties) = NM.value

"""
    ARL(value)

Value of the in-control average run length of the control chart, i.e. if ``RL = \\inf\\{t > 0: \\text{Chart detects OC}\\}`` is the run length, then the average run length ``ARL`` is

``ARL = \\mathbb{E}[RL|\\tau = +\\infty]``,

where ``\\{\\tau = +\\infty\\}`` represents the process being in-control.
"""
@with_kw struct ARL <: NominalProperties
    value::Float64   
    @assert value > 0
end
export ARL


"""
    QRL(value, qtl)

Value of the in-control quantile of the run length of the control chart, i.e. if ``RL = \\inf\\{t > 0 : \\text{Chart detects OC}\\}`` is the run length, then `value` is the value of the `qtl`-level quantile of the distribution of ``RL`` if the process is in-control.
"""
@with_kw struct QRL <: NominalProperties
    value::Float64
    qtl::Float64 = 0.5
    @assert value > 0
    @assert 0.0 < qtl < 1.0
end
export QRL
