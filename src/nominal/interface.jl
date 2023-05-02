abstract type NominalProperties end

get_value(A::NominalProperties) = A.value

@with_kw struct ARL <: NominalProperties
    value::Float64   
    @assert value > 0
end
export ARL


@with_kw struct QRL <: NominalProperties
    value::Float64
    qtl::Float64 = 0.5
    @assert value > 0
    @assert 0.0 < qtl < 1.0
end
export QRL
