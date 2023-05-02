abstract type NominalProperties end
#TODO: testing

get_value(A::NominalProperties) = A.value

struct ARL <: NominalProperties
    value::Float64   
end
export ARL

@with_kw struct QRL <: NominalProperties
    value::Float64
    qtl::Float64 = 0.5
end
export QRL
