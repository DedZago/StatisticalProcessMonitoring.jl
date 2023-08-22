@with_kw struct FunctionalObservation{A, B}
    x::Vector{A}
    y::Vector{B}
    @assert length(x) == length(y)
end
export FunctionalObservation

FunctionalObservation(x::Float64, y::Float64) = FunctionalObservation([x], [y])


const FunctionalData{A,B} = Vector{FunctionalObservation{A,B}} where {A,B}
export FunctionalData

function FunctionalData(X::AbstractVector, Y::AbstractVector)
    @assert size(X) == size(Y) "Size of X <$(size(X))> is different from size of Y <$(size(Y))>"
    return [FunctionalObservation(X[i], Y[i]) for i in 1:size(X,1)]
end

function FunctionalData(X::AbstractMatrix, Y::AbstractMatrix)
    @assert size(X) == size(Y) "Size of X <$(size(X))> is different from size of Y <$(size(Y))>"
    return [FunctionalObservation(X[i,:], Y[i,:]) for i in 1:size(X,1)]
end
export FunctionalData