abstract type AbstractPhase1 end

struct Phase1Data{T} <: AbstractPhase1 where T <: AbstractVecOrMat
    x::T
end
export Phase1Data


"""
    new_data(P1::Phase1Data)

Generates a new observation via bootstrap, based on the observed Phase I data.
"""
function new_data(P1::Phase1Data{T}) where T <: AbstractVector
    return P1.x[rand(1:length(P1.x))]
end

function new_data(P1::Phase1Data{T}) where T <: AbstractMatrix
    return view(P1.x, rand(1:size(P1.x)[1]), :)
end
export new_data