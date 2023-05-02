abstract type AbstractPhase1 end
#TODO: testing

struct Phase1Data{T} <: AbstractPhase1 where T <: AbstractVecOrMat
    x::T
end
export Phase1Data

function new_data(P1::Phase1Data{T}) where T <: AbstractVector
    return P1.x[rand(1:length(P1.x))]
end

function new_data(P1::Phase1Data{T}) where T <: AbstractMatrix
    return P1.x[rand(1:size(P1.x)[1]), :]
end
export new_data