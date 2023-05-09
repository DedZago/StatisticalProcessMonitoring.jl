abstract type AbstractPhase1{T} end

"""
    new_data(P1::AbstractPhase1{T})
    new_data(P1::AbstractPhase1{AbstractVector})
    new_data(P1::AbstractPhase1{AbstractMatrix})

Generates a new observation based on the observed Phase I (in-control) data.
If it is not overloaded, then it defaults to generating data using a nonparametric bootstrap.
"""
function new_data(P1::AbstractPhase1{T}) where T <: AbstractVector
    return P1.x[rand(1:length(P1.x))]
end

function new_data(P1::AbstractPhase1{T}) where T <: AbstractMatrix
    return view(P1.x, rand(1:size(P1.x)[1]), :)
end
export new_data


struct Phase1Data{T} <: AbstractPhase1{T}
    x::T
end
export Phase1Data



################# TEST TRUE DGP #################
struct Phase2Distribution{T} <: AbstractPhase1{T}
    dist::T
end
export Phase2Distribution

new_data(DGP::Phase2Distribution) = rand(DGP.dist)