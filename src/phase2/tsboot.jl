using Distributions 

"""
    BlockBootstrap{T} <: AbstractSampling

Represents a (circular) block bootstrap sampling method.

# Fields
- `block::T`: The current block of data being sampled from.
- `blocksize::Int`: The size of each block.
- `t::Int`: The current index within the block.

# Constructors
- `BlockBootstrap(blocksize::Int, data::Vector{T}) where T`: Constructs a `BlockBootstrap` object for a vector of data.
- `BlockBootstrap(blocksize::Int, data::Matrix{T}) where T`: Constructs a `BlockBootstrap` object for a matrix of data.
"""
mutable struct BlockBootstrap{T} <: AbstractSampling
    block::T
    const blocksize::Int
    t::Int

    # Initialized with t=0 so that update_block! updates the block the first time it is called
    function BlockBootstrap(blocksize::Int, data::Vector{T}) where T
        @assert blocksize > 0
        out = new{Vector{T}}(similar(data)[1:blocksize], blocksize, 0)
        return out
    end

    function BlockBootstrap(blocksize::Int, data::Matrix{T}) where T 
        @assert blocksize > 0
        out = new{Matrix{T}}(zeros(T, blocksize, size(data)[2]), blocksize, 0)
        return out
    end
end
export BlockBootstrap

get_block(B::BlockBootstrap) = B.block
export get_block
get_blocksize(B::BlockBootstrap) = B.blocksize
export get_blocksize
get_counter(B::BlockBootstrap) = B.t
export get_counter
set_counter!(B::BlockBootstrap, t) = B.t = t
export set_counter!

new_data(B::BlockBootstrap, data::AbstractVector) = get_block(B)[get_counter(B)]
new_data(B::BlockBootstrap, data::AbstractMatrix) = view(get_block(B), get_counter(B), :)

function update_block!(B::BlockBootstrap, data::AbstractVector)
    set_counter!(B, (get_counter(B) % get_blocksize(B)) + 1)
    if get_counter(B) <= 1
        bb = sample(1:length(data))
        for i in 1:get_blocksize(B)
            # wrap around the circle
            get_block(B)[i] = data[(bb + i - 1) % length(data) + 1]
        end
    end
end

function update_block!(B::BlockBootstrap, data::AbstractMatrix)
    set_counter!(B, get_counter(B) % get_blocksize(B) + 1)
    if get_counter(B) <= 1
        bb = sample(1:length(data))
        for i in 1:get_blocksize(B)
            # wrap around the circle
            get_block(B)[i, :] = data[(bb + i - 1) % size(data)[1] + 1, :]
        end
    end
end

function new_data!(B::BlockBootstrap, data::AbstractVecOrMat)
    update_block!(B, data)
    return new_data(B, data) 
end


"""
    StationaryBootstrap{T} <: AbstractSampling

Represents a stationary block bootstrap sampling method, where the block length is sampled from a `Geometric` random variable.

# Fields
- `block::T`: The current block of data being sampled from.
- `blocksize::Int`: The average size of each block.
- `t::Int`: The current index within the block.

# Constructors
- `StationaryBootstrap(blocksize::Int, data::Vector{T}) where T`: Constructs a `StationaryBootstrap` object for a vector of data.
- `StationaryBootstrap(blocksize::Int, data::Matrix{T}) where T`: Constructs a `StationaryBootstrap` object for a matrix of data.
"""
mutable struct StationaryBootstrap{T} <: AbstractSampling
    block::T
    blocksize::Int
    t::Int

    # Initialized with t=0 so that update_block! updates the block the first time it is called
    function StationaryBootstrap(blocksize::Int, data::Vector{T}) where T
        @assert blocksize > 0
        out = new{Vector{T}}(zeros(T, blocksize), blocksize, 0)
        return out
    end

    function StationaryBootstrap(blocksize::Int, data::Matrix{T}) where T 
        @assert blocksize > 0
        out = new{Matrix{T}}(zeros(T, blocksize, size(data)[2]), blocksize, 0)
        return out
    end
end
export StationaryBootstrap

set_block!(B::StationaryBootstrap{T}, x::T) where T = B.block = x 
get_block(B::StationaryBootstrap) = B.block
get_blocksize(B::StationaryBootstrap) = B.blocksize
get_counter(B::StationaryBootstrap) = B.t
set_counter!(B::StationaryBootstrap, t) = B.t = t

new_data(B::StationaryBootstrap, data::AbstractVector) = get_block(B)[get_counter(B)]
new_data(B::StationaryBootstrap, data::AbstractMatrix) = view(get_block(B), get_counter(B), :)

function new_data!(B::StationaryBootstrap, data::AbstractVecOrMat)
    update_block!(B, data)
    return new_data(B, data) 
end

function update_block!(B::StationaryBootstrap, data::AbstractVector)
    set_counter!(B, (get_counter(B) % length(get_block(B))) + 1)
    if get_counter(B) <= 1
        ll = rand(1 + Geometric(1.0/get_blocksize(B)))
        set_block!(B, zeros(ll))
        bb = sample(1:length(data))
        for i in 1:ll
            # wrap around the circle
            get_block(B)[i] = data[(bb + i - 1) % length(data) + 1]
        end
    end
end

function update_block!(B::StationaryBootstrap, data::AbstractMatrix)
    set_counter!(B, (get_counter(B) % size(get_block(B))[1]) + 1)
    if get_counter(B) <= 1
        ll = rand(1 + Geometric(1.0/get_blocksize(B)))
        set_block!(B, zeros(ll, size(data)[2]))
        bb = sample(1:length(data))
        for i in 1:ll
            # wrap around the circle
            get_block(B)[i, :] = data[(bb + i - 1) % size(data)[2] + 1, :]
        end
    end
end

