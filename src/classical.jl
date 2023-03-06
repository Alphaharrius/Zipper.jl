import Base: +, -, ==

module classical

using LinearAlgebra

"""
Represents a classical space with basis that spans the space.
"""
struct ClassicalSpace
    """ The basis vectors of this space stacked in columns. """
    basis::Matrix{Float64}
    """ The dimension of this space. """
    dimension::Int8
end

"""
A general point in a [`ClassicalSpace`](@ref).
"""
struct Point
    """ The position in spacial basis in the `parent_space`. """
    position::Vector{Float64}
    """ The parent space that includes this point. """
    parent_space::ClassicalSpace
end

"""
    create_classical_space(basis)

Create a new [`ClassicalSpace`](@ref) with `basis`.
"""
function create_classical_space(basis::AbstractArray{Float64})::ClassicalSpace
    return ClassicalSpace(Matrix(basis), size(basis, 1))
end

function euclidean(dimension::Int8)::ClassicalSpace
    return create_classical_space(Matrix{Float64}(I, dimension, dimension))
end

function phys_operator(space::ClassicalSpace)::Function
    basis::Matrix{Float64} = space.basis
    euclidean_space::ClassicalSpace = euclidean(space.dimension)
    return function (point::Point)
        return Point(basis * point.position, euclidean_space)
    end
end

"""
    create_point(position[, parent_space])

Create a new [`Point`](@ref) of `position` in the [`ClassicalSpace`](@ref) of `parent_space`.
"""
function create_point(position::Vector{Float64}, parent_space::ClassicalSpace)::Point
    @assert(length(position) == parent_space.dimension)
    return Point(position, parent_space)
end

function Base.:+(a::Point, b::Point)::Point
    @assert(typeof(a.parent_space) == typeof(b.parent_space))
    return Point(a.position + b.position, a.parent_space)
end

function Base.:-(a::Point, b::Point)::Point
    @assert(typeof(a.parent_space) == typeof(b.parent_space))
    return Point(a.position - b.position, a.parent_space)
end

function Base.:(==)(a::Point, b::Point)::Bool
    return a.parent_space == b.parent_space && isapprox(a.position, b.position)
end

function __init__()
end

end
