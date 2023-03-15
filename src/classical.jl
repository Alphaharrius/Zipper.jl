import Base: +, -, ==

module classical

using LinearAlgebra

"""
    Classical(basis)

Create a new `Space` with `basis`.
"""
struct Space
    """ The basis vectors of this space stacked in columns. """
    basis::Matrix{Float64}
    """ The dimension of this space. """
    dimension::Int8

    Space(basis::AbstractArray{Float64}) = new(Matrix(basis), size(basis, 1))
end

"""
    Point(position, parent_space)

Create a `Point` of `position` in the `Space` of `parent_space`.
"""
struct Point
    """ The position in spacial basis in the `parent_space`. """
    position::Vector{Rational{Int64}} # Julia can compare rational numbers.
    """ The parent space that includes this point. """
    parent_space::Space

    Point(position::Vector{Rational{Int64}}, parent_space::Space) = (
        length(position) == parent_space.dimension ? new(position, parent_space) : error("Dimension mismatch."))
end

function Base.:(==)(a::Space, b::Space)
    return isapprox(a.basis, b.basis)
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
    return a.parent_space == b.parent_space && a.position == b.position
end

function __init__()
end

end
