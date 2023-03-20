module spaces

using LinearAlgebra

abstract type AbstractSpace end

"""
    basis(space::AbstractSpace)

Retrieve the basis of an `AbstractSpace`, this method assumed the existance of the attribute
`representation::Matrix{Float64}` within all subtype of `AbstractSpace`.
"""
basis(space::AbstractSpace)::Matrix{Float64} = space.representation

"""
    dimension(space::AbstractSpace)

Retieve the dimension spanned by the basis vectors of `space`.
"""
dimension(space::AbstractSpace)::Int64 = size(basis(space), 1)

function Base.:(==)(a::AbstractSpace, b::AbstractSpace)::Bool
    return typeof(a) == typeof(b) && isapprox(basis(a), basis(b))
end

"""
    RealSpace(basis::AbstractArray{Float64})

Create a new `RealSpace` with `basis`.
"""
struct RealSpace <: AbstractSpace
    """ The basis vectors of this space stacked in columns. """
    representation::Matrix{Float64}

    RealSpace(basis::AbstractArray{Float64}) = new(Matrix(basis))
end

"""
    MomentumSpace(real_space::RealSpace)

Create a new `MomentumSpace` from a `RealSpace`.
"""
struct MomentumSpace <: AbstractSpace
    """ The basis vectors of this space stacked in columns. """
    representation::Matrix{Float64}

    # MomentumSpace(real_space::RealSpace) = new(Matrix(transpose(inv(basis(real_space)))))
end

Base.:convert(::Type{T}, source::T) where {T <: AbstractSpace} = source
Base.:convert(::Type{T}, source::F) where {T <: MomentumSpace, F <: RealSpace} = MomentumSpace(2.0 * π * Matrix(transpose(inv(basis(source)))))
Base.:convert(::Type{T}, source::F) where {T <: RealSpace, F <: MomentumSpace} = RealSpace(Matrix(transpose(inv(basis(source) / (2.0 * π)))))

abstract type AbstractSpaceSubset end

"""
    Point(position, parent_space)

Create a `Point` of `position` in the `Space` of `parent_space`.
"""
struct Point <: AbstractSpaceSubset
    """ The position in spacial basis in the `parent_space`. """
    position::Vector{Rational{Int64}} # Julia can compare rational numbers.
    """ The parent space that includes this point. """
    parent_space::AbstractSpace

    Point(position::Vector{Rational{Int64}}, parent_space::AbstractSpace) = (
        length(position) == dimension(parent_space) ? new(position, parent_space) : error("Dimension mismatch."))
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

struct SpaceSubset <: AbstractSpaceSubset
    representation::Set{AbstractSpaceSubset}
    parent_space::AbstractSpace

    function SpaceSubset(elements::Set{AbstractSpaceSubset})::SpaceSubset
        parent_space::AbstractSpace = first(elements).parent_space
        @assert(all(element.parent_space == parent_space for element in elements))
        return new(elements, parent_space)
    end
end

""" Reflexive. """
Base.:convert(::Type{T}, source::T) where {T <: SpaceSubset} = source
Base.:convert(::Type{T}, source::F) where {T <: SpaceSubset, F <: Point} = SpaceSubset(Set{AbstractSpaceSubset}([source]))

""" Chowning `SpaceSubset` as the common base for all subtypes of `AbstractSpaceSubset`. """
Base.promote_rule(::Type{T}, ::Type{F}) where {T <: SpaceSubset, F <: AbstractSpaceSubset} = SpaceSubset
Base.promote_rule(::Type{T}, ::Type{F}) where {T <: AbstractSpaceSubset, F <: SpaceSubset} = SpaceSubset

Base.:(==)(a::SpaceSubset, b::SpaceSubset)::Bool = a.parent_space == b.parent_space && a.representation == b.representation
Base.:union(a::SpaceSubset, b::SpaceSubset)::SpaceSubset = SpaceSubset(union(a.representation, b.representation))
Base.:intersect(a::SpaceSubset, b::SpaceSubset)::SpaceSubset = SpaceSubset(intersect(a.representation, b.representation))

Base.:(==)(a::AbstractSpaceSubset, b::AbstractSpaceSubset)::Bool = (==)(promote(a, b) ...)
Base.:union(a::AbstractSpaceSubset, b::AbstractSpaceSubset)::SpaceSubset = union(promote(a, b) ...)
Base.:intersect(a::AbstractSpaceSubset, b::AbstractSpaceSubset)::SpaceSubset = intersect(promote(a, b) ...)

function __init__()
end

end
