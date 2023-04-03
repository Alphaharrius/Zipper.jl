module Spaces
export representation, basis, dimension, euclidean
export AbstractSpace, ClassicalSpace, RealSpace, MomentumSpace, AbstractSpaceSubset, SpaceSubset

using LinearAlgebra

abstract type AbstractSpace{T} end

representation(space::AbstractSpace{T}) where {T <: Any} = space.representation

abstract type ClassicalSpace <: AbstractSpace{Matrix{Float64}} end

"""
    basis(space::ClassicalSpace)

Retrieve the basis of an `ClassicalSpace`, this method assumed the existance of the attribute
`representation::Matrix{Float64}` within all subtype of `ClassicalSpace`.
"""
basis(space::ClassicalSpace)::Matrix{Float64} = representation(space)

"""
    dimension(space::AbstractSpace)

Retieve the dimension spanned by the basis vectors of `space`.
"""
dimension(space::ClassicalSpace)::Int64 = size(basis(space), 1)

function Base.:(==)(a::ClassicalSpace, b::ClassicalSpace)::Bool
    return typeof(a) == typeof(b) && isapprox(basis(a), basis(b))
end

"""
    RealSpace(basis::AbstractArray{Float64})

Create a new `RealSpace` with `basis`.
"""
struct RealSpace <: ClassicalSpace
    """ The basis vectors of this space stacked in columns. """
    representation::Matrix{Float64}

    RealSpace(basis::AbstractArray{Float64}) = new(Matrix(basis))
end

euclidean(dimension::Int64) = RealSpace(Matrix{Float64}(I(dimension)))

"""
    MomentumSpace(real_space::RealSpace)

Create a new `MomentumSpace` from a `RealSpace`.
"""
struct MomentumSpace <: ClassicalSpace
    """ The basis vectors of this space stacked in columns. """
    representation::Matrix{Float64}
end

Base.:convert(::Type{T}, source::T) where {T <: AbstractSpace} = source

Base.:convert(::Type{MomentumSpace}, source::RealSpace) = MomentumSpace(2.0 * π * Matrix(transpose(inv(basis(source)))))
Base.:convert(::Type{RealSpace}, source::MomentumSpace) = RealSpace(Matrix(transpose(inv(basis(source) / (2.0 * π)))))

abstract type AbstractSpaceSubset end

struct SpaceSubset <: AbstractSpaceSubset
    representation::Set{AbstractSpaceSubset}
    parent_space::AbstractSpace

    function SpaceSubset(elements::Set{T}, space::AbstractSpace)::SpaceSubset where {T <: AbstractSpaceSubset}
        reference_space::AbstractSpace = first(elements).parent_space
        @assert(reference_space == space)
        @assert(all(element.parent_space == reference_space for element in elements))
        return new(elements, reference_space)
    end
end

Base.:length(source::SpaceSubset)::Int64 = length(source.representation)

""" Reflexive. """
Base.:convert(::Type{T}, source::T) where {T <: SpaceSubset} = source

""" Chowning `SpaceSubset` as the common base for all subtypes of `AbstractSpaceSubset`. """
Base.promote_rule(::Type{T}, ::Type{F}) where {T <: AbstractSpaceSubset, F <: AbstractSpaceSubset} = SpaceSubset

Base.:(==)(a::SpaceSubset, b::SpaceSubset)::Bool = a.parent_space == b.parent_space && a.representation == b.representation

function __union(a::SpaceSubset, b::SpaceSubset)::SpaceSubset
    @assert(a.parent_space == b.parent_space)
    return SpaceSubset(union(a.representation, b.representation), a.parent_space)
end

function __intersect(a::SpaceSubset, b::SpaceSubset)::SpaceSubset
    @assert(a.parent_space == b.parent_space)
    return SpaceSubset(intersect(a.representation, b.representation), a.parent_space)
end

Base.:union(a::SpaceSubset, b::SpaceSubset)::SpaceSubset = __union(a, b)
Base.:intersect(a::SpaceSubset, b::SpaceSubset)::SpaceSubset = __intersect(a, b)

Base.:(==)(a::AbstractSpaceSubset, b::AbstractSpaceSubset)::Bool = (==)(promote(a, b) ...)
Base.:union(a::AbstractSpaceSubset, b::AbstractSpaceSubset)::SpaceSubset = union(promote(a, b) ...)
Base.:intersect(a::AbstractSpaceSubset, b::AbstractSpaceSubset)::SpaceSubset = intersect(promote(a, b) ...)
Base.:isempty(a::AbstractSpaceSubset)::Bool = isempty(convert(SpaceSubset, a).representation)

function __init__()
end

end
