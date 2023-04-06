module Spaces

using LinearAlgebra

abstract type Element{R <: Any} end

""" Reflexive relation. """
Base.:convert(::Type{T}, source::T) where {T <: Element} = source
Base.:convert(::Type{Set{T}}, source::T) where {T <: Element} = Set{T}([source])

rep(source::Element{R}) where {R <: Any} = convert(R, source)

abstract type AbstractSpace{T} <: Element{T} end

abstract type AffineSpace <: AbstractSpace{Matrix{Float64}} end

Base.:convert(::Type{Matrix{Float64}}, source::AffineSpace) = source.rep

euclidean(S::Type{<: AffineSpace}, n::Int64) = S(Matrix{Float64}(I(n)))

basis(space::AffineSpace)::Matrix{Float64} = rep(space)

dimension(space::AffineSpace) = size(basis(space), 1)

abstract type AbstractSubset{T} <: Element{Set{T}} end

space_of(subset::AbstractSubset) = subset.space

struct RealSpace <: AffineSpace
    rep::Matrix{Float64}
end

struct MomentumSpace <: AffineSpace
    rep::Matrix{Float64}
end

Base.:convert(::Type{MomentumSpace}, source::RealSpace)::MomentumSpace = MomentumSpace(2.0 * π * Matrix(transpose(inv(basis(source)))))
Base.:convert(::Type{RealSpace}, source::MomentumSpace)::RealSpace = RealSpace(Matrix(transpose(inv(basis(source) / (2.0 * π)))))

struct Point <: AbstractSubset{Point}
    pos::Vector{Rational{Int64}}
    space::AbstractSpace
end

center(space::AffineSpace)::Point = Point(zeros(Rational{Int64}, dimension(space)), space)

pos(point::Point)::Vector{Rational{Int64}} = point.pos

linear_transform(new_space::AffineSpace, point::Point)::Point = Point(collect(Rational{Int64}, inv(basis(new_space)) * basis(space_of(point)) * pos(point)), new_space)

Base.:(==)(a::Point, b::Point)::Bool = space_of(a) == space_of(b) && pos(a) == pos(b)
LinearAlgebra.:norm(point::Point) = norm(pos(point))

function Base.:+(a::Point, b::Point)::Point
    @assert(typeof(space_of(a)) == typeof(space_of(b)))
    return Point(pos(a) + pos(b), space_of(a))
end

function Base.:-(a::Point, b::Point)::Point
    @assert(typeof(space_of(a)) == typeof(space_of(b)))
    return Point(pos(a) - pos(b), space_of(a))
end

distance(a::Point, b::Point)::Float64 = sqrt(norm(a - b))

struct Subset{T <: AbstractSubset} <: AbstractSubset{T}
    rep::Set{T}
    space::AbstractSpace
end

function flatten(subset::Subset{T})::Subset where {T <: AbstractSubset}
    if T <: Subset
        return union([flatten(element) for element in rep(subset)]...)
    end
    return subset
end

members(subset::Subset)::Tuple = (rep(subset)...,)

Base.:(==)(a::Subset, b::Subset)::Bool = space_of(a) == space_of(b) && rep(a) == rep(b)

Base.:convert(::Type{Set{T}}, source::Subset{T}) where {T <: AbstractSubset} = source.rep
""" Reflexive relation. """
Base.:convert(::Type{Subset}, source::Subset{T}) where {T <: AbstractSubset} = source
Base.:convert(::Type{Subset{T}}, source::Subset{T}) where {T <: AbstractSubset} = source
""" Element-wise conversions. """
Base.:convert(::Type{Subset{A}}, source::Subset{B}) where {A <: AbstractSubset, B <: AbstractSubset} = (
    Subset(Set{A}([convert(A, el) for el in rep(source)]), space_of(source)))
""" Convert to the generalized type of `AbstractSubset`. """
Base.:convert(::Type{Subset}, source::T) where {T <: AbstractSubset} = Subset{T}(rep(source), space_of(source))
Base.:convert(::Type{Subset{T}}, source::T) where {T <: AbstractSubset} = convert(Subset, source)

function Base.:union(input::Subset...)::Subset
    @assert(length(Set{AbstractSpace}([space_of(subset) for subset in input])) == 1)
    return Subset(union([rep(subset) for subset in input]...), space_of(first(input)))
end

Base.:union(input::T...) where {T <: AbstractSubset} = union((convert(Subset{T}, element) for element in input)...)

function Base.:intersect(input::Subset...)::Subset
    @assert(length(Set{AbstractSpace}([space_of(subset) for subset in input])) == 1)
    return Subset(intersect([rep(subset) for subset in input]...), space_of(first(input)))
end

Base.:intersect(input::T...) where {T <: AbstractSubset} = intersect((convert(Subset{T}, element) for element in input)...)

export Element, AbstractSpace, AffineSpace, RealSpace, MomentumSpace, AbstractSubset, Point, Subset
export rep, euclidean, basis, dimension, space_of, center, pos, linear_transform, distance, flatten, members

end
