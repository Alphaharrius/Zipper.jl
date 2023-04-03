module Spaces

using LinearAlgebra

abstract type Element{R <: Any} end

Base.:convert(::Type{T}, source::T) where {T <: Element} = source

representation(source::Element{R}) where {R <: Any} = convert(R, source)

abstract type AbstractSpace{T} <: Element{T} end

abstract type AffineSpace <: AbstractSpace{Matrix{Float64}} end

Base.:convert(::Type{Matrix{Float64}}, source::T) where {T <: AffineSpace} = source.rep

euclidean(S::Type{<: AffineSpace}, n::Int64) = S(Matrix{Float64}(I(n)))

basis(space::AffineSpace)::Matrix{Float64} = representation(space)

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

position(point::Point)::Vector{Rational{Int64}} = point.pos

linear_transform(new_space::AffineSpace, point::Point)::Point = Point(collect(Rational{Int64}, basis(new_space) * inv(basis(space_of(point))) * position(point)), new_space)

Base.:(==)(a::Point, b::Point)::Bool = space_of(a) == space_of(b) && position(a) == position(b)
LinearAlgebra.:norm(point::Point) = norm(position(point))

function Base.:+(a::Point, b::Point)::Point
    @assert(typeof(space_of(a)) == typeof(space_of(b)))
    return Point(position(a) + position(b), space_of(a))
end

function Base.:-(a::Point, b::Point)::Point
    @assert(typeof(space_of(a)) == typeof(space_of(b)))
    return Point(position(a) - position(b), space_of(a))
end

distance(a::Point, b::Point)::Float64 = sqrt(norm(a - b))

struct Subset <: AbstractSubset{AbstractSubset}
    rep::Set{AbstractSubset}
    space::AbstractSpace
end

function Base.:union(input::Subset...)::Subset
    @assert(length(Set{AbstractSpace}([space_of(subset) for subset in input])) == 1)
    return Subset(union([representation(subset) for subset in input]...), space_of(input[1]))
end

function Base.:intersect(input::Subset...)::Subset
    @assert(length(Set{AbstractSpace}([space_of(subset) for subset in input])) == 1)
    return Subset(intersect([representation(subset) for subset in input]...), space_of(input[1]))
end

Base.:convert(::Type{Subset}, source::Point) = Subset(Set{Point}([source]), space_of(source))
Base.:convert(::Type{Set{AbstractSubset}}, source::Subset) = source.rep

export Element, AbstractSpace, AffineSpace, RealSpace, MomentumSpace, AbstractSubset, Point
export representation, euclidean, basis, dimension, space_of, center, position, linear_transform, distance

end
