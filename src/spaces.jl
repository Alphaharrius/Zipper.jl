module Spaces

using LinearAlgebra

abstract type Element{R <: Any} end

""" Reflexive relation. """
Base.:convert(::Type{T}, source::T) where {T <: Element} = source
Base.:convert(::Type{Set{T}}, source::T) where {T <: Element} = Set{T}([source])

rep(source::Element{R}) where {R <: Any} = convert(R, source)

abstract type AbstractSpace{T} <: Element{T} end

dimension(space::T) where {T <: AbstractSpace} = length(space)

abstract type AffineSpace <: AbstractSpace{Matrix{Float64}} end

Base.:convert(::Type{Matrix{Float64}}, source::AffineSpace) = source.rep

euclidean(S::Type{<: AffineSpace}, n::Int64) = S(Matrix{Float64}(I(n)))

basis(space::AffineSpace)::Matrix{Float64} = rep(space)

Base.:length(space::T) where {T <: AffineSpace} = size(basis(space), 1)

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

# Since it is not possible to hash the AffineSpace with Matrix{Float64}, we will hash only the Vector representation and leave to the Base.isequal for identification.
Base.:hash(point::Point)::UInt = hash(pos(point))
Base.:isequal(a::Point, b::Point)::Bool = a == b

center(space::AffineSpace)::Point = Point(zeros(Rational{Int64}, dimension(space)), space)

pos(point::Point)::Vector{Rational{Int64}} = point.pos

linear_transform(new_space::AffineSpace, point::Point)::Point = Point(collect(Rational{Int64}, inv(basis(new_space)) * basis(space_of(point)) * pos(point)), new_space)

function fourier_coef(momentum::Point, point::Point, vol::Integer)::ComplexF64
    m_space::MomentumSpace = space_of(momentum)
    r_space::RealSpace = space_of(point)
    phys_momentum::Point = linear_transform(euclidean(MomentumSpace, dimension(m_space)), momentum)
    phys_point::Point = linear_transform(euclidean(RealSpace, dimension(r_space)), point)
    return exp(-1im * dot(collect(Float64, pos(phys_momentum)), collect(Float64, pos(phys_point)))) / sqrt(vol)
end

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

Base.:*(point::Point, val::T) where {T <: Real} = Point(collect(Float64, pos(point)) * val, space_of(point))
Base.:/(point::Point, val::T) where {T <: Real} = Point(collect(Float64, pos(point)) / val, space_of(point))

distance(a::Point, b::Point)::Float64 = sqrt(norm(a - b))

function interpolate(from::Point, to::Point, count::T)::Array{Point} where {T <: Integer}
    @assert(space_of(from) == space_of(to))
    march::Point = (to - from) / (count + 1)
    points::Array{Point} = [from + march * n for n in 0:(count + 1)]
    points[end] = to
    return points
end

struct Subset{T <: AbstractSubset} <: AbstractSubset{T}
    rep::Set{T}
    space::AbstractSpace

    Subset(elements::Set{T}) where {T <: AbstractSubset} = new{T}(elements, space_of(first(elements)))
end

# Since it is not possible to hash the AffineSpace with Matrix{Float64}, we will hash only the Set representation and leave to the Base.isequal for identification.
Base.:hash(subset::Subset)::UInt = hash(rep(subset))
Base.:isequal(a::Subset, b::Subset)::Bool = a == b

Base.:iterate(subset::T, i...) where {T <: Subset} = Base.iterate(rep(subset), i...)
Base.:length(subset::T) where {T <: Subset} = Base.length(rep(subset))

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
    Subset(Set{A}([convert(A, el) for el in rep(source)])))
""" Convert to the generalized type of `AbstractSubset`. """
Base.:convert(::Type{Subset}, source::T) where {T <: AbstractSubset} = Subset(rep(source))
Base.:convert(::Type{Subset{T}}, source::T) where {T <: AbstractSubset} = convert(Subset, source)

function Base.:union(input::Subset...)::Subset
    @assert(length(Set{AbstractSpace}([space_of(subset) for subset in input])) == 1)
    return Subset(union([rep(subset) for subset in input]...))
end

Base.:union(input::T...) where {T <: AbstractSubset} = union((convert(Subset{T}, element) for element in input)...)

function Base.:intersect(input::Subset...)::Subset
    @assert(length(Set{AbstractSpace}([space_of(subset) for subset in input])) == 1)
    return Subset(intersect([rep(subset) for subset in input]...))
end

Base.:intersect(input::T...) where {T <: AbstractSubset} = intersect((convert(Subset{T}, element) for element in input)...)

export Element, AbstractSpace, AffineSpace, RealSpace, MomentumSpace, AbstractSubset, Point, Subset
export rep, euclidean, basis, dimension, space_of, center, pos, linear_transform, fourier_coef, distance, interpolate, flatten, members

end
