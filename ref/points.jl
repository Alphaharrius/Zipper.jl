module Points
export map_to, center, abs, distance
export Point

include("spaces.jl")

using LinearAlgebra
using ..Spaces

"""
    Point(position, parent_space)

Create a `Point` of `position` in the `Space` of `parent_space`.
"""
struct Point <: AbstractSpaceSubset
    """ The position in spacial basis in the `parent_space`. """
    position::Vector{Rational{Int64}} # Julia can compare rational numbers.
    """ The parent space that includes this point. """
    parent_space::ClassicalSpace

    function Point(position::Vector{Rational{Int64}}, parent_space::ClassicalSpace)::Point
        @assert(length(position) == dimension(parent_space))
        return new(position, parent_space)
    end
end

function Base.:+(a::Point, b::Point)::Point
    @assert(typeof(a.parent_space) == typeof(b.parent_space))
    return Point(a.position + b.position, a.parent_space)
end

function Base.:-(a::Point, b::Point)::Point
    @assert(typeof(a.parent_space) == typeof(b.parent_space))
    return Point(a.position - b.position, a.parent_space)
end

map_to(point::Point, space::ClassicalSpace)::Point = Point(collect(Rational{Int64}, basis(space) * inv(basis(point.parent_space)) * point.position), space)
center(space::ClassicalSpace)::Point = Point(zeros(Rational, dimension(space)), space)
abs(point::Point)::Float64 = norm(point.position)
distance(from::Point, to::Point) = sqrt(abs(to - from))

Base.:(==)(a::Point, b::Point)::Bool = a.parent_space == b.parent_space && a.position == b.position

Base.:convert(::Type{SpaceSubset}, source::Point) = SpaceSubset(Set{AbstractSpaceSubset}([source]), source.parent_space)

end