module Geometries
export radius, span_subset

include("spaces.jl")
include("points.jl")

using ..Spaces, ..Points

function radius(center::Point, region::SpaceSubset)::Float64
    @assert(!isempty(intersect(center, region)))
    return maximum(distance(center, point) for point in region.representation)
end

function span_subset(dimensions::Tuple{Int64, Vararg{Int64}}, space::ClassicalSpace)
    raw_positions::Array{Tuple{Int64, Vararg{Int64}}} = collect(Iterators.product((0:(n - 1) for n in dimensions)...))
    representation::Set{Point} = Set{Point}([Point(collect(Rational{Int64}, raw), space) for raw in raw_positions])
    return SpaceSubset(representation, space)
end

end