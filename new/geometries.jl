include("spaces.jl")

module Geometries

using ..Spaces: Subset, Point, distance, representation

radius(region::Subset, center::Point)::Float64 = maximum(distance(center, point) for point in representation(region))

export radius

end