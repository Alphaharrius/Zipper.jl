if !isdefined(Main, :Spaces) include("spaces.jl") end

module Geometries

using ..Spaces

radius(region::Subset, center::Point)::Float64 = maximum(distance(center, point) for point in representation(region))

export radius

end