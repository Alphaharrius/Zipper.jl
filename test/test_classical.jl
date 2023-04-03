include("../new/geometries.jl")

using ..Spaces: RealSpace, Subset, Point, center
using ..Geometries: radius

space = RealSpace([1. 1.5; 0.5 0.5])
point0 = Point([1//1, 3//1], space)
point1 = Point([-8//1, 3//1], space)
point2 = Point([0//1, 0//1], space)

region = union(convert(Subset, point0), convert(Subset, point1), convert(Subset, point2))

radius(region, center(space))
