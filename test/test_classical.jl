include("../src/spaces.jl")

using .spaces: RealSpace, MomentumSpace, Point, AbstractSpaceSubset, SpaceSubset

space = RealSpace([1. 0. 1.; 0. 1. 0.; 1. 0. -1.]')

point = Point([1//1, 1//1, 0//1], space)

new_point = Point([1//1, 1//1, 0//1], space) + Point([1//1, 0//1, 0//1], space)

point == convert(SpaceSubset, point)

union(point, convert(SpaceSubset, point))