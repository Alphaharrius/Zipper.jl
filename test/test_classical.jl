include("../src/classical.jl")

using .classical: ClassicalSpace, Point

space = ClassicalSpace([1. 0. 1.; 0. 1. 0.; 1. 0. -1.]')

point = Point([1//1, 1//1, 0//1], space)

new_point = Point([1//1, 1//1, 0//1], space) + Point([1//1, 0//1, 0//1], space)
