include(joinpath("..", "src", "classical.jl"))

using .classical: ClassicalSpace, Point, classical_space, point

space = classical_space([1. 0. 1.; 0. 1. 1.]')
new_point = point([1., 1., 0.], space) + point([1., 0., 0.], space)
