include("../src/classical.jl")

using .classical: ClassicalSpace, Point, create_classical_space, create_point, phys_operator

space = create_classical_space([1. 0. 1.; 0. 1. 0.; 1. 0. -1.]')
new_point = create_point([1., 1., 0.], space) + create_point([1., 0., 0.], space)
phys_operator(space)(new_point)
