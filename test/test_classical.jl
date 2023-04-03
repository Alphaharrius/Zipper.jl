include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/plotting.jl")

using ..Spaces
using ..Geometries
using ..Quantum
using ..Plotting

space = RealSpace([1. 1.5; 0.5 0.5])
point0 = Point([1//1, 3//1], space)
point1 = Point([-8//1, 3//1], space)
point2 = Point([0//1, 0//1], space)

region::Subset{Point} = union(point0, point1, point2)

visualize_region(region, euclidean(RealSpace, 2))

modes::Subset{Mode} = quantize(region, 2)
m0, m1, m2 = members(modes)
fock_space::FockSpace = FockSpace(modes)

fock_map::FockMap = FockMap(
    Dict((m0, m0) => 1. + 0im,
         (m1, m1) => 1. + 0im,
         (m2, m2) => 1. + 0im),
    fock_space, fock_space)

cut::FockMap = columns(fock_map, union(m0, m1))
