include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/plotting.jl")

using ..Spaces
using ..Geometries
using ..Quantum
using ..Plotting

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
honeycomb_unit_cell = union(Point([1/3, 2/3], triangular), Point([2/3, 1/3], triangular))
honeycomb_crystal = Crystal(honeycomb_unit_cell, [32, 32])
real_zone::Subset{Point} = points(honeycomb_crystal)
reciprocal_zone::Subset{Point} = brillouin_zone(honeycomb_crystal)

visualize_region("Honeycomb lattice", real_zone, euclidean(RealSpace, 2))
visualize_region("Honeycomb lattice", reciprocal_zone, euclidean(RealSpace, 2))

space = RealSpace([1. 1.5; 0.5 0.5])
point0 = Point([1//1, 3//1], space)
point1 = Point([-8//1, 3//1], space)
point2 = Point([0//1, 0//1], space)

region::Subset{Point} = union(point0, point1, point2)

# visualize_region(region, euclidean(RealSpace, 2))

modes::Subset{Mode} = quantize("physical", region, 2)
m0, m1, m2 = members(modes)
fock_space::FockSpace = FockSpace(modes)

fock_map::FockMap = FockMap(
    Dict((m0, m0) => 1. + 0im,
         (m1, m1) => 1. + 0im,
         (m2, m2) => 1. + 0im),
    fock_space, fock_space)

cut::FockMap = columns(fock_map, union(m0))
