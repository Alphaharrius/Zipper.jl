include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")

using ..Spaces
using ..Geometries
using ..Quantum
using ..Physical
using ..Plotting

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')

momentum_space = convert(MomentumSpace, triangular)

honeycomb_unit_cell = union(Point([1/3, 2/3], triangular), Point([2/3, 1/3], triangular))
honeycomb_crystal = Crystal(honeycomb_unit_cell, [256, 256])
real_zone::Subset{Point} = points(honeycomb_crystal)
reciprocal_zone::Subset{Point} = brillouin_zone(honeycomb_crystal)

modes::Subset{Mode} = quantize("physical", honeycomb_unit_cell, 1)
fock_space::FockSpace = FockSpace(modes)

m0, m1 = members(modes)

visualize_region("Honeycomb lattice", real_zone, euclidean(RealSpace, 2))
visualize_region("Honeycomb lattice", reciprocal_zone, euclidean(RealSpace, 2))
