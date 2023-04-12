include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")

using LinearAlgebra, PlotlyJS
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')

k_space = convert(MomentumSpace, triangular)

unit_cell = union(Point([1/3, 2/3], triangular), Point([2/3, 1/3], triangular))
crystal = Crystal(unit_cell, [32, 32])

modes::Subset{Mode} = quantize("physical", unit_cell, 1)
fock_space::FockSpace = FockSpace(modes)

m0, m1 = members(modes)
tₙ = -1.
bonds::Set{Bond} = Set([
    Bond((m0, m1), Point([0, 0], triangular), tₙ),
    Bond((m0, m1), Point([-1, 0], triangular), tₙ),
    Bond((m0, m1), Point([0, 1], triangular), tₙ)
])

sampling = interpolate(Point([-1/2, 0], k_space), Point([1/2, 0], k_space), 300)


# visualize_region("Honeycomb lattice", real_zone, euclidean(RealSpace, 2))
# visualize_region("Honeycomb lattice", reciprocal_zone, euclidean(RealSpace, 2))
