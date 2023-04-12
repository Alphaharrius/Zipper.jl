include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")

using LinearAlgebra, PlotlyJS
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting

fcc_space = RealSpace([1/2 1/2 0; 0. 1/2 1/2; 1/2 0. 1/2]')

r_space = convert(MomentumSpace, fcc_space)
unit_cell = union(Point([0, 0, 0], fcc_space), Point([1/4, 1/4, 1/4], fcc_space))
crystal = Crystal(unit_cell, [4, 4, 4])
zone::Subset{Point} = points(crystal)
k_zone::Subset{Point} = brillouin_zone(crystal)

modes::Subset{Mode} = quantize("physical", unit_cell, 1)
fock::FockSpace = FockSpace(modes)

m0, m1 = members(modes)
tₙ = -1.
bonds::Set{Bond} = Set([
    Bond((m0, m1), Point([0, 0, 0], fcc_space), tₙ),
    Bond((m0, m1), Point([1, 0, 0], fcc_space), tₙ),
    Bond((m0, m1), Point([0, 1, 0], fcc_space), tₙ),
    Bond((m0, m1), Point([0, 0, 1], fcc_space), tₙ)
])

ΓX = interpolate(Point([0, 0, 0], r_space), Point([1/2, 0, 0], r_space), 200)
XW = interpolate(Point([1/2, 0, 0], r_space), Point([1/2, 1/4, 0], r_space), 200)
WL = interpolate(Point([1/2, 1/4, 0], r_space), Point([1/4, 1/4, 1/4], r_space), 200)
LΓ = interpolate(Point([1/4, 1/4, 1/4], r_space), Point([0, 0, 0], r_space), 200)
ΓK = interpolate(Point([0, 0, 0], r_space), Point([3/8, 3/8, 0], r_space), 200)
KW = interpolate(Point([3/8, 3/8, 0], r_space), Point([1/2, 1/4, 0], r_space), 200)
WU = interpolate(Point([1/2, 1/4, 0], r_space), Point([1/2, 1/8, 1/8], r_space), 200)
UX = interpolate(Point([1/2, 1/8, 1/8], r_space), Point([1/2, 0, 0], r_space), 200)
line = vcat(ΓX, XW, WL, LΓ, ΓK, KW, WU, UX)
spectrum = hcat([eigvalsh(bloch(bonds, k, crystal)) for k in line]...)
top = map(c -> c.value, spectrum[1, :])
bottom = map(c -> c.value, spectrum[2, :])
plot([scatter(y=top), scatter(y=bottom)])

visualize_region("Honeycomb lattice", zone, euclidean(RealSpace, 3))
visualize_region("Honeycomb lattice", k_zone, euclidean(MomentumSpace, 3))
