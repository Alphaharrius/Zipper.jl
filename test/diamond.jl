include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")

using LinearAlgebra, PlotlyJS
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting

fcc = RealSpace([1/2 1/2 0; 0. 1/2 1/2; 1/2 0. 1/2]')

r_space = convert(MomentumSpace, fcc)
unit_cell = union(Point([0, 0, 0], fcc), Point([1/4, 1/4, 1/4], fcc))
crystal = Crystal(unit_cell, [4, 4, 4])
zone::Subset{Point} = points(crystal)
k_zone::Subset{Point} = brillouin_zone(crystal)

modes::Subset{Mode} = quantize("physical", :pos, unit_cell, 1)
fock::FockSpace = FockSpace(modes)

m0, m1 = members(modes)
tₙ = ComplexF64(-1.)
bmap::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([1, 0, 0], fcc))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1, 0], fcc))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 0, 1], fcc))) => tₙ])

ΓX = interpolate(Point([0, 0, 0], r_space), Point([1/2, 0, 0], r_space), 200)
XW = interpolate(Point([1/2, 0, 0], r_space), Point([1/2, 1/4, 0], r_space), 200)
WL = interpolate(Point([1/2, 1/4, 0], r_space), Point([1/4, 1/4, 1/4], r_space), 200)
LΓ = interpolate(Point([1/4, 1/4, 1/4], r_space), Point([0, 0, 0], r_space), 200)
ΓK = interpolate(Point([0, 0, 0], r_space), Point([3/8, 3/8, 0], r_space), 200)
KW = interpolate(Point([3/8, 3/8, 0], r_space), Point([1/2, 1/4, 0], r_space), 200)
WU = interpolate(Point([1/2, 1/4, 0], r_space), Point([1/2, 1/8, 1/8], r_space), 200)
UX = interpolate(Point([1/2, 1/8, 1/8], r_space), Point([1/2, 0, 0], r_space), 200)
line = vcat(ΓX, XW, WL, LΓ, ΓK, KW, WU, UX)
spectrum = espec(bmap, line)

tspec = filter(p -> getattr(p.first, :index) == 1, spectrum)
bspec = filter(p -> getattr(p.first, :index) == 2, spectrum)

trace0 = scatter(y=map(p -> p.second, tspec))
trace1 = scatter(y=map(p -> p.second, bspec))
plot([trace0, trace1])
