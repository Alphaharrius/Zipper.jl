include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")

using LinearAlgebra, PlotlyJS, BenchmarkTools
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')

k_space = convert(MomentumSpace, triangular)

unitcell = union(Point([1/3, 2/3], triangular), Point([2/3, 1/3], triangular))
crystal = Crystal(unitcell, [128, 128])
modes::Subset{Mode} = quantize("physical", :pos, unitcell, 1)
m0, m1 = members(modes)

subset::Subset{Mode} = Subset([m0, m1, setattr(m1, :offset => Point([0, -1], triangular))])
bzone = @time brillouin_zone(crystal)

Fₖ::FockMap = @time fourier(bzone, subset)

tₙ = -1.
bonds::Set{Bond} = Set([
    Bond((m0, m1), Point([0, 0], triangular), tₙ),
    Bond((m0, m1), Point([-1, 0], triangular), tₙ),
    Bond((m0, m1), Point([0, 1], triangular), tₙ)
])

spec = [hcat([eigvalsh(bloch(bonds, Point([h, k], k_space), crystal), :offset => Point([h, k], k_space)) for h in -1:0.1:1]...) for k in -1:0.1:1]

top = map(p -> linear_transform(euclidean(MomentumSpace, 2), getattr(p.first, :offset)) => p.second, hcat([v[1, :] for v in spec]...))
bottom = map(p -> linear_transform(euclidean(MomentumSpace, 2), getattr(p.first, :offset)) => p.second, hcat([v[2, :] for v in spec]...))

topvals = hcat([map(p -> [pos(p.first)..., p.second], top)...]...)
botvals = hcat([map(p -> [pos(p.first)..., p.second], bottom)...]...)

layout = Layout(title="", scene=attr(aspectmode="data"))
plot([scatter3d(x=topvals[1,:], y=topvals[2,:], z=topvals[3,:]*100), scatter3d(x=botvals[1,:], y=botvals[2,:], z=botvals[3,:]*100)], layout)



# visualize_region("Honeycomb lattice", real_zone, euclidean(RealSpace, 2))
# visualize_region("Honeycomb lattice", reciprocal_zone, euclidean(RealSpace, 2))
