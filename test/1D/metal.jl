include("../../src/spaces.jl")
include("../../src/geometries.jl")
include("../../src/quantum.jl")
include("../../src/physical.jl")
include("../../src/plotting.jl")

using LinearAlgebra, PlotlyJS
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting

real_space = euclidean(RealSpace, 1)
unit_cell = union(Point([1/2], real_space))
chain = Crystal(unit_cell, [256])

k_space = convert(MomentumSpace, real_space)

modes::Subset{Mode} = quantize("physical", unit_cell, 1)
m0 = first(rep(modes))
t_n = -1
bonds::Set{Bond} = Set([
    Bond((m0, m0), Point([1], real_space), t_n)
])
sampling = interpolate(Point([-4], k_space), Point([4], k_space), 1000)
spectrum = hcat([eigvalsh(bloch(bonds, k, chain)) for k in sampling]...)
top = map(p -> p.second, spectrum[1, :])

trace0 = scatter(y=top)
plot([trace0])
