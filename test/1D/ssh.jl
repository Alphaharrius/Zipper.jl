include("../../src/spaces.jl")
include("../../src/geometries.jl")
include("../../src/quantum.jl")
include("../../src/physical.jl")
include("../../src/plotting.jl")

using LinearAlgebra, PlotlyJS
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting

chain = euclidean(RealSpace, 1)
unit_cell = union(Point([1/4], chain), Point([3/4], chain))
ssh_chain = Crystal(unit_cell, [64])

k_space = convert(MomentumSpace, chain)

modes::Subset{Mode} = quantize("physical", unit_cell, 1)
m0, m1 = members(modes)
t_n = -0.4
t_nn = -0.6
bonds::Set{Bond} = Set([
    Bond((m0, m1), Point([0], chain), t_n),
    Bond((m1, m0), Point([1], chain), t_nn)
])
sampling = interpolate(Point([-4], k_space), Point([4], k_space), 1000)
spectrum = real(hcat([eigvals(test(bonds, k, ssh_chain)) for k in sampling]...))
top = spectrum[1, :]
bottom = spectrum[2, :]

trace0 = scatter(y=top)
trace1 = scatter(y=bottom)
plot([trace0, trace1])
