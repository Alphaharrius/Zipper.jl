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

modes::Subset{Mode} = quantize("physical", :pos, unit_cell, 1)
m0, m1 = members(modes)
t_n = ComplexF64(-0.4)
t_nn = ComplexF64(-0.6)

bonds::FockMap = bondmap([
    (m0, m1) => t_n,
    (m1, setattr(m0, :offset => Point([1], chain))) => t_nn])

sampling::Vector{Point} = interpolate(Point([-4], k_space), Point([4], k_space), 1000)
spectrum = espec(bonds, sampling)

tspec = filter(p -> getattr(p.first, :index) == 1, spectrum)
bspec = filter(p -> getattr(p.first, :index) == 2, spectrum)

trace0 = scatter(y=map(p -> p.second, tspec))
trace1 = scatter(y=map(p -> p.second, bspec))
plot([trace0, trace1])
