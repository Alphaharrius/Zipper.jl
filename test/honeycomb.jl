include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")

using LinearAlgebra, PlotlyJS, BenchmarkTools
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')

k_space = convert(MomentumSpace, triangular)

unit_cell = union(Point([1/3, 2/3], triangular), Point([2/3, 1/3], triangular))
crystal = Crystal(unit_cell, [32, 32])

modes::Subset{Mode} = quantize("physical", :pos, unit_cell, 1)
fock_space::FockSpace = FockSpace(modes)

m0, m1 = members(modes)

k::Point = Point([1/5, 1/3], k_space)
k1::Point = Point([1/7, 1/8], k_space)
subset::Subset{Mode} = Subset(map(m -> identify(m, :offset),[m0, m1, setattr(m1, :offset => Point([0, -1], triangular))]))
kk = brillouin_zone(crystal)
Fₖ::FockMap = fourier(Subset([k]), subset)

momentums = kk
inmodes = subset
dict::Dict{Tuple{Mode, Mode}, ComplexF64} = Dict()
# Disable the identification by :offset so that all inmodes collapes to the basis mode. 
basismodes::Set{Mode} = Set(unidentify(inmode, :offset) for inmode in inmodes)
# Enable the identification by :offset when adding momentum as :offset for each basis mode.
outmodes = Iterators.map(tup -> identify(setattr(tup[1], :offset => tup[2]), :offset), Iterators.product(basismodes, momentums))
@time for (outmode, inmode) in Iterators.product(outmodes, inmodes)
    momentum::Point = getattr(outmode, :offset)
    euc_momentum::Point = linear_transform(euclidean(MomentumSpace, length(momentum)), momentum)
    inoffset::Point = getattr(inmode, :offset)
    euc_inoffset::Point = linear_transform(euclidean(RealSpace, length(inoffset)), inoffset)
    dict[(outmode, inmode)] = exp(-1im * dot(euc_momentum, euc_inoffset))
end

@time outspace = FockSpace(Subset(Set(outmodes)))
@time inspace = FockSpace(inmodes)

@time rep::SparseMatrixCSC{ComplexF64, Int64} = spzeros(dimension(outspace), dimension(inspace))
@time for ((out_mode::Mode, in_mode::Mode), value::ComplexF64) in dict
    rep[1, inspace.ordering[in_mode]] = value
end

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
