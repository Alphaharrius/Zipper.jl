include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")

using LinearAlgebra, PlotlyJS
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')

k_space = convert(MomentumSpace, triangular)

unitcell = union(Point([1/3, 2/3], triangular), Point([2/3, 1/3], triangular))
crystal = Crystal(unitcell, [128, 128])
modes::Subset{Mode} = quantize("physical", :pos, unitcell, 1)
m0, m1 = members(modes)

tₙ = ComplexF64(-1.)
bm = bondmap([
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => tₙ])

rng = -1:0.1:1
tspec::Matrix{Float64} = zeros(Float64, length(rng), length(rng))
bspec::Matrix{Float64} = zeros(Float64, length(rng), length(rng))
for (j, k) in enumerate(rng)
    for (i, h) in enumerate(rng)
        fkk::FockMap = fourier(Subset([Point([h, k], k_space)]), Subset(ordered_modes(bm.outspace)))
        bloch::FockMap = fkk * bm * dagger(fkk)
        spec = eigvalsh(bloch)
        tspec[j, i] = real(spec[1].second)
        bspec[j, i] = real(spec[2].second)
    end
end

# spec = [hcat([eigvalsh(bloch(bonds, Point([h, k], k_space), crystal), :offset => Point([h, k], k_space)) for h in -1:0.1:1]...) for k in -1:0.1:1]

# top = map(p -> linear_transform(euclidean(MomentumSpace, 2), getattr(p.first, :offset)) => p.second, hcat([v[1, :] for v in spec]...))
# bottom = map(p -> linear_transform(euclidean(MomentumSpace, 2), getattr(p.first, :offset)) => p.second, hcat([v[2, :] for v in spec]...))

# topvals = hcat([map(p -> [pos(p.first)..., p.second], top)...]...)
# botvals = hcat([map(p -> [pos(p.first)..., p.second], bottom)...]...)

plot([surface(z=tspec), surface(z=bspec)])



# visualize_region("Honeycomb lattice", real_zone, euclidean(RealSpace, 2))
# visualize_region("Honeycomb lattice", reciprocal_zone, euclidean(RealSpace, 2))
