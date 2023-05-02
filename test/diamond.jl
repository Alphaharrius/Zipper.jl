include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")

using LinearAlgebra, PlotlyJS
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting

fcc = RealSpace([1/2 1/2 0; 0. 1/2 1/2; 1/2 0. 1/2]')

kspace = convert(MomentumSpace, fcc)
unitcell = union(Point([0, 0, 0], fcc), Point([1/4, 1/4, 1/4], fcc))
crystal = Crystal(unitcell, [4, 4, 4])
zone::Subset{Point} = points(crystal)
k_zone::Subset{Point} = brillouinzone(crystal)

modes::Subset{Mode} = quantize("physical", :pos, unitcell, 1)
fock::FockSpace = FockSpace(modes)

m0, m1 = members(modes)
tₙ = ComplexF64(-1.)

bmap::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([1, 0, 0], fcc))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1, 0], fcc))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 0, 1], fcc))) => tₙ])

ΓX = interpolate(Point([0, 0, 0], kspace), Point([1/2, 0, 0], kspace), 200)
XW = interpolate(Point([1/2, 0, 0], kspace), Point([1/2, 1/4, 0], kspace), 200)
WL = interpolate(Point([1/2, 1/4, 0], kspace), Point([1/4, 1/4, 1/4], kspace), 200)
LΓ = interpolate(Point([1/4, 1/4, 1/4], kspace), Point([0, 0, 0], kspace), 200)
ΓK = interpolate(Point([0, 0, 0], kspace), Point([3/8, 3/8, 0], kspace), 200)
KW = interpolate(Point([3/8, 3/8, 0], kspace), Point([1/2, 1/4, 0], kspace), 200)
WU = interpolate(Point([1/2, 1/4, 0], kspace), Point([1/2, 1/8, 1/8], kspace), 200)
UX = interpolate(Point([1/2, 1/8, 1/8], kspace), Point([1/2, 0, 0], kspace), 200)
line = vcat(ΓX, XW, WL, LΓ, ΓK, KW, WU, UX)
spectrum = espec(bmap, line)

tspec = filter(p -> getattr(p.first, :index) == 1, spectrum)
bspec = filter(p -> getattr(p.first, :index) == 2, spectrum)

trace0 = scatter(y=map(p -> p.second, tspec))
trace1 = scatter(y=map(p -> p.second, bspec))
plot([trace0, trace1])

rng = -2:0.1:2
tspecs::Matrix{Float64} = zeros(Float64, length(rng), length(rng))
bspecs::Matrix{Float64} = zeros(Float64, length(rng), length(rng))
for (j, k) in enumerate(rng)
    for (i, h) in enumerate(rng)
        fkk::FockMap = fourier(Subset([Point([h, k, 0], kspace)]), Subset(orderedmodes(bmap.outspace)))
        bloch::FockMap = fkk * bmap * fkk'
        spec = eigvalsh(bloch)
        tspecs[j, i] = real(spec[1].second)
        bspecs[j, i] = real(spec[2].second)
    end
end

plot([surface(z=tspecs), surface(z=bspecs)])
