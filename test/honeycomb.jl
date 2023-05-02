include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")
include("../src/zer.jl")

using LinearAlgebra, PlotlyJS, OrderedCollections, SparseArrays, ColorTypes 
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting, ..Zer

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')

k_space = convert(MomentumSpace, triangular)

unitcell = union(Point([1/3, 2/3], triangular), Point([2/3, 1/3], triangular))
crystal = Crystal(unitcell, [32, 32])
modes::Subset{Mode} = quantize("physical", :pos, unitcell, 1)
m0, m1 = members(modes)
tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => tₙ])

be = eigvecsh(bonds)

𝐻::FockMap = @time hamiltonian(crystal, bonds)
filled::FockMap = @time groundstates(𝐻)
𝐶::FockMap = filled * filled'

plot(heatmap(z=real(rep(𝐶))))

genmodes = [setattr(m, :offset => Point([x, y], triangular)) for (x, y) in Iterators.product(-4:1:4, -4:1:4) for m in modes]
function incirc(mode::Mode)::Bool
    𝑝 = getattr(mode, :offset) + getattr(mode, :pos)
    𝑝ₑ = linear_transform(euclidean(RealSpace, dimension(𝑝)), 𝑝)
    return sqrt(norm(𝑝ₑ)) < 1.3
end
circular = Subset(filter(m -> incirc(m), genmodes))

distill = frozenisometries(crystal, 𝐶, DistillRegion(Point([0, 0], triangular), circular))

function globaldistiller_(crystal::Crystal, correlations::FockMap, region::Zer.DistillRegion; frozen_threshold::Float64 = 1e-3)
    𝐹ₖ::FockMap = fourier(brillouin_zone(crystal), FockSpace(region.modes)) / sqrt(vol(crystal))
    𝐶ᵣ::FockMap = 𝐹ₖ' * correlations * 𝐹ₖ
    corrspec::Vector{Pair{Mode, Float64}} = eigvalsh(𝐶ᵣ)
    plot(scatter(y=map(p -> p.second, corrspec), mode="markers"))
    filledmodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second < frozen_threshold, corrspec)))
    emptymodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second > 1.0 - frozen_threshold, corrspec)))
    𝔘ᵣ::FockMap = eigvecsh(𝐶ᵣ)
    filledisometry::FockMap = columns(𝔘ᵣ, FockSpace(filledmodes))
    emptyisometry::FockMap = columns(𝔘ᵣ, FockSpace(emptymodes))
    # TODO: Handle symmetry.
    return emptyisometry * emptyisometry' - filledisometry * filledisometry'
end

gd = globaldistiller_(crystal, 𝐶, Zer.DistillRegion(Point([0, 0], triangular), circular))
plot(heatmap(z=real(rep(gd))))
gdspec = eigvalsh(gd)
plot(scatter(y=map(p -> p.second, gdspec), mode="markers"))
gdvecs = eigvecsh(gd)
plot(heatmap(z=real(rep(gdvecs))))
gdm = gdspec[24].first
gdmm = columns(gdvecs, FockSpace(Subset([gdm])))
values = columnspec(gdmm)
visualize_spectrum("Mode", values)

cvs = eigvalsh(𝐶)
plot(scatter(y=map(p -> p.second, cvs), mode="markers"))
𝔘::FockMap = eigvecsh(𝐶)
plot(heatmap(z=real(rep(𝔘))))

genmodes = [setattr(m, :offset => Point([x, y], triangular)) for (x, y) in Iterators.product(-4:1:4, -4:1:4) for m in modes]
function incirc(mode::Mode)::Bool
    𝑝 = getattr(mode, :offset) + getattr(mode, :pos)
    𝑝ₑ = linear_transform(euclidean(RealSpace, dimension(𝑝)), 𝑝)
    return sqrt(norm(𝑝ₑ)) < 1.3
end
circular = Subset(filter(m -> incirc(m), genmodes))

# small_crystal = Crystal(unitcell, [6, 6])
# restricted_fockspace = FockSpace(Subset([setattr(m, :offset => p) for p in latticepoints(small_crystal) for m in modes]))
circle_fockspace = FockSpace(circular)

𝐹::FockMap = @time fourier(brillouin_zone(crystal), circle_fockspace) / sqrt(vol(crystal))
𝐶ᵣ::FockMap = 𝐹' * 𝐶 * 𝐹
plot(heatmap(z=real(rep(𝐶ᵣ))))
𝑈ᵣ::FockMap = eigvecsh(𝐶ᵣ)
plot(heatmap(z=real(rep(𝑈ᵣ))))
crvs = eigvalsh(𝐶ᵣ)
plot(scatter(y=map(p -> p.second, crvs), mode="markers"))
emode1::Mode = orderedmodes(𝑈ᵣ.inspace)[12]
emode2::Mode = orderedmodes(𝑈ᵣ.inspace)[13]
moderep1 = columns(𝑈ᵣ, FockSpace(Subset([emode1])))
moderep2::FockMap = columns(𝑈ᵣ, FockSpace(Subset([emode2])))
moderep = moderep1 + (FockMap(moderep2.outspace, moderep1.inspace, rep(moderep2)) * 1im)
values = columnspec(moderep)

visualize_spectrum("Mode", values)

rng = -1:0.1:1
tspec::Matrix{Float64} = zeros(Float64, length(rng), length(rng))
bspec::Matrix{Float64} = zeros(Float64, length(rng), length(rng))
for (j, k) in enumerate(rng)
    for (i, h) in enumerate(rng)
        fkk::FockMap = fourier(Subset([Point([h, k], k_space)]), Subset(orderedmodes(bm.outspace)))
        bloch::FockMap = fkk * bm * fkk'
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
