include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")
include("../src/transformations.jl")
include("../src/zer.jl")

using LinearAlgebra, PlotlyJS, OrderedCollections, SparseArrays, ColorTypes, SmithNormalForm
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting, ..Zer, ..Transformations

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')

k_space = convert(MomentumSpace, triangular)

unitcell = union(Point([1/3, 2/3], triangular), Point([2/3, 1/3], triangular))
unitcell[1]
crystal = Crystal(unitcell, [32, 32])

modes::Subset{Mode} = quantize(:pos, unitcell, 1)

m0, m1 = members(modes)
tₙ = ComplexF64(-1.)

reciprocalfock::FockSpace = crystalfock(modes, crystal)

bonds::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => tₙ])

𝐻::FockMap = hamiltonian(crystal, bonds)
visualize(𝐻, title="Hamiltonian", rowrange=1:64, colrange=1:64)

𝐶::FockMap = groundstatecorrelations(𝐻)
visualize(𝐶, title="Correlation", rowrange=1:64, colrange=1:64)

crystalpoints::Subset{Point} = latticepoints(crystal)
physicalmodes::Subset{Mode} = spanoffset(modes, crystalpoints)
restrictedregion::Subset{Mode} = filter(circularfilter(origin(euclidean(RealSpace, 2)), 2.0), physicalmodes)
restrictedfock::FockSpace = FockSpace(restrictedregion)

distill = frozenisometries(crystal, 𝐶, restrictedfock)

fillediso::FockMap = distill[:filled]
𝐹ₖ = fourier(reciprocalfock, fillediso.outspace) / sqrt(vol(crystal))
test = 𝐹ₖ * fillediso
plot(heatmap(z=real(rep(fillediso))))
# FockSpace(setattr(orderedmodes(fillediso.outspace), :offset => getattr(first(partition)
∑𝐹ₖ::Vector{FockMap} = [fourier(FockSpace(partition), fillediso.outspace) / sqrt(vol(crystal)) for partition in rep(reciprocalfock)]
filledglobalstates = [𝐹ₖ * fillediso for 𝐹ₖ in ∑𝐹ₖ]
test1 = focksum(filledglobalstates)
visualize(test1, title="Filled Global State")
filledglobalstate::FockMap = focksum([FockMap(state.outspace, FockSpace(setattr(orderedmodes(state.inspace), :offset => k)), rep(state)) for (k, state) in zip(brillouinzone(crystal), filledglobalstates)])
visualize(filledglobalstate, title="Filled Global State", rowrange=1:32, colrange=1:32)
filledprojector::FockMap = filledglobalstate * filledglobalstate'
visualize(filledprojector, title="Filled Projector", rowrange=1:32, colrange=1:32)

emptyiso::FockMap = distill[:empty]
test2 = 𝐹ₖ * emptyiso
emptyglobalstates = [𝐹ₖ * emptyiso for 𝐹ₖ in ∑𝐹ₖ]
emptyglobalstate::FockMap = focksum([FockMap(state.outspace, FockSpace(setattr(orderedmodes(state.inspace), :offset => k)), rep(state)) for (k, state) in zip(brillouinzone(crystal), emptyglobalstates)])
emptyprojector::FockMap = emptyglobalstate * emptyglobalstate'
plot(heatmap(z=real(rep(emptyiso))))

projector::FockMap = emptyprojector - filledprojector

globalvals = eigvalsh(projector)
plot(scatter(y=map(p -> p.second, globalvals), mode="markers"))
plot(heatmap(z=real(rep(projector))))

# small_crystal = Crystal(unitcell, [6, 6])
# restricted_fockspace = FockSpace(Subset([setattr(m, :offset => p) for p in latticepoints(small_crystal) for m in modes]))
circle_fockspace = FockSpace(region)

[m.attrs for m in orderedmodes(restrictedfock)]

𝐹::FockMap = fourier(brillouinzone(crystal), restrictedfock) / sqrt(vol(crystal))
𝐶ᵣ::FockMap = 𝐹' * 𝐶 * 𝐹
visualize(𝐶ᵣ, title="Restricted Correlations")
𝑈ᵣ::FockMap = eigvecsh(𝐶ᵣ)
visualize(𝑈ᵣ, title="Restricted Unitary")
crvs = eigvalsh(𝐶ᵣ)
plot(scatter(y=map(p -> p.second, crvs), mode="markers"))
emode1::Mode = orderedmodes(𝑈ᵣ.inspace)[12]
emode2::Mode = orderedmodes(𝑈ᵣ.inspace)[13]
moderep1 = columns(𝑈ᵣ, FockSpace(Subset(emode1)))
moderep2::FockMap = columns(𝑈ᵣ, FockSpace(Subset(emode2)))
moderep = moderep1 + (FockMap(moderep2.outspace, moderep1.inspace, rep(moderep2)) * 1im)
values = columnspec(moderep1)

visualize(values)

rng = -1:0.1:1
tspec::Matrix{Float64} = zeros(Float64, length(rng), length(rng))
bspec::Matrix{Float64} = zeros(Float64, length(rng), length(rng))
for (j, k) in enumerate(rng)
    for (i, h) in enumerate(rng)
        fkk::FockMap = fourier(Subset(Point([h, k], k_space)), Subset(orderedmodes(bm.outspace)))
        bloch::FockMap = fkk * bm * fkk'
        spec = eigvalsh(bloch)
        tspec[j, i] = real(spec[1].second)
        bspec[j, i] = real(spec[2].second)
    end
end

# spec = [hcat([eigvalsh(bloch(bonds, Point([h, k], k_space), crystal), :offset => Point([h, k], k_space)) for h in -1:0.1:1]...) for k in -1:0.1:1]

# top = map(p -> lineartransform(euclidean(MomentumSpace, 2), getattr(p.first, :offset)) => p.second, hcat([v[1, :] for v in spec]...))
# bottom = map(p -> lineartransform(euclidean(MomentumSpace, 2), getattr(p.first, :offset)) => p.second, hcat([v[2, :] for v in spec]...))

# topvals = hcat([map(p -> [pos(p.first)..., p.second], top)...]...)
# botvals = hcat([map(p -> [pos(p.first)..., p.second], bottom)...]...)

plot([surface(z=tspec), surface(z=bspec)])



# visualize_region("Honeycomb lattice", real_zone, euclidean(RealSpace, 2))
# visualize_region("Honeycomb lattice", reciprocal_zone, euclidean(RealSpace, 2))
