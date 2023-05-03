include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")
include("../src/zer.jl")
include("../src/transformations.jl")

using LinearAlgebra, PlotlyJS, OrderedCollections, SparseArrays, ColorTypes, SmithNormalForm
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting, ..Zer, ..Transformations

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
scale = Scale([2 0; 0 2])

k_space = convert(MomentumSpace, triangular)

unitcell = union(Point([1/3, 2/3], triangular), Point([2/3, 1/3], triangular))
pp = scale * unitcell[1]
invscale = Scale(inv(rep(scale)))
ip = invscale * pp
crystal = Crystal(unitcell, [32, 32])

scaledcrystal = scale * (scale * crystal)
visualize_region("space", scaledcrystal.unitcell, euclidean(RealSpace, 2))

modes::Subset{Mode} = quantize("physical", :pos, unitcell, 1)
m0, m1 = members(modes)
tₙ = ComplexF64(-1.)

reciprocalfock::FockSpace = sparsefock(modes, brillouinzone(crystal))

bonds::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => tₙ])
plot(heatmap(z=real(rep(bonds))))


𝐻::FockMap = hamiltonian(crystal, bonds)
plot(heatmap(z=real(rep(𝐻)[1:64,1:64])))
filled::FockMap = groundstates(𝐻)

𝐶::FockMap = trivial(filled.outspace) - filled * filled'
plot(heatmap(z=real(rep(𝐶)[1:128,1:128])))

crystalpoints::Subset{Point} = latticepoints(crystal)
physicalmodes::Subset{Mode} = spanbasis(modes, crystalpoints)
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
plot(heatmap(z=real(rep(test1)[1:8, :])))
filledglobalstate::FockMap = focksum([FockMap(state.outspace, FockSpace(setattr(orderedmodes(state.inspace), :offset => k)), rep(state)) for (k, state) in zip(brillouinzone(crystal), filledglobalstates)])
filledprojector::FockMap = filledglobalstate * filledglobalstate'
plot(heatmap(z=real(rep(filledprojector))))

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

𝐹::FockMap = @time fourier(brillouinzone(crystal), restrictedfock) / sqrt(vol(crystal))
𝐶ᵣ::FockMap = 𝐹' * 𝐶 * 𝐹
plot(heatmap(z=real(rep(𝐶ᵣ))))
𝑈ᵣ::FockMap = eigvecsh(𝐶ᵣ)
plot(heatmap(z=real(rep(𝑈ᵣ))))
crvs = eigvalsh(𝐶ᵣ)
plot(scatter(y=map(p -> p.second, crvs), mode="markers"))
emode1::Mode = orderedmodes(𝑈ᵣ.inspace)[23]
emode2::Mode = orderedmodes(𝑈ᵣ.inspace)[24]
moderep1 = columns(𝑈ᵣ, FockSpace(Subset([emode1])))
moderep2::FockMap = columns(𝑈ᵣ, FockSpace(Subset([emode2])))
moderep = moderep1 + (FockMap(moderep2.outspace, moderep1.inspace, rep(moderep2)) * 1im)
values = columnspec(moderep1)

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

# top = map(p -> lineartransform(euclidean(MomentumSpace, 2), getattr(p.first, :offset)) => p.second, hcat([v[1, :] for v in spec]...))
# bottom = map(p -> lineartransform(euclidean(MomentumSpace, 2), getattr(p.first, :offset)) => p.second, hcat([v[2, :] for v in spec]...))

# topvals = hcat([map(p -> [pos(p.first)..., p.second], top)...]...)
# botvals = hcat([map(p -> [pos(p.first)..., p.second], bottom)...]...)

plot([surface(z=tspec), surface(z=bspec)])



# visualize_region("Honeycomb lattice", real_zone, euclidean(RealSpace, 2))
# visualize_region("Honeycomb lattice", reciprocal_zone, euclidean(RealSpace, 2))
