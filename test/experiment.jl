include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")
include("../src/transformations.jl")
include("../src/zer.jl")

using PlotlyJS, LinearAlgebra, OrderedCollections, SparseArrays, ColorTypes, SmithNormalForm
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Transformations, ..Zer
using ..Plotting

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

unitcell = Subset([triangular & [1/3, 2/3], triangular & [2/3, 1/3]])
crystal = Crystal(unitcell, [32, 32])
modes::Subset{Mode} = quantize(:pos, unitcell, 1)

tₙ = ComplexF64(-1.)
m0, m1 = members(modes)

c3 = Symmetry(
    group=:c3,
    pointgrouprep=[cos(2π/3) -sin(2π/3); sin(2π/3) cos(2π/3)],
    order=3,
    shift=(euclidean(RealSpace, 2) |> origin),
    irreps=Dict(
        :s => Irrep(exp(0im)),
        :pplus => Irrep(exp(1im * 2π/3)),
        :pminus => Irrep(exp(2im * 2π/3))
    ))

c6 = Symmetry(
    group=:c6,
    pointgrouprep=[cos(2π/6) -sin(2π/6); sin(2π/6) cos(2π/6)],
    order=6,
    shift=(euclidean(RealSpace, 2) |> origin),
    irreps=Dict(:s => Irrep(exp(0im))))

bonds::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => tₙ])

𝐻::FockMap = hamiltonian(crystal, bonds)
𝐶::FockMap = groundstatecorrelations(𝐻)

crystalfock = 𝐻.inspace

blocked = blocking(:scale => Scale([2. 0.; 0. 2.]), :correlations => 𝐶, :crystal => crystal)
blockedcorrelations::FockMap = blocked[:correlations]

visualize(𝐻, title="Hamiltonian", rowrange=1:64, colrange=1:64)
visualize(𝐶, title="Correlation", rowrange=1:64, colrange=1:64)
visualize(blockedcorrelations, title="Correlation", rowrange=1:64, colrange=1:64)

newcrystal = blocked[:crystal]

crystalpoints::Subset{Point} = latticepoints(newcrystal)
newmodes::Subset{Mode} = quantize(:pos, newcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(newmodes, crystalpoints)
restrictedregion::Subset{Mode} = filter(circularfilter(origin(euclidean(RealSpace, 2)), 2.0), physicalmodes)
restrictedfock::FockSpace = FockSpace(restrictedregion)

isometries = localfrozenisometries(blocked[:correlations], restrictedfock, selectionstrategy=frozenselectionbycount(3))

function fourierisometries(; localisometry::FockMap, crystalfock::FockSpace{Crystal})::Dict{Point, FockMap}
    crystal::Crystal = crystalof(crystalfock)
    fouriermap::FockMap = fourier(crystalfock, localisometry.outspace) / (crystal |> vol |> sqrt)
    momentumfouriers::Vector{FockMap} = rowsubmaps(fouriermap)
    bz::Subset{Point} = brillouinzone(crystal)
    return Dict(k => fourier_k * localisometry for (k, fourier_k) in zip(bz, momentumfouriers))
end

function isometryglobalprojector(; localisometry::FockMap, crystalfock::FockSpace{Crystal})
    momentuisometries::Dict{Point, FockMap} = fourierisometries(localisometry=localisometry, crystalfock=crystalfock)
    crystal::Crystal = crystalof(crystalfock)
    bz::Subset{Point} = brillouinzone(crystal)
    globalprojector::FockMap = focksum(map(k -> momentuisometries[k] * momentuisometries[k]', bz))
    return FockMap(
        globalprojector,
        outspace=FockSpace(globalprojector.outspace, reflected=crystal),
        inspace=FockSpace(globalprojector.inspace, reflected=crystal))
end

filledglobalprojector = isometryglobalprojector(localisometry=isometries[:filled], crystalfock=blockedcorrelations.inspace)
visualize(filledglobalprojector, rowrange=1:128, colrange=1:128)
emptyglobalprojector = isometryglobalprojector(localisometry=isometries[:empty], crystalfock=blockedcorrelations.inspace)
visualize(emptyglobalprojector, rowrange=1:128, colrange=1:128)

globalprojector = emptyglobalprojector - filledglobalprojector
visualize(globalprojector, rowrange=1:64, colrange=1:64)

visualize(groupelement(c6, 2) * globalprojector.outspace, rowrange=1:64, colrange=1:64)

symmetrytransformers::Vector{FockMap} = [element * globalprojector.outspace for element in groupelements(c6)]

symmetrizedglobalprojector = focksum([transformer * globalprojector * transformer' for transformer in symmetrytransformers])

visualize(symmetrizedglobalprojector, rowrange=1:64, colrange=1:64)

globaldistillspec = eigvalsh(symmetrizedglobalprojector)
filter(p -> -1e-3 < p.second < 1e-3, globaldistillspec)
plot(scatter(y=map(p -> p.second, globaldistillspec), mode="markers"))

fouriersubspaces = crystalsubspaces(symmetrizedglobalprojector.inspace)

kprojectors = [(k, restrict(symmetrizedglobalprojector, subspace, subspace)) for (k, subspace) in fouriersubspaces]

spectrum = Dict()

for (k, proj) in kprojectors
    kspec = map(p -> p.second, eigvalsh(proj)) |> sort
    spectrum[k] = kspec
end

spectrum

using Compat

mesh = brillouinmesh(newcrystal)

spec3d = stack(map(k -> spectrum[k], mesh))

plot([surface(z=spec3d[n, :, :]) for n in 1:size(spec3d, 1)])

specmat = hcat([v for (k, v) in spectrum]...)
plot([scatter(y=specmat[n, :]) for n in 1:size(specmat, 1)])

function Base.:*(symmetry::Symmetry, crystalfock::FockSpace{Crystal})::Vector{FockMap}
    homefock::FockSpace = unitcellfock(crystalfock)
    homemaps::Vector{FockMap} = symmetry * homefock
    momentumsubspaces::Dict{Point, FockSpace} = crystalsubspaces(crystalfock)
    symmetrizedbzs::Subset{Subset{Point}} = symmetry * (crystalfock |> crystalof |> brillouinzone)

    function computesymmetrizetransformer(homemap::FockMap, symmetrizedbz::Subset{Point})::FockMap
        fouriermap::FockMap = fourier(crystalfock, homefock)
        symmetrizedfouriermap::FockMap = fourier(crystalfock, homemap.outspace)
        symmetrizedfockpairs::Vector{Pair{FockSpace, FockSpace}} = [subspace => momentumsubspaces[symmetrizedbz[n]] for (n, subspace) in enumerate(crystalfock |> subspaces)]
        fockmap::FockMap = focksum(
            [rows(symmetrizedfouriermap, subspace) * homemap * rows(fouriermap, subspace)' for subspace in crystalfock |> subspaces])
        return FockMap(
            FockSpace(fockmap.outspace, reflected=crystalof(crystalfock)),
            FockSpace(fockmap.inspace, reflected=crystalof(crystalfock)),
            rep(fockmap))
    end

    return map(computesymmetrizetransformer, homemaps)
end

Base.:*(symmetry::Symmetry, region::Subset{Point})::Subset{Point} = Subset(
    Point((p |> spaceof |> rep |> inv) * sym * ((p - symmetry.center) |> euclidean |> pos), spaceof(p)) + symmetry.center for p in region
    for sym in Transformations.groupreps(symmetry))

c3 * (crystal |> brillouinzone)

symmetrymaps = c3 * globalprojector.outspace
visualize(symmetrymaps[2], rowrange=1:64, colrange=1:64)

ss = symmetrymaps[2]

sg = ss * globalprojector * ss'

symmetrizedglobalprojector = focksum([s * globalprojector * s' for s in symmetrymaps])
visualize(sg, rowrange=1:64, colrange=1:64)

globaldistillspec = eigvalsh(globalprojector)
filter(p -> -1e-3 < p.second < 1e-3, globaldistillspec)
plot(scatter(y=map(p -> p.second, globaldistillspec), mode="markers"))

emptyprojector::FockMap = isometries[:empty] * isometries[:empty]'
filledprojector::FockMap = isometries[:filled] * isometries[:filled]'
projector::FockMap = emptyprojector - filledprojector
pvs, pu = eigh(projector)
plot(scatter(y=map(p -> p.second, pvs), mode="markers"))
visualize(pu)
cmode::Mode = orderedmodes(pu.inspace)[4]
cmr::FockMap = columns(pu, FockSpace(Subset(cmode)))
cspec = columnspec(cmr)
visualize(cspec)

restrictfourier::FockMap = fourier(blockedcorrelations.outspace, restrictedfock) / sqrt(vol(newcrystal))
visualize(restrictfourier, rowrange=1:32)
restrictedcorrelations::FockMap = restrictfourier' * blockedcorrelations * restrictfourier
visualize(restrictedcorrelations)
crvs = eigvalsh(restrictedcorrelations)
plot(scatter(y=map(p -> p.second, crvs), mode="markers"))

restrictedunitary::FockMap = eigvecsh(restrictedcorrelations)
visualize(restrictedunitary)

emode1::Mode = orderedmodes(restrictedunitary.inspace)[13]
emode2::Mode = orderedmodes(restrictedunitary.inspace)[14]
mr1::FockMap = columns(restrictedunitary, FockSpace(Subset(emode1)))
mr2::FockMap = columns(restrictedunitary, FockSpace(Subset(emode2)))
mr = FockMap(mr1.outspace, mr1.inspace, rep(mr1) + 1im * rep(mr2))
values = columnspec(mr)

visualize(values)
