using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Revise
using Zipper

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [48, 48])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = Subset(m for m in quantize(unitcell, 1))
m0, m1 = members(modes)

tₙ = -1 + 0im
tₕ = 0.1im

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:r => [-1, 0] ∈ triangular)) => tₙ,
    (m0, m1 |> setattr(:r => [0, 1] ∈ triangular)) => tₙ]

haldane = [
    (m0, setattr(m0, :r => Point([1, 1], triangular))) => tₕ,
    (m0, setattr(m0, :r => Point([-1, 0], triangular))) => tₕ,
    (m0, setattr(m0, :r => Point([0, -1], triangular))) => tₕ,
    (m1, setattr(m1, :r => Point([1, 1], triangular))) => -tₕ,
    (m1, setattr(m1, :r => Point([-1, 0], triangular))) => -tₕ,
    (m1, setattr(m1, :r => Point([0, -1], triangular))) => -tₕ]

bonds::FockMap = bondmap([nearestneighbor..., haldane...])

bondmodes::Subset{Mode} = bonds|>getoutspace|>orderedmodes
bondmodes|>collect
basismodes::Subset{Mode} = bondmodes|>removeattr(:r)
fockspace::CrystalFock = getcrystalfock(basismodes, crystal)
transform::FockMap = fourier(fockspace, bondmodes|>RegionFock)

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize

groundstateprojector = groundstates |> crystalprojector
C = idmap(groundstateprojector.outspace) - groundstateprojector

# inplaceadd(a::FockMap, b::FockMap)::FockMap = a + FockMap(b, outspace=a|>getoutspace, inspace=a|>getinspace, performpermute=false)

# function computequantumgeometrictensor2d(state::CrystalSpectrum)
#     statevectors::Dict{Momentum, FockMap} = state|>geteigenvectors
#     Δkx = (1 / size(state|>getcrystal)[1], 0) ∈ (convert(MomentumSpace, state|>getcrystal|>getspace))
#     Δky = (0, 1 / size(state|>getcrystal)[2]) ∈ (convert(MomentumSpace, state|>getcrystal|>getspace))

#     kprojectors::Dict{Momentum, FockMap} = Dict(k => kstate * kstate' for (k, kstate) in statevectors)
#     kcorrelations::Dict{Momentum, FockMap} = Dict(k => idmap(projector|>getoutspace) - projector for (k, projector) in kprojectors)
#     ΔCkxs::Dict{Momentum, FockMap} = Dict(k => inplaceadd(-kcorr, kcorrelations[k + Δkx |> basispoint]) for (k, kcorr) in kcorrelations)
#     ΔCkys::Dict{Momentum, FockMap} = Dict(k => inplaceadd(-kcorr, kcorrelations[k + Δky |> basispoint]) for (k, kcorr) in kcorrelations)
    
#     function kqgtensor(k::Momentum)::FockMap
#         gxx = statevectors[k]' * (ΔCkxs[k]' * kcorrelations[k] * ΔCkxs[k]) * statevectors[k]
#         gxy = statevectors[k]' * (ΔCkxs[k]' * kcorrelations[k] * ΔCkys[k]) * statevectors[k]
#         gyx = statevectors[k]' * (ΔCkys[k]' * kcorrelations[k] * ΔCkxs[k]) * statevectors[k]
#         gyy = statevectors[k]' * (ΔCkys[k]' * kcorrelations[k] * ΔCkys[k]) * statevectors[k]

#         xmode::Mode = Mode(:offset=>k, :axis=>(:x))
#         ymode::Mode = Mode(:offset=>k, :axis=>(:y))

#         fockspace::FockSpace = FockSpace([xmode, ymode])
#         return FockMap(fockspace, fockspace, [rep(gxx)[1,1] rep(gxy)[1,1]; rep(gyx)[1,1] rep(gyy)[1,1]])
#     end

#     return Dict(k=>k|>kqgtensor for (k, _) in statevectors)
# end

# Base.:real(fockmap::FockMap)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, fockmap|>rep|>real)
# Base.:imag(fockmap::FockMap)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, fockmap|>rep|>imag)

# function computequantummetric(state::CrystalSpectrum)
#     gtensors = computequantumgeometrictensor2d(state)
#     barecrystal = Crystal(state|>getcrystal|>getspace|>getorigin, state|>getcrystal|>size)

#     gxx = directsum(tensor[1,1] for (k, tensor) in gtensors)
#     gxxbarecrystalfock = FockSpace(gxx|>getoutspace, reflected=barecrystal)
#     gxx = FockMap(gxx, inspace=gxxbarecrystalfock, outspace=gxxbarecrystalfock, performpermute=false)|>real

#     gyy = directsum(tensor[2,2] for (k, tensor) in gtensors)
#     gyybarecrystalfock = FockSpace(gyy|>getoutspace, reflected=barecrystal)
#     gyy = FockMap(gyy, inspace=gyybarecrystalfock, outspace=gyybarecrystalfock, performpermute=false)|>real

#     gxy = directsum(tensor[1,2] for (k, tensor) in gtensors)
#     gxybarecrystalfock = FockSpace(gxy|>getoutspace, reflected=barecrystal)
#     gxy = FockMap(gxy, inspace=gxybarecrystalfock, outspace=gxybarecrystalfock, performpermute=false)|>real

#     gyx = directsum(tensor[2,1] for (k, tensor) in gtensors)
#     gyxbarecrystalfock = FockSpace(gyx|>getoutspace, reflected=barecrystal)
#     gyx = FockMap(gyx, inspace=gyxbarecrystalfock, outspace=gyxbarecrystalfock, performpermute=false)|>real

#     return Dict(:gxx => gxx, :gyy => gyy, :gxy => gxy, :gyx => gyx)
# end



# gsqm = computequantummetric(groundstates)
# visualize(gsqm[:gxx]|>crystalspectrum, usecontour=true, title="gxx-GS")
# visualize(gsqm[:gyy]|>crystalspectrum, usecontour=true, title="gyy-GS")

# gsqmdet = quantummetricdeterminants(groundstates)
# visualize(gsqmdet|>crystalspectrum, usecontour=true, title="qgtdet-GS")

# kqgtensors = computequantumgeometrictensor2d(groundstates)
# bc = directsum(tensor[1,2]|>imag for (k, tensor) in kqgtensors)
# barecrystal = Crystal(groundstates|>getcrystal|>getspace|>getorigin, groundstates|>getcrystal|>size)
# barecrystalfock = FockSpace(bc|>getoutspace, reflected=barecrystal)
# berrycurvature = 2 * FockMap(bc, inspace=barecrystalfock, outspace=barecrystalfock, performpermute=false)

# (berrycurvature|>rep|>sum) / (2π)

# visualize(berrycurvature|>crystalspectrum, usecontour=true, title="berrycurvature-GS")

# groundstateeigenvectors = groundstates |> geteigenvectors
# δkx = (1 / size(crystal)[1], 0) ∈ kspace
# δky = (0, 1 / size(crystal)[2]) ∈ kspace

# δIxs = Dict(k => FockMap(groundstateeigenvectors[k + δkx |> basispoint], inspace=I|>getinspace, outspace=I|>getoutspace, performpermute=false) - I for (k, I) in groundstateeigenvectors)
# δIys = Dict(k => FockMap(groundstateeigenvectors[k + δky |> basispoint], inspace=I|>getinspace, outspace=I|>getoutspace, performpermute=false) - I for (k, I) in groundstateeigenvectors)
# δIx = directsum(I for (_, I) in δIxs)
# δIy = directsum(I for (_, I) in δIys)

# Base.:real(fockmap::FockMap)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, fockmap|>rep|>real)

# qgtcrystal = Crystal(triangular |> getorigin, crystal|>size)
# qgt = δIx' * groundstateprojector * δIy
# qgtoutspace = FockSpace(qmt|>getoutspace, reflected=qmtcrystal)
# quantumgeometrictensor = idmap(qgtoutspace) - FockMap(qmt, inspace=qmtoutspace, outspace=qmtoutspace)
# quantummetrictensor = quantumgeometrictensor |> real
# quantummetrictensor |> crystalspectrum |> visualize

# C |> crystalspectrum |> visualize

entanglemententropy(C |> crystalspectrum)

function zer(correlations::FockMap)
    crystalfock = correlations.outspace

    scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)
    @time blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

    blockedcrystal::Crystal = blockresult[:crystal]
    blockedcorrelations::FockMap = blockresult[:correlations]
    blocker = blockresult[:transformer]

    function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
        currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthospace
        physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
        return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
    end

    crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
    samplepoints::Subset{Offset} = crystalpoints + c6^2 * crystalpoints + c6^4 * crystalpoints
    blockedmodes::Subset{Mode} = quantize(blockedcrystal.unitcell, 1)|>orderedmodes
    physicalmodes::Subset{Mode} = spanoffset(blockedmodes, samplepoints)


    frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 2)
    # visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
    frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)

    @time globaldistiller = globaldistillerhamiltonian(
        correlations=blockresult[:correlations],
        restrictspace=frozenseedingfock,
        localisometryselectionstrategy=frozenselectionbycount(3))

    globaldistillerspectrum = globaldistiller |> crystalspectrum
    # visualize(globaldistillerspectrum, title="Global Distiller")

    @time distillresult = distillation(globaldistillerspectrum, :courier => v -> abs(v) < 1e-5, :empty => v -> v > 1e-5, :filled => v -> v < -1e-5)

    courierseedingcenter::Offset = [2/3, 1/3] ∈ (blockedmodes |> getspace)
    courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
    courierseedingregion::Subset{Offset} = Subset(m |> getpos for m in courierseedingmodes)
    # visualize(courierseedingregion, courierseedingcenter |> Subset, title="Courier Seeding Region", visualspace=euclidean(RealSpace, 2))
    courierseedingfock::FockSpace{Region} = FockSpace{Region}(courierseedingmodes)

    c3 = c6^2 |> recenter(courierseedingcenter)

    blockedcourierprojector = distillresult[:courier] |> crystalprojector
    blockedcouriercorrelation = idmap(blockedcourierprojector.outspace, blockedcourierprojector.outspace) - blockedcourierprojector

    # return globaldistillerspectrum

    findlocalspstates(statecorrelations=blockedcouriercorrelation, regionfock=courierseedingfock, symmetry=c3, spectrumextractpredicate=v -> v < 5e-2)|>println
    @time localcourierseed = findlocalspstates(statecorrelations=blockedcouriercorrelation, regionfock=courierseedingfock, symmetry=c3, spectrumextractpredicate=v -> v < 5e-2)[1]
    fullcourierseed = localcourierseed + (c6 * localcourierseed.outspace) * localcourierseed * (c6 * localcourierseed.inspace)'

    crystalcourierseeds = crystalisometries(localisometry=fullcourierseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)

    wanniercourierisometry = wannierprojection(
        crystalisometries=distillresult[:courier].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds)

    # regionalrestriction(wanniercourierisometry, courierseedingfock) |> visualize

    couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
    # couriercorrelations |> crystalspectrum |> visualize
    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    # purifiedcorrelationspectrum |> visualize
    unpurecorrelations = couriercorrelations
    couriercorrelations = purifiedcorrelationspectrum |> FockMap

    blockedfilledprojector = distillresult[:filled] |> crystalprojector
    blockedfilledcorrelation = idmap(blockedfilledprojector.outspace, blockedfilledprojector.outspace) - blockedfilledprojector
    @time filledseed = findlocalspstates(statecorrelations=blockedfilledcorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalfilledseeds = crystalisometries(localisometry=filledseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannierfilledisometry = wannierprojection(
        crystalisometries=distillresult[:filled].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledseeds)

    blockedemptyprojector = distillresult[:empty] |> crystalprojector
    blockedemptycorrelation = idmap(blockedemptyprojector.outspace, blockedemptyprojector.outspace) - blockedemptyprojector
    @time emptyseed = findlocalspstates(statecorrelations=blockedemptycorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalemptyseeds = crystalisometries(localisometry=emptyseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannieremptyisometry = wannierprojection(
        crystalisometries=distillresult[:empty].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyseeds)

    filledcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
    emptycorrelations = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry

    Ft = fourier(blockedcorrelations.outspace, frozenseedingfock)

    return Dict(
        :blocker => blocker,
        :correlations => couriercorrelations,
        :frozenregioncorrelations => Ft' * blockedcorrelations * Ft,
        :unpurecorrelations => unpurecorrelations,
        :courierzipper => wanniercourierisometry' * blocker, 
        :filledzipper => wannierfilledisometry' * blocker,
        :emptyzipper => wannieremptyisometry' * blocker,
        :emptyisometry => wannieremptyisometry,
        :filledisometry => wannierfilledisometry,
        :courierisometry => wanniercourierisometry,
        :globaldistiller => globaldistiller,
        :filledcorrelations => filledcorrelations, 
        :emptycorrelations => emptycorrelations,
        :frozenseedingfock => frozenseedingfock,
        :entanglemententropy => entanglemententropy(couriercorrelationspectrum) / (couriercorrelations|>getoutspace|>getcrystal|>vol))
end

rg1 = zer(C)
rg1[:blocker]|>visualize
scale = Scale([2 0; 0 2], crystal|>getspace)
blocker = scale * getoutspace(rg1[:correlations])
blocker|>visualize
rg1bc = blocker * rg1[:correlations] * blocker'
rg1bc|>crystalspectrum|>visualize

    crystalfock = rg1[:correlations]|>getoutspace
    crystal = crystalfock |> getcrystal
    scaledcrystal::Crystal = scale * crystal
    unscaledblockedregion::Subset{Offset} = (scale |> inv) * scaledcrystal.unitcell
    bz::Subset{Momentum} = crystal |> brillouinzone
    basismodes::Subset{Mode} = crystalfock|>unitcellfock|>orderedmodes
    scaledbz::Subset{Momentum} = scaledcrystal |> brillouinzone

    momentummappings::Base.Generator = (basispoint(scale * p) => p for p in bz)
    mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point, Vector{Point}}()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end

    ksubsets::Dict{Momentum, Subset{Mode}} = crystalfock |> crystalsubsets
    scaledfock::FockSpace = ((ksubsets[k] for k in mappingpartitions[kscaled]) |> subsetunion |> FockSpace
        for kscaled in scaledbz) |> fockspaceunion

    # TODO: Use of latticeoff requires revisit, since latticeoff only truncates the decimals.
    unscaledblockedlatticeoffsets::Region = Subset(pbc(crystal, p) |> latticeoff for p in unscaledblockedregion)
    unscaledblockedunitcellfock::RegionFock = spanoffset(basismodes, unscaledblockedlatticeoffsets)|>RegionFock

    restrictedfourier::FockMap = fourier(crystalfock, unscaledblockedunitcellfock)
    crystalfock|>unitcellfock|>modeattrs
    unscaledblockedunitcellfock|>unitcellfock|>modeattrs
    restrictedfourier|>visualize

    volumeratio::Number = (crystal |> vol) / (scaledcrystal |> vol)
    permutedfourier::FockMap = restrictedfourier / sqrt(volumeratio)

    function repackfourierblocks(source::FockMap, kscaled::Momentum, partition::Subset{Mode})::FockMap
        partitionrows::FockMap = rows(source, partition |> FockSpace)
        inspace::FockSpace = FockSpace(mode|>setattr(:k=>kscaled, :b=>*(scale, mode|>getpos))|>removeattr(:r) for mode in partitionrows|>getinspace)
        return FockMap(partitionrows.outspace, inspace, partitionrows |> rep)
    end

    repackedblocks::Base.Generator = (
        repackfourierblocks(permutedfourier, kscaled, partition)
        for (kscaled, partition) in Iterators.zip(scaledbz, scaledfock |> rep))
    repackedblocks|>first|>getinspace|>modeattrs
    blocker = directsum(repackedblocks)
    blocker = FockMap(blocker, inspace=FockSpace(blocker |> getinspace, reflected=scaledcrystal), outspace=crystalfock)'
    blocker * rg1[:correlations] * blocker'|>crystalspectrum|>visualize

scale * (rg1[:correlations]|>getoutspace|>getcrystal)
scale * (rg1[:correlations]|>getoutspace)|>visualize

rg2 = zer(rg1[:correlations])

rg2|>visualize

rg3 = zer(rg2[:correlations])

rg4 = zer(rg3[:correlations])

rg5 = zer(rg4[:correlations])

rg1[:entanglemententropy] / (rg1[:correlations]|>getoutspace|>getcrystal|>vol)
rg2[:entanglemententropy] / (rg2[:correlations]|>getoutspace|>getcrystal|>vol)
rg3[:entanglemententropy] / (rg3[:correlations]|>getoutspace|>getcrystal|>vol)
rg4[:entanglemententropy] / (rg4[:correlations]|>getoutspace|>getcrystal|>vol)

rg1filledC = rg1[:filledzipper]' * rg1[:filledzipper]
rg2filledC = rg1[:courierzipper]' * rg2[:filledzipper]' * rg2[:filledzipper] * rg1[:courierzipper]
rg3filledC = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:filledzipper]' * rg3[:filledzipper] * rg2[:courierzipper] * rg1[:courierzipper]
rg4filledC = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:filledzipper]' * rg4[:filledzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]
rg4P = groundstatespectrum(rg4[:correlations]|>crystalspectrum, perunitcellfillings=1)|>crystalprojector
rg4C = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:courierzipper]' * rg4P * rg4[:courierzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]

rgP = rg1filledC + rg2filledC + rg3filledC + rg4filledC + rg4C
rgC = idmap(rgP|>getoutspace) - rgP
visualize(rgC, colrange=1:512, rowrange=1:512)
rgC|>crystalspectrum|>visualize

visualize(rgC - C, colrange=1:512, rowrange=1:512)
M = rgC - C
(M * M' |> tr) / 

function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = blockedcrystal |> getspace |> orthospace
    physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
end

# blockedcrystal = rg1[:correlations]|>getoutspace|>getcrystal
# crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
# samplepoints::Subset{Offset} = crystalpoints + c6^2 * crystalpoints + c6^4 * crystalpoints
# blockedmodes::Subset{Mode} = quantize(:b, blockedcrystal.unitcell, 1)
# physicalmodes::Subset{Mode} = spanoffset(blockedmodes, samplepoints)


# frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 8)
# frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)


# regionalrestriction(rg1[:filledisometry], frozenseedingfock) |> visualize

# rg1[:courierzipper]' * rg2[:blocker]' * rg2[:filledisometry]

# regionalrestriction(rg1[:courierisometry] * rg2[:blocker]' * rg2[:filledisometry], frozenseedingfock) |> visualize
# rg1[:courierisometry] * rg2[:courierzipper]' * rg3[:blocker]' * rg3[:filledisometry]
# regionalrestriction(rg1[:courierisometry] * rg2[:courierzipper]' * rg3[:blocker]' * rg3[:filledisometry], frozenseedingfock) |> visualize

# rg3[:correlations] |> crystalspectrum |> visualize

# regionalrestriction(rg3[:blocker]' * rg3[:filledisometry], rg2[:frozenseedingfock]) |> visualize

# rg1[:entanglemententropy]

# rg1[:globaldistiller] |> crystalspectrum |> visualize
# rg1[:filledisometry] |> getinspace |> unitcellfock |> modeattrs

# rg1[:filledcorrelations] |> crystalspectrum |> visualize

H = energyspectrum |> FockMap

function Chern_number_multiband(state::CrystalSpectrum)
    statevectors::Dict{Momentum, FockMap} = state |> geteigenvectors
    statecrystal = state |> getcrystal
    Δkx = (1 / size(statecrystal)[1], 0) ∈ kspace
    Δky = (0, 1 / size(statecrystal)[2]) ∈ kspace
    Δk = [Δkx, Δky]
    
    function LinkVariable(k::Momentum, direction)
        delta_k = Δk[direction]
        overlap = (statevectors[k] |> rep)' * (statevectors[k+delta_k|> basispoint]|> rep)
        U_mu = det(overlap)
        U_mu /= abs(U_mu)
        return U_mu
    end
    
    function Berry_phase(k::Momentum)
        return LinkVariable(k, 1) * LinkVariable(k + Δk[1]|> basispoint, 2) / (LinkVariable(k + Δk[2]|> basispoint, 1) * LinkVariable(k, 2))
    end

    function Berry_curvature(k::Momentum)
        return -1im*log(Berry_phase(k))
    end

    Chern = 0 
    for kpoints in statevectors
        k = kpoints[1]
        Chern += Berry_curvature(k)
    end
    return Chern/(2*pi)
end

function berrycurvaturemultiband(state::CrystalSpectrum)
    statevectors::Dict{Momentum, FockMap} = state |> geteigenvectors
    statecrystal = state |> getcrystal
    Δkx = (1 / size(statecrystal)[1], 0) ∈ kspace
    Δky = (0, 1 / size(statecrystal)[2]) ∈ kspace
    Δk = [Δkx, Δky]
    
    function LinkVariable(k::Momentum, direction)
        delta_k = Δk[direction]
        overlap = (statevectors[k] |> rep)' * (statevectors[k+delta_k|> basispoint]|> rep)
        U_mu = det(overlap)
        U_mu /= abs(U_mu)
        return U_mu
    end
    
    function Berry_phase(k::Momentum)
        return LinkVariable(k, 1) * LinkVariable(k + Δk[1]|> basispoint, 2) / (LinkVariable(k + Δk[2]|> basispoint, 1) * LinkVariable(k, 2))
    end

    function Berry_curvature(k::Momentum)
        return -1im*log(Berry_phase(k))
    end

    function kberrycurvature(k)
        mode = Mode(:offset => k, :pos => statecrystal|>getspace|>getorigin)
        fockspace = FockSpace(mode)
        return FockMap(fockspace, fockspace, [Berry_curvature(k)][:, :])
    end

    fockmap = directsum(kberrycurvature(k) for (k, _) in statevectors)
    barecrystal = Crystal(statecrystal|>getspace|>getorigin, statecrystal|>size)
    crystalfock = FockSpace(fockmap|>getoutspace, reflected=barecrystal)
    return FockMap(fockmap, inspace=crystalfock, outspace=crystalfock)
end
:b
scale

Chern_number_multiband(groundstates)
berrycurvaturemap = berrycurvaturemultiband(groundstates)
visualize(berrycurvaturemap|>crystalspectrum, usecontour=true, title="berrycurvature-conventional-GS")

rg1courierprojector = rg1[:courierzipper]' * rg1[:courierzipper]
rg1H = rg1courierprojector * H * rg1courierprojector'
rg1gs = groundstatespectrum(rg1H |> crystalspectrum, perunitcellfillings=1)
visualize(rg1gs, usecontour=true, title="gsH-RG1")

rg1H = rg1[:courierzipper] * H * rg1[:courierzipper]'
rg1gs = groundstatespectrum(rg1H|>crystalspectrum, perunitcellfillings=1)
visualize(rg1gs, usecontour=true, title="gsH-RG1-RG1space")

using SparseArrays
function quantummetricdeterminants(state::CrystalSpectrum)
    spatialorigin = state|>getcrystal|>getspace|>getorigin
    barecrystal = Crystal(spatialorigin|>Subset, state|>getcrystal|>size)
    
    function makedetblock(k, kgeometrictensor)
        mode = Mode(:offset => k, :pos => spatialorigin)
        fockspace = FockSpace(mode)
        return FockMap(fockspace, fockspace, [kgeometrictensor|>rep|>imag|>det][:, :]|>SparseMatrixCSC)
    end

    kqgtensors = computequantumgeometrictensor2d(state)
    product = directsum(makedetblock(k, kgeometrictensor) for (k, kgeometrictensor) in kqgtensors)
    fockspace = FockSpace(product|>getoutspace, reflected=barecrystal)
    return FockMap(product, inspace=fockspace, outspace=fockspace)
end
:b
rg1det = quantummetricdeterminants(rg1gs)
visualize(rg1det, colrange=512:1024, rowrange=512:1024)
visualize(rg1det|>crystalspectrum, usecontour=true, title="qgtdet-RG1-RG1space")

rg1occ = momentumoccupations(rg1[:correlations])
rg1occ|>crystalspectrum|>geteigenvalues

rg1qm = computequantummetric(rg1gs)
visualize(rg1qm[:gxx], colrange=512:1024, rowrange=512:1024)

visualize(rg1qm[:gxx]|>crystalspectrum, usecontour=true, title="gxx-RG1-RG1space")
visualize(rg1qm[:gyy]|>crystalspectrum, usecontour=true, title="gyy-RG1-RG1space")
visualize(rg1qm[:gxy]|>crystalspectrum, usecontour=true, title="gxy-RG1-RG1space")
visualize(rg1qm[:gyx]|>crystalspectrum, usecontour=true, title="gyx-RG1-RG1space")

rg2courierisometry = rg2[:courierzipper] * rg1[:courierzipper]
rg2courierprojector = rg2courierisometry' * rg2courierisometry
rg2H = rg2courierprojector * H * rg2courierprojector'
rg2gs = groundstatespectrum(rg2H |> crystalspectrum, perunitcellfillings=1)
visualize(rg2gs, usecontour=true, title="gsH-RG2")

rg2H = rg2[:courierzipper] * rg1H * rg2[:courierzipper]'
rg2gs = groundstatespectrum(rg2H|>crystalspectrum, perunitcellfillings=1)
visualize(rg2gs, usecontour=true, title="gsH-RG2-RG2space")

rg2det = quantummetricdeterminants(rg2gs)
visualize(rg2det|>crystalspectrum, usecontour=true, title="qgtdet-RG2-RG2space")

rg2qm = computequantummetric(rg2gs)
visualize(rg2qm[:gxx]|>crystalspectrum, usecontour=true, title="gxx-RG2-RG2space")
visualize(rg2qm[:gyy]|>crystalspectrum, usecontour=true, title="gyy-RG2-RG2space")

rg3courierisometry = rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]
rg3courierprojector = rg3courierisometry' * rg3courierisometry
rg3H = rg3courierprojector * H * rg3courierprojector'
rg3gs = groundstatespectrum(rg3H |> crystalspectrum, perunitcellfillings=1)
visualize(rg3gs, usecontour=true, title="gsH-RG3")

rg3H = rg3[:courierzipper] * rg2H * rg3[:courierzipper]'
rg3gs = groundstatespectrum(rg3H|>crystalspectrum, perunitcellfillings=1)
visualize(rg3gs, usecontour=true, title="gsH-RG3-RG3space")

rg3det = quantummetricdeterminants(rg3gs)
visualize(rg3det|>crystalspectrum, usecontour=true, title="qgtdet-RG3-RG3space")

rg3qm = computequantummetric(rg3gs)
visualize(rg3qm[:gxx]|>crystalspectrum, usecontour=true, title="gxx-RG3-RG3space")
visualize(rg3qm[:gyy]|>crystalspectrum, usecontour=true, title="gyy-RG3-RG3space")

rg1[:courierzipper] * H * rg1[:courierzipper]' |> crystalspectrum |> visualize
rg1[:filledzipper] * H * rg1[:filledzipper]' |> crystalspectrum |> visualize
rg1[:emptyzipper] * H * rg1[:emptyzipper]' |> crystalspectrum |> visualize

courierH = rg1[:courierzipper] * H * rg1[:courierzipper]'
rg2[:blocker] * courierH * rg2[:blocker]' |> crystalspectrum |> visualize

rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' |> crystalspectrum |> visualize
rg2[:filledzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:filledzipper]' |> crystalspectrum |> visualize
rg2[:emptyzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:emptyzipper]' |> crystalspectrum |> visualize

courierH = rg2[:courierzipper] * courierH * rg2[:courierzipper]'
rg3[:blocker] * courierH * rg3[:blocker]' |> crystalspectrum |> visualize

rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' |> crystalspectrum |> visualize
rg3[:filledzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:filledzipper]' |> crystalspectrum |> visualize
rg3[:emptyzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:emptyzipper]' |> crystalspectrum |> visualize

courierH = rg3[:courierzipper] * courierH * rg3[:courierzipper]'
rg4[:blocker] * courierH * rg4[:blocker]' |> crystalspectrum |> visualize

rg4[:courierzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:courierzipper]' |> crystalspectrum |> visualize
rg4[:filledzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:filledzipper]' |> crystalspectrum |> visualize
rg4[:emptyzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:emptyzipper]' |> crystalspectrum |> visualize

filled = rg1[:filledzipper]' * rg1[:filledzipper]
empty = rg1[:emptyzipper]' * rg1[:emptyzipper]

filled1 = rg1[:courierzipper]' * rg2[:filledzipper]' * rg2[:filledzipper] * rg1[:courierzipper]
empty1 = rg1[:courierzipper]' * rg2[:emptyzipper]' * rg2[:emptyzipper] * rg1[:courierzipper]

filled2 = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:filledzipper]' * rg3[:filledzipper] * rg2[:courierzipper] * rg1[:courierzipper]
empty2 = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:emptyzipper]' * rg3[:emptyzipper] * rg2[:courierzipper] * rg1[:courierzipper]

filled3 = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:filledzipper]' * rg4[:filledzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]
empty3 = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:emptyzipper]' * rg4[:emptyzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]

core = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:courierzipper]' * rg4[:courierzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]

occ1 = momentumoccupations(filled) + momentumoccupations(empty)
occ2 = momentumoccupations(filled1) + momentumoccupations(empty1)
occ3 = momentumoccupations(filled2) + momentumoccupations(empty2)
occ4 = momentumoccupations(filled3) + momentumoccupations(empty3)
occcore = momentumoccupations(core)

occ1 |> crystalspectrum |> visualize
occ2 |> crystalspectrum |> visualize
occ3 |> crystalspectrum |> visualize
occ4 |> crystalspectrum |> visualize
occcore |> crystalspectrum |> visualize

visualize(occ1|>crystalspectrum, usecontour=true, title="occ1")
visualize(occ2|>crystalspectrum, usecontour=true, title="occ2")
visualize(occ3|>crystalspectrum, usecontour=true, title="occ3")
visualize(occ4|>crystalspectrum, usecontour=true, title="occ4")
visualize(occcore|>crystalspectrum, usecontour=true, title="occcore")

occ1 + occ2 + occ3 + occ4 + occcore |> crystalspectrum |> visualize

occ1 + occ2 + occ3 + occcore |> crystalspectrum |> visualize

frozenocc = momentumoccupations(rg1[:filledcorrelations] + rg1[:emptycorrelations]) |> crystalspectrum
frozenocc |> FockMap |> eigspech |> visualize

occ0 = momentumoccupations(C)
occ1 = momentumoccupations(rg1[:courierzipper]' * rg1[:correlations] * rg1[:courierzipper])
occ2 = momentumoccupations(rg1[:courierzipper]' * rg2[:courierzipper]' * rg2[:correlations] * rg2[:courierzipper] * rg1[:courierzipper]) - occ1
occ3 = momentumoccupations(rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg3[:correlations] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]) - occ2 - occ1

rg1C = rg1[:courierzipper]' * rg1[:correlations] * rg1[:courierzipper]
rg2C = rg1[:courierzipper]' * rg2[:courierzipper]' * rg2[:correlations] * rg2[:courierzipper] * rg1[:courierzipper]
rg3C = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg3[:correlations] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]

visualize(rg1C, colrange=1:512, rowrange=1:512)
visualize(rg2C, colrange=1:512, rowrange=1:512)
visualize(rg3C, colrange=1:512, rowrange=1:512)

occ0 |> crystalspectrum |> visualize
occ1 |> crystalspectrum |> visualize
occ2 |> crystalspectrum |> visualize
occ3 |> crystalspectrum |> visualize

occ1 + occ2 + occ3 |> crystalspectrum |> visualize
