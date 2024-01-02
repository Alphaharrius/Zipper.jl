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
crystal = Crystal(unitcell, [24, 24])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

tₙ = -1 + 0im
tₕ = 0.1im

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:offset => [-1, 0] ∈ triangular)) => tₙ,
    (m0, m1 |> setattr(:offset => [0, 1] ∈ triangular)) => tₙ]

haldane = [
    (m0, setattr(m0, :offset => Point([1, 1], triangular))) => tₕ,
    (m0, setattr(m0, :offset => Point([-1, 0], triangular))) => tₕ,
    (m0, setattr(m0, :offset => Point([0, -1], triangular))) => tₕ,
    (m1, setattr(m1, :offset => Point([1, 1], triangular))) => -tₕ,
    (m1, setattr(m1, :offset => Point([-1, 0], triangular))) => -tₕ,
    (m1, setattr(m1, :offset => Point([0, -1], triangular))) => -tₕ]

bonds::FockMap = bondmap([nearestneighbor..., haldane...])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize

groundstateprojector = groundstates |> crystalprojector
C = idmap(groundstateprojector.outspace) - groundstateprojector

CC = CrystalFockMap(C)

function _globaldistillerhamiltonian(;
    correlations::FockMap, restrictspace::FockSpace,
    localisometryselectionstrategy::Function, manualeigenenergies::Dict{Symbol, <:Number} = Dict(:filled => -1, :empty => 1))

    localisometries::Dict{Symbol} = localfrozenisometries(correlations, restrictspace, selectionstrategy=localisometryselectionstrategy)
    crystalprojectors::Dict{Symbol, FockMap} = Dict(
        name => crystalprojector(localisometry=localisometries[name], crystalfock=correlations.inspace)
        for (name, isometry) in localisometries)
    return crystalprojectors
    return reduce(+, manualeigenenergies[name] * crystalprojector for (name, crystalprojector) in crystalprojectors)
end

    correlations = C
    crystalfock = correlations.outspace

    scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)
    @time blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

    blockedcrystal::Crystal = blockresult[:crystal]
    blockedcorrelations::FockMap = blockresult[:correlations]
    blocker = blockresult[:transformer]

    function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
        currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
        physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
        return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
    end

    crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
    samplepoints::Subset{Offset} = crystalpoints + c6^2 * crystalpoints + c6^4 * crystalpoints
    blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
    physicalmodes::Subset{Mode} = spanoffset(blockedmodes, samplepoints)


    frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 2)
    # visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
    frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)

    @time globaldistiller = _globaldistillerhamiltonian(
        correlations=blockresult[:correlations],
        restrictspace=frozenseedingfock,
        localisometryselectionstrategy=frozenselectionbycount(3))

    globaldistiller[:filled] + globaldistiller[:empty]
    k0 = (0.5,0.5)∈convert(MomentumSpace, blockedcrystal|>getspace)
    b = globaldistiller[:empty]
    bsubmaps::Dict{Tuple{Momentum, Momentum}, FockMap} = Dict((k, b.momentummappings[k]) => b.subfockmaps[k] for k in b.outbz)
    kk = Dict(k => k for k in b.outbz)[k0]

    Base.:hash(v::Vector{Momentum}) = reduce(+, hash(e) for e in v)

    (kk,kk) == (k0,k0)
    (kk,kk) |> hash
    (k0,k0)|>hash
    [kk,kk]|>hash
    [k0,k0]|>hash    



inplaceadd(a::FockMap, b::FockMap)::FockMap = a + FockMap(b, outspace=a|>getoutspace, inspace=a|>getinspace, performpermute=false)

function computequantumgeometrictensor2d(state::CrystalSpectrum)
    statevectors::Dict{Momentum, FockMap} = state|>geteigenvectors
    Δkx = (1 / size(crystal)[1], 0) ∈ kspace
    Δky = (0, 1 / size(crystal)[2]) ∈ kspace

    kprojectors::Dict{Momentum, FockMap} = Dict(k => kstate * kstate' for (k, kstate) in statevectors)
    kcorrelations::Dict{Momentum, FockMap} = Dict(k => idmap(projector|>getoutspace) - projector for (k, projector) in kprojectors)
    ΔCkxs::Dict{Momentum, FockMap} = Dict(k => inplaceadd(-kcorr, kcorrelations[k + Δkx |> basispoint]) for (k, kcorr) in kcorrelations)
    ΔCkys::Dict{Momentum, FockMap} = Dict(k => inplaceadd(-kcorr, kcorrelations[k + Δky |> basispoint]) for (k, kcorr) in kcorrelations)
    
    function kqgtensor(k::Momentum)::FockMap
        gxx = statevectors[k]' * (ΔCkxs[k]' * kcorrelations[k] * ΔCkxs[k]) * statevectors[k]
        gxy = statevectors[k]' * (ΔCkxs[k]' * kcorrelations[k] * ΔCkys[k]) * statevectors[k]
        gyx = statevectors[k]' * (ΔCkys[k]' * kcorrelations[k] * ΔCkxs[k]) * statevectors[k]
        gyy = statevectors[k]' * (ΔCkys[k]' * kcorrelations[k] * ΔCkys[k]) * statevectors[k]

        xmode::Mode = Mode(:offset=>k, :axis=>(:x))
        ymode::Mode = Mode(:offset=>k, :axis=>(:y))

        fockspace::FockSpace = FockSpace([xmode, ymode])
        return FockMap(fockspace, fockspace, [rep(gxx)[1,1] rep(gxy)[1,1]; rep(gyx)[1,1] rep(gyy)[1,1]])
    end

    return Dict(k=>k|>kqgtensor for (k, _) in statevectors)
end

Base.:real(fockmap::FockMap)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, fockmap|>rep|>real)
Base.:imag(fockmap::FockMap)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, fockmap|>rep|>imag)

kqgtensors = computequantumgeometrictensor2d(groundstates)
gxx = directsum(tensor[1,2] for (k, tensor) in kqgtensors)
barecrystal = Crystal(triangular|>getorigin, crystal|>size)
barecrystalfock = FockSpace(gxx|>getoutspace, reflected=barecrystal)
gxx = FockMap(gxx, inspace=barecrystalfock, outspace=barecrystalfock, performpermute=false)|>real
gxx|>crystalspectrum|>visualize

bc = directsum(tensor[1,2]|>imag for (k, tensor) in kqgtensors)
barecrystalfock = FockSpace(bc|>getoutspace, reflected=barecrystal)
berrycurvature = 2 * FockMap(bc, inspace=barecrystalfock, outspace=barecrystalfock, performpermute=false)

(berrycurvature|>rep|>sum) / (2π)

berrycurvature|>crystalspectrum|>visualize

groundstateeigenvectors = groundstates |> geteigenvectors
δkx = (1 / size(crystal)[1], 0) ∈ kspace
δky = (0, 1 / size(crystal)[2]) ∈ kspace

δIxs = Dict(k => FockMap(groundstateeigenvectors[k + δkx |> basispoint], inspace=I|>getinspace, outspace=I|>getoutspace, performpermute=false) - I for (k, I) in groundstateeigenvectors)
δIys = Dict(k => FockMap(groundstateeigenvectors[k + δky |> basispoint], inspace=I|>getinspace, outspace=I|>getoutspace, performpermute=false) - I for (k, I) in groundstateeigenvectors)
δIx = directsum(I for (_, I) in δIxs)
δIy = directsum(I for (_, I) in δIys)

Base.:real(fockmap::FockMap)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, fockmap|>rep|>real)

qgtcrystal = Crystal(triangular |> getorigin, crystal|>size)
qgt = δIx' * groundstateprojector * δIy
qgtoutspace = FockSpace(qmt|>getoutspace, reflected=qmtcrystal)
quantumgeometrictensor = idmap(qgtoutspace) - FockMap(qmt, inspace=qmtoutspace, outspace=qmtoutspace)
quantummetrictensor = quantumgeometrictensor |> real
quantummetrictensor |> crystalspectrum |> visualize

C |> crystalspectrum |> visualize

entanglemententropy(C |> crystalspectrum)

function zer(correlations::FockMap)        
    crystalfock = correlations.outspace

    scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)
    @time blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

    blockedcrystal::Crystal = blockresult[:crystal]
    blockedcorrelations::FockMap = blockresult[:correlations]
    blocker = blockresult[:transformer]

    function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
        currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
        physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
        return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
    end

    crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
    samplepoints::Subset{Offset} = crystalpoints + c6^2 * crystalpoints + c6^4 * crystalpoints
    blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
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

    return Dict(
        :blocker => blocker,
        :correlations => couriercorrelations,
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
        :entanglemententropy => entanglemententropy(couriercorrelationspectrum))
end

rg1 = zer(C)
rg2 = zer(rg1[:correlations])
(rg2[:unpurecorrelations] |> crystalspectrum |> entanglemententropy) / (rg1[:correlations]|>getoutspace|>getcrystal|>vol)
rg3 = zer(rg2[:correlations])
(rg3[:unpurecorrelations] |> crystalspectrum |> entanglemententropy) / (rg2[:correlations]|>getoutspace|>getcrystal|>vol)
rg4 = zer(rg3[:correlations])

visualize(frozenseedingregion, Subset(m |> getpos for m in rg2[:frozenseedingfock]))


function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = blockedcrystal |> getspace |> orthogonalspace
    physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
end

blockedcrystal = rg1[:correlations]|>getoutspace|>getcrystal
crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
samplepoints::Subset{Offset} = crystalpoints + c6^2 * crystalpoints + c6^4 * crystalpoints
blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, samplepoints)


frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 8)
frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)


regionalrestriction(rg1[:filledisometry], frozenseedingfock) |> visualize

rg1[:courierzipper]' * rg2[:blocker]' * rg2[:filledisometry]

regionalrestriction(rg1[:courierisometry] * rg2[:blocker]' * rg2[:filledisometry], frozenseedingfock) |> visualize
rg1[:courierisometry] * rg2[:courierzipper]' * rg3[:blocker]' * rg3[:filledisometry]
regionalrestriction(rg1[:courierisometry] * rg2[:courierzipper]' * rg3[:blocker]' * rg3[:filledisometry], frozenseedingfock) |> visualize

rg3[:correlations] |> crystalspectrum |> visualize

regionalrestriction(rg3[:blocker]' * rg3[:filledisometry], rg2[:frozenseedingfock]) |> visualize

rg1[:entanglemententropy]

rg1[:globaldistiller] |> crystalspectrum |> visualize
rg1[:filledisometry] |> getinspace |> unitcellfock |> modeattrs

rg1[:filledcorrelations] |> crystalspectrum |> visualize

H = energyspectrum |> FockMap
energyspectrum |> visualize

rg1[:blocker] * H * rg1[:blocker]' |> crystalspectrum |> visualize

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

occ1 + occ2 + occ3 + occ4 + occcore |> crystalspectrum |> visualize

occ1 + occ2 + occ3 + occcore |> crystalspectrum |> visualize

frozenocc = momentumoccupations(rg1[:filledcorrelations] + rg1[:emptycorrelations]) |> crystalspectrum
frozenocc |> FockMap |> eigspech |> visualize

occ0 = momentumoccupations(C)
occ1 = momentumoccupations(rg1[:courierzipper]' * rg1[:correlations] * rg1[:courierzipper])
occ2 = momentumoccupations(rg1[:courierzipper]' * rg2[:courierzipper]' * rg2[:correlations] * rg2[:courierzipper] * rg1[:courierzipper]) - occ1
occ3 = momentumoccupations(rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg3[:correlations] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]) - occ2 - occ1

occ0 |> crystalspectrum |> visualize
occ1 |> crystalspectrum |> visualize
occ2 |> crystalspectrum |> visualize
occ3 |> crystalspectrum |> visualize

occ1 + occ2 + occ3 |> crystalspectrum |> visualize
