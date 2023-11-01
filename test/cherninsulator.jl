using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Revise
using Zipper

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = triangular & [1/3, 2/3]
pb = triangular & [2/3, 1/3]
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [96, 96])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

tₙ = -1 + 0im
tₕ = 0.1im

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:offset => triangular & [-1, 0])) => tₙ,
    (m0, m1 |> setattr(:offset => triangular & [0, 1])) => tₙ]

haldane = [
    (m0, m0 |> setattr(:offset => triangular & [1, 1])) => tₕ,
    (m0, m0 |> setattr(:offset => triangular & [-1, 0])) => tₕ,
    (m0, m0 |> setattr(:offset => triangular & [0, -1])) => tₕ,
    (m1, m1 |> setattr(:offset => triangular & [1, 1])) => -tₕ,
    (m1, m1 |> setattr(:offset => triangular & [-1, 0])) => -tₕ,
    (m1, m1 |> setattr(:offset => triangular & [0, -1])) => -tₕ]

bonds::FockMap = bondmap([nearestneighbor..., haldane...])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize
groundstateprojector = groundstates |> crystalprojector

C = idmap(groundstateprojector.outspace) - groundstateprojector

C |> crystalspectrum |> visualize

function entanglemententropy(correlationspectrum::CrystalSpectrum)
    sumresult = 0
    for (_, v) in correlationspectrum |> geteigenvalues
        if isapprox(v, 0, atol=1e-7) || isapprox(v, 1, atol=1e-7)
            continue
        end
        sumresult += v * log(v) + (1 - v) * log(1 - v)
    end
    return sumresult
end

entanglemententropy(C |> crystalspectrum)

function zer(correlations::FockMap)
    crystalfock = correlations.outspace

    scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)
    blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

    blockedcrystal::Crystal = blockresult[:crystal]
    blockedcorrelations::FockMap = blockresult[:correlations]
    blocker = blockresult[:transformer]

    function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
        currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
        physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
        return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
    end

    crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
    blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
    physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)

    frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 2.0)
    frozenseedingregion::Subset{Offset} = Subset(m |> getpos for m in frozenseedingmodes)
    # visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
    frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)

    globaldistiller = globaldistillerhamiltonian(
        correlations=blockresult[:correlations],
        restrictspace=frozenseedingfock,
        localisometryselectionstrategy=frozenselectionbycount(3),
        symmetry=c6)

    globaldistillerspectrum = globaldistiller |> crystalspectrum
    # visualize(globaldistillerspectrum, title="Global Distiller")

    distillresult = distillation(globaldistillerspectrum, :courier => v -> abs(v) < 1e-5, :filled => v -> v > 1e-5, :empty => v -> v < -1e-5)

    courierseedingcenter::Offset = (blockedmodes |> getspace) & [2/3, 1/3]
    courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
    courierseedingregion::Subset{Offset} = Subset(m |> getpos for m in courierseedingmodes)
    # visualize(courierseedingregion, courierseedingcenter |> Subset, title="Courier Seeding Region", visualspace=euclidean(RealSpace, 2))
    courierseedingfock::FockSpace{Region} = FockSpace{Region}(courierseedingmodes)

    c3 = c6^2 |> recenter(courierseedingcenter)

    blockedcourierprojector = distillresult[:courier] |> crystalprojector
    blockedcouriercorrelation = idmap(blockedcourierprojector.outspace, blockedcourierprojector.outspace) - blockedcourierprojector

    localcourierseed = findlocalspstates(statecorrelations=blockedcouriercorrelation, regionfock=courierseedingfock, symmetry=c3, spectrumextractpredicate=v -> v < 5e-2)[1]
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
    couriercorrelations = purifiedcorrelationspectrum |> FockMap

    blockedfilledprojector = distillresult[:filled] |> crystalprojector
    blockedfilledcorrelation = idmap(blockedfilledprojector.outspace, blockedfilledprojector.outspace) - blockedfilledprojector
    filledseed = findlocalspstates(statecorrelations=blockedfilledcorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalfilledseeds = crystalisometries(localisometry=filledseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannierfilledisometry = wannierprojection(
        crystalisometries=distillresult[:filled].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledseeds)

    blockedemptyprojector = distillresult[:empty] |> crystalprojector
    blockedemptycorrelation = idmap(blockedemptyprojector.outspace, blockedemptyprojector.outspace) - blockedemptyprojector
    emptyseed = findlocalspstates(statecorrelations=blockedemptycorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalemptyseeds = crystalisometries(localisometry=emptyseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannieremptyisometry = wannierprojection(
        crystalisometries=distillresult[:empty].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyseeds)

    filledcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
    emptycorrelations = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry

    return Dict(
        :blocker => blocker,
        :correlations => couriercorrelations, 
        :courierzipper => wanniercourierisometry' * blocker, 
        :filledzipper => wannierfilledisometry' * blocker,
        :emptyzipper => wannieremptyisometry' * blocker,
        :globaldistiller => globaldistiller, 
        :filledcorrelations => filledcorrelations, 
        :emptycorrelations => emptycorrelations,
        :entanglemententropy => entanglemententropy(couriercorrelationspectrum))
end

rg1 = zer(C)
rg2 = zer(rg1[:correlations])
rg3 = zer(rg2[:correlations])
rg4 = zer(rg3[:correlations])

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

function momentumoccupations(correlations::FockMap)::FockMap
    kcorrelations::Base.Generator = correlations |> crystalsubmaps
    crystal::Crystal = correlations |> getoutspace |> getcrystal
    center::Offset = crystal |> getspace |> getorigin
    tracecrystal::Crystal = Crystal(center |> Subset, crystal |> size)
    mode::Mode = Mode(:pos => center)

    function tracing(k::Momentum, corr::FockMap)::FockMap
        space::FockSpace = mode |> setattr(:offset => k) |> FockSpace
        return FockMap(space, space, [corr |> tr][:, :] |> SparseMatrixCSC) / dimension(corr |> getoutspace)
    end

    occupations::FockMap = directsum(tracing(k, corr) for (k, corr) in kcorrelations)
    fockspace::FockSpace = FockSpace(occupations |> getoutspace, reflected=tracecrystal)
    return FockMap(occupations, inspace=fockspace, outspace=fockspace)
end

using Compat
function Zipper.visualize(spectrum::CrystalSpectrum{2}; title="", toppadding::Bool = true)
    kspectrum::Dict{Momentum} = Dict(k => ([spectrum.eigenvalues[m] for m in modes] |> sort) for (k, modes) in spectrum.eigenmodes)
    mesh::Matrix{Momentum} = spectrum.crystal |> brillouinmesh
    plottingdata::Matrix{Vector} = map(k -> haskey(kspectrum, k) ? kspectrum[k] : [], mesh)
    bandcount::Integer = map(v -> v |> length, plottingdata) |> maximum
    function padding(v::Vector)::Vector
        return toppadding ? vcat(v, repeat([NaN], bandcount - length(v))) : vcat(repeat([NaN], bandcount - length(v)), v)
    end
    paddeddata::Matrix{Vector} = map(v -> v |> padding, plottingdata)
    plottingspectrum::Array = paddeddata |> stack
    layout::Layout = Layout(title=title)
    plot([contour(z=plottingspectrum[n, :, :]) for n in axes(plottingspectrum, 1)], layout)
end

function Zipper.visualize(spectrum::CrystalSpectrum{2}; title="", toppadding::Bool = true)
    kspectrum::Dict{Momentum} = Dict(k => ([spectrum.eigenvalues[m] for m in modes] |> sort) for (k, modes) in spectrum.eigenmodes)
    mesh::Matrix{Momentum} = spectrum.crystal |> brillouinmesh
    plottingdata::Matrix{Vector} = map(k -> haskey(kspectrum, k) ? kspectrum[k] : [], mesh)
    bandcount::Integer = map(v -> v |> length, plottingdata) |> maximum
    function padding(v::Vector)::Vector
        return toppadding ? vcat(v, repeat([NaN], bandcount - length(v))) : vcat(repeat([NaN], bandcount - length(v)), v)
    end
    paddeddata::Matrix{Vector} = map(v -> v |> padding, plottingdata)
    plottingspectrum::Array = paddeddata |> stack
    layout::Layout = Layout(title=title)
    plot([surface(z=plottingspectrum[n, :, :]) for n in axes(plottingspectrum, 1)], layout)
end

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
