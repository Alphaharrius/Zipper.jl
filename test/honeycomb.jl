include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/transformations.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")
include("../src/quantumtransformations.jl")
include("../src/renormalization.jl")

using PlotlyJS, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using ..Spaces, ..Geometries, ..Quantum, ..Transformations, ..Plotting, ..QuantumTransformations, ..Physical, ..Renormalization

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = triangular & [1/3, 2/3]
pb = triangular & [2/3, 1/3]
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [48, 48])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:offset => triangular & [-1, 0])) => tₙ,
    (m0, m1 |> setattr(:offset => triangular & [0, 1])) => tₙ])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstateprojector = groundstates |> crystalprojector

C = idmap(groundstateprojector.outspace) - groundstateprojector

correlations = C

crystalfock = correlations.outspace

scale = Scale([2 0; 0 2])
blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

blockedcrystal::Crystal = blockresult[:crystal]
blockedcorrelations::FockMap = blockresult[:correlations]

function circularregionmodes(origin::Position, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
    physicalnorm = m -> lineartransform(currentspace, m |> pos) |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
end

crystalpoints::Subset{Position} = latticepoints(blockedcrystal)
blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)

frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> origin, physicalmodes, 2.0)
frozenseedingregion::Subset{Position} = Subset(m |> pos for m in frozenseedingmodes)
visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
frozenseedingfock::FockSpace = FockSpace(frozenseedingmodes)

globaldistiller = globaldistillerhamiltonian(
    correlations=blockresult[:correlations],
    restrictspace=frozenseedingfock,
    localisometryselectionstrategy=frozenselectionbycount(3),
    symmetry=c6)

visualize(globaldistiller |> crystalspectrum, title="Global Distiller")
distillresult = distillation(globaldistiller, courierenergythreshold=1e-5)

courierseedingcenter::Position = (blockedmodes |> getspace) & [2/3, 1/3]
courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
courierseedingregion::Subset{Position} = Subset(m |> pos for m in courierseedingmodes)
visualize(courierseedingregion, courierseedingcenter |> Subset, title="Courier Seeding Region", visualspace=euclidean(RealSpace, 2))
courierseedingfock::FockSpace = FockSpace(courierseedingmodes)

c3 = c6^2 |> recenter(courierseedingcenter)

blockedcourierprojector = distillresult[:courier] |> crystalprojector
blockedcouriercorrelation = idmap(blockedcourierprojector.outspace, blockedcourierprojector.outspace) - blockedcourierprojector

localcourierseed = [regionalwannierseeding(blockedcouriercorrelation, courierseedingfock, symmetry=c3)...][1]
fullcourierseed = localcourierseed + (c6 * localcourierseed.outspace) * localcourierseed * (c6 * localcourierseed.inspace)'

crystalcourierseeds = crystalisometries(localisometry=fullcourierseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)

wanniercourierisometry = wannierprojection(
    crystalisometries=distillresult[:courier].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds)

couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
couriercorrelations |> crystalspectrum |> visualize
couriercorrelationspectrum = couriercorrelations |> crystalspectrum
purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
purifiedcorrelationspectrum |> visualize
couriercorrelations = purifiedcorrelationspectrum |> FockMap

commutation(c6 * couriercorrelations.outspace, couriercorrelations) |> maximum

blockedfilledprojector = distillresult[:filled] |> crystalprojector
blockedfilledcorrelation = idmap(blockedfilledprojector.outspace, blockedfilledprojector.outspace) - blockedfilledprojector
filledseed = [regionalwannierseeding(blockedfilledcorrelation, frozenseedingfock, symmetry=c6, seedsgroupingprecision=1e-3)...][1]

crystalfilledseeds = crystalisometries(localisometry=filledseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
wannierfilledisometry = wannierprojection(
    crystalisometries=distillresult[:filled].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledseeds)

blockedemptyprojector = distillresult[:empty] |> crystalprojector
blockedemptycorrelation = idmap(blockedemptyprojector.outspace, blockedemptyprojector.outspace) - blockedemptyprojector
emptyseed = [regionalwannierseeding(blockedemptycorrelation, frozenseedingfock, symmetry=c6, seedsgroupingprecision=1e-3)...][1]

crystalemptyseeds = crystalisometries(localisometry=emptyseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
wannieremptyisometry = wannierprojection(
    crystalisometries=distillresult[:empty].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyseeds)

filledC = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
filledC |> crystalspectrum |> visualize
emptyC = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry
emptyC |> crystalspectrum |> visualize

filledP = wannierfilledisometry * wannierfilledisometry'
courierP = wanniercourierisometry * wanniercourierisometry'

filledC = idmap(filledP |> getinspace) - filledP
courierC = idmap(courierP |> getinspace) - courierP

evals = regioncorrelations(filledC, frozenseedingfock) |> eigvalsh
plot(scatter(y=[v for (_, v) in evals] |> sort, mode="markers"))

somemodes = spanoffset(blockedhamiltonian |> getoutspace |> unitcellfock |> orderedmodes, blockedcrystal |> latticepoints)

Ft = fourier(filledC.outspace, somemodes |> FockSpace)

realfilledC = Ft' * filledC * Ft

U = regioncorrelations(filledC, frozenseedingfock) |> eigvecsh

M = idmap(somemodes |> FockSpace) - idmap(frozenseedingfock) + FockMap(U, inspace=U |> getoutspace, performpermute=false)

R = M * realfilledC * M'

visualize(R - realfilledC, rowrange=1:64, colrange=1:64)
R - realfilledC |> maximum
courierR = restrict(R, courierseedingfock, courierseedingfock)
plot(scatter(y=[v for (_, v) in courierR |> eigvalsh] |> sort, mode="markers"))

translatedfrozenmodes = circularregionmodes(triangular & [1, 0], physicalmodes, 2.0)
visualize(frozenseedingregion, Subset(m |> pos for m in translatedfrozenmodes), visualspace=euclidean(RealSpace, 2))

TfrozenR = restrict(R, FockSpace(translatedfrozenmodes), FockSpace(translatedfrozenmodes))
plot(scatter(y=[v for (_, v) in TfrozenR |> eigvalsh] |> sort, mode="markers"))

frozenR = restrict(R, frozenseedingfock, frozenseedingfock)
plot(scatter(y=[v for (_, v) in frozenR |> eigvalsh] |> sort, mode="markers"))

localcourierinfilled = regioncorrelations(courierC, courierseedingfock)
plot(scatter(y=[v for (_, v) in localcourierinfilled |> eigvalsh] |> sort, mode="markers"))
