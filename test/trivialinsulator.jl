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
c3 = c6^2

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [32, 32])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

t_a = ComplexF64(-0.3)
t_b = ComplexF64(-0.7)
t_n = ComplexF64(-0.5)
bonds::FockMap = bondmap([
    (m0, m0) => t_a,
    (m1, m1) => t_b,
    (m0, m1) => t_n,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => t_n,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => t_n])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize
groundstateprojector = groundstates |> crystalprojector

C = idmap(groundstateprojector.outspace) - groundstateprojector
correlationspectrum = C |> crystalspectrum
visualize(correlationspectrum, title="Physical Correlation")

correlations = C

crystalfock = correlations.outspace

scale = Scale([2 0; 0 2])
blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

blockedcrystal::Crystal = blockresult[:crystal]
blockedcorrelations::FockMap = blockresult[:correlations]

function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
    physicalnorm = m -> lineartransform(currentspace, m |> pos) |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
end

crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)

frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> origin, physicalmodes, 2.0)
frozenseedingregion::Subset{Offset} = Subset(m |> pos for m in frozenseedingmodes)
visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
frozenseedingfock::FockSpace = FockSpace(frozenseedingmodes)

globaldistiller = globaldistillerhamiltonian(
    correlations=blockresult[:correlations],
    restrictspace=frozenseedingfock,
    localisometryselectionstrategy=frozenselectionbycount(6),
    symmetry=c3)

globaldistillerspectrum = globaldistiller |> crystalspectrum
globaldistillerspectrum |> visualize

distillresult = distillation(globaldistillerspectrum, :courier => v -> abs(v) < 1e-5, :filled => v -> v > 1e-5, :empty => v -> v < -1e-5)

filledseedingcenter::Offset = (blockedmodes |> getspace) & [2/3, 1/3]
courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
courierseedingregion::Subset{Offset} = Subset(m |> pos for m in courierseedingmodes)
visualize(courierseedingregion, courierseedingcenter |> Subset, title="Courier Seeding Region", visualspace=euclidean(RealSpace, 2))
courierseedingfock::FockSpace = FockSpace(courierseedingregion)

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
couriercorrelations = purifiedcorrelationspectrum |> FockMap

commutation(c6 * couriercorrelations.outspace, couriercorrelations) |> maximum

