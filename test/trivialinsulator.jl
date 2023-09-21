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
crystal = Crystal(unitcell, [64, 64])
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

bonds |> visualize
H::FockMap = hamiltonian(crystal, bonds)
C::FockMap = groundstatecorrelations(H)

energyspectrum = crystalspectrum(H)
correlationspectrum = C |> crystalspectrum

visualize(energyspectrum, title="Physical Hamiltonian")
visualize(correlationspectrum, title="Physical Correlation")

correlations = C

crystalfock = correlations.outspace

scale = Scale([2 0; 0 2])
blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

blockedcrystal::Crystal = blockresult[:crystal]
blockedcorrelations::FockMap = blockresult[:correlations]

blockedhamiltonian = blockresult[:transformer] * H * blockresult[:transformer]'
blockedhamiltonian |> crystalspectrum |> visualize

commutation(c3 * blockedhamiltonian.outspace, blockedhamiltonian) |> maximum

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
    localisometryselectionstrategy=frozenselectionbycount(6),
    symmetry=c3)

commutation(c3 * globaldistiller.outspace, globaldistiller) |> maximum

visualize(globaldistiller |> crystalspectrum, title="Global Distiller")
distillresult = distillation(globaldistiller, courierenergythreshold=1e-5)

filledseedingcenter::Position = (blockedmodes |> getspace) & [2/3, 1/3]
courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
courierseedingregion::Subset{Position} = Subset(m |> pos for m in courierseedingmodes)
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

