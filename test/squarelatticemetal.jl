using Revise, LinearAlgebra
using Zipper

square = euclidean(RealSpace, 2)
point = [1/2, 1/2] ∈ square
spatialsnappingcalibration([point])

kspace = convert(MomentumSpace, square)

c4 = pointgrouptransform([0 -1; 1 0])

unitcell = Subset(point)
crystal = Crystal(unitcell, [12, 12])
reciprocalhashcalibration(crystal.sizes)

m = quantize(:pos, unitcell, 1) |> first

subsetunion((Subset(m), Subset(m)))

tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m, m |> setattr(:offset => [1, 0] ∈ square)) => tₙ,
    (m, m |> setattr(:offset => [0, 1] ∈ square)) => tₙ])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=0.5)
groundstates |> visualize
groundstateprojector = groundstates |> crystalprojector
groundstateprojector |> crystalspectrum |> visualize

C = idmap(groundstateprojector.outspace) - groundstateprojector

correlations = C

correlations |> crystalspectrum |> visualize

crystalfock = correlations |> getoutspace

scale = Scale([2 0; 0 2])
blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

H = energyspectrum |> FockMap
blockresult[:transformer] * H * blockresult[:transformer]' |> crystalspectrum |> visualize

function crosssectionalrestriction(normalvector::Offset, center::Offset, radius::Real)

end
