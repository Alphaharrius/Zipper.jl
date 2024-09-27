using Zipper
using LinearAlgebra,Plots

triangular = RealSpace([1 0; 1/2 sqrt(3)/2]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
# spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

unitcell = Subset(pa, pb)
visualize(unitcell)
crystal = Crystal(unitcell, [32, 32])
space = crystal|>getspace
# reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1 = members(modes)

t_a = 0 + 0im
t_b = -0 + 0im
tₙ = -1 + 0im

onsite = [
    (m1, m1) => t_b,
    (m0, m0) => t_a
]

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:r => [-1, 0] ∈ triangular)) => tₙ,
    (m0, m1 |> setattr(:r => [0, 1] ∈ triangular)) => tₙ]

bonds::FockMap = bondmap([onsite..., nearestneighbor...])

energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)

groupedeigenvalues = groundstatespectrum(energyspectrum, perunitcellfillings=1)

length(groupedeigenvalues)

for (eval, modeset) in groupedeigenvalues
    println(eval)
end