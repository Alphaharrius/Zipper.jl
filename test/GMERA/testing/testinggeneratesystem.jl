using LinearAlgebra
using Zipper, Plots

# usecrystaldensemap()
setmaxthreads(Threads.nthreads())

triangular = RealSpace([1 0; 1/2 sqrt(3)/2]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 0] ∈ triangular
pb = [0, 1/3] ∈ triangular
pc = [0, 2/3] ∈ triangular
pd = [1/3, 2/3] ∈ triangular
pe = [2/3, 1/3] ∈ triangular
pf = [2/3, 0] ∈ triangular

pg = (pa + pb + pc + pd + pe + pf) / 6
pg
spatialsnappingcalibration((pa, pb, pc, pd, pe, pf, pg))

unitcell = Subset(pa, pb, pc, pd, pe, pf)
visualize(unitcell)
crystal = Crystal(unitcell, [32, 32])
space = crystal|>getspace
modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1, m2, m3, m4, m5 = members(modes)
reciprocalhashcalibration(crystal.sizes)

t_a = 0 + 0im
t_b = 0 + 0im
tₙ = -1 + 0im

onsite = [
    (m0, m0) => t_a,
    (m1, m1) => t_b,
    (m2, m2) => t_a,
    (m3, m3) => t_b,
    (m4, m4) => t_a,
    (m5, m5) => t_b
]

nearestneighbor = [
        (m1, m0) => tₙ,
        (m2, m5|> setattr(:r => [-1, 1] ∈ triangular)) => -tₙ,
        (m1, m2) => tₙ,
        (m3, m2) => tₙ,
        (m4, m1|> setattr(:r => [1, 0] ∈ triangular)) => -tₙ,
        (m3, m4) => tₙ,
        (m5, m4) => tₙ,
        (m0, m3|> setattr(:r => [0, -1] ∈ triangular)) => -tₙ,
        (m5, m0) => tₙ]

bonds::FockMap = bondmap([onsite..., nearestneighbor...])
energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=3)
groundstateprojector = groundstates|>crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector
H = CrystalFockMap(energyspectrum)


seedingfock::RegionFock = quantize(unitcell, 1)
rscorrelations = regioncorrelations(correlations,seedingfock)
rscorrelations[[m for m in seedingfock][6],[m for m in seedingfock][6]]|>eigspech|>geteigenvalues
