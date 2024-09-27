using Zipper
using LinearAlgebra,Plots

# function basispoint(point::Point)::Point
#     rationalized::Vector = [hashablereal(v) for v in point |> vec]
#     return Point([mod(v |> numerator, v |> denominator) / denominator(v) for v in rationalized], point |> getspace)
# end

# triangular = RealSpace([1 0; 1/2 sqrt(3)/2]')
# kspace = convert(MomentumSpace, triangular)

# # pa = [1/3, 0] ∈ triangular
# # pb = [0, 1/3] ∈ triangular
# # pc = [0, 2/3] ∈ triangular
# # pd = [1/3, 2/3] ∈ triangular
# # pe = [2/3, 1/3] ∈ triangular
# # pf = [2/3, 0] ∈ triangular

# pa = [1/3, 0] ∈ triangular
# pb = [0, 1/3] ∈ triangular
# pc = [-1/3, 1/3] ∈ triangular
# pd = [-1/3, 0] ∈ triangular
# pe = [0, -1/3] ∈ triangular
# pf = [1/3, -1/3] ∈ triangular

# pg = (pa + pb + pc + pd + pe + pf) / 6
# spatialsnappingcalibration((pa, pb, pc, pd, pe, pf, pg))

# unitcell = Subset(pa, pb, pc, pd, pe, pf)
# crystal = Crystal(unitcell, [systemsize, systemsize])

# reciprocalhashcalibration(crystal.sizes)
    
# modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
# pc|>basispoint

# [m|>getattr(:r) for m in modes]
# m0, m1, m2, m3, m4, m5 = members(modes)

# onsite = [
#     (m0, m0) => t_a,
#     (m1, m1) => t_b,
#     (m2, m2) => t_a,
#     (m3, m3) => t_b,
#     (m4, m4) => t_a,
#     (m5, m5) => t_b
# ]

triangular = RealSpace([1 0; 1/2 sqrt(3)/2]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 0] ∈ triangular
pb = [0, 1/3] ∈ triangular
pc = [0, 2/3] ∈ triangular
pd = [1/3, 2/3] ∈ triangular
pe = [2/3, 1/3] ∈ triangular
pf = [2/3, 0] ∈ triangular

pg = (pa + pb + pc + pd + pe + pf) / 6
spatialsnappingcalibration((pa, pb, pc, pd, pe, pf, pg))

unitcell = Subset(pa, pb, pc, pd, pe, pf)
crystal = Crystal(unitcell, [systemsize, systemsize])
reciprocalhashcalibration(crystal.sizes)

unitcell|>visualize
    
modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1, m2, m3, m4, m5 = members(modes)

t_a = 0.00001+0.0im
t_b = -0.00001+0.0im
tₙ = -1+0.0im
tₕ = 0

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
        (m1, m4|> setattr(:r => [-1, 0] ∈ triangular)) => tₙ,
        (m1, m2) => tₙ,
        (m3, m2) => tₙ,
        (m3, m0|> setattr(:r => [0, 1] ∈ triangular)) => tₙ,
        (m3, m4) => tₙ,
        (m5, m4) => tₙ,
        (m5, m2|> setattr(:r => [1, -1] ∈ triangular)) => tₙ,
        (m5, m0) => tₙ]

haldaneterm = [
        (m1, m3|> setattr(:r => [-1, 0] ∈ triangular)) => tₕ,
        (m1, m3) => tₕ,
        (m1, m3|> setattr(:r => [0, -1] ∈ triangular)) => tₕ,
        (m3, m5|> setattr(:r => [-1, 1] ∈ triangular)) => tₕ,
        (m3, m5|> setattr(:r => [0, 1] ∈ triangular)) => tₕ,
        (m3, m5) => tₕ,
        (m5, m1) => tₕ,
        (m5, m1|> setattr(:r => [1, 0] ∈ triangular)) => tₕ,
        (m5, m1|> setattr(:r => [1, -1] ∈ triangular)) => tₕ,
        (m0, m4|> setattr(:r => [-1, 0] ∈ triangular)) => -tₕ,
        (m0, m4) => -tₕ,
        (m0, m4|> setattr(:r => [0, -1] ∈ triangular)) => -tₕ,
        (m2, m0|> setattr(:r => [-1, 1] ∈ triangular)) => -tₕ,
        (m2, m0|> setattr(:r => [0, 1] ∈ triangular)) => -tₕ,
        (m2, m0) => -tₕ,
        (m4, m2) => -tₕ,
        (m4, m2|> setattr(:r => [1, 0] ∈ triangular)) => -tₕ,
        (m4, m2|> setattr(:r => [1, -1] ∈ triangular)) => -tₕ]

bonds::FockMap = bondmap([onsite..., nearestneighbor...,haldaneterm...])

energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
# groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=3)
groupedeigenvalues = groundstatespectrum(energyspectrum, perunitcellfillings=3)

length(groupedeigenvalues)

for (eval, modeset) in groupedeigenvalues
    if round(eval,digits=10)==0
        mlist = [m for m in modeset]
        println("length of mlist ",length(mlist))
    end
    # println(eval)
end

# [eval for (eval, modeset) in groupedeigenvalues if round(eval,digits=10)==0]

groundstateprojector = groundstates|>crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector
H = CrystalFockMap(energyspectrum) |> CrystalDenseMap

crystalfock = correlations|>getoutspace
crystal = crystalfock|>getcrystal

scaling = 2

blockedcorrelations,blockedH,blocker = blocking(correlations,H,scaling)

blockedcrystalfock = blockedcorrelations|>getoutspace
blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal
blockedspace::RealSpace = blockedcrystal|>getspace
    
refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

noofflavourpermode=1

firstcenter = [0,0] ∈ blockedspace
firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=firstcenter, metricspace=blockedspace)
firsthexagonalregionfock = quantize(firsthexagonalregion,noofflavourpermode)

firstinnerhexagonalregion = rankedandgroupoffsets(firsthexagonalregion,1)
offstlist = [offset for offset in firstinnerhexagonalregion]
firstinnerhexagonalregionhalf1 = Subset([offstlist[1],c3*offstlist[1],(c3)^2*offstlist[1]])
firstinnerhexagonalregionhalf2 = Subset([c6*offstlist[1],c6*c3*offstlist[1],c6*(c3)^2*offstlist[1]])
firstinnerhexagonalregionfock = quantize(firstinnerhexagonalregion,1)
firstinnerhexagonalregionhalf1fock = quantize(firstinnerhexagonalregionhalf1,1)
firstinnerhexagonalregionhalf2fock = quantize(firstinnerhexagonalregionhalf2,1)

firstboundaryhexagonalregion = firsthexagonalregion-firstinnerhexagonalregion
firstboundaryhexagonalregionfock = quantize(firstboundaryhexagonalregion,1)

# @info "Computing local correlations..."
localrestrict = fourier(blockedcrystalfock, firsthexagonalregionfock) / (blockedcrystal|>vol|>sqrt)
localcorrelations = localrestrict'*blockedcorrelations*localrestrict
localspectrum = localcorrelations|>eigspec
display(localspectrum|>visualize)
for val in [(localcorrelations[mode,mode]|>rep)[1,1]|>abs for mode in firsthexagonalregionfock]
    println(val)
end

ref = [1,2,3,4]
ref[2:2:4]