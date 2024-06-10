using Zipper
using LinearAlgebra

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)
[sqrt(3)/2 -1/2; 0. 1.]'*[2/(3*sqrt(3)), 1/(3*sqrt(3))]
[sqrt(3)/2 -1/2; 0. 1.]'*[1/(3*sqrt(3)), 4/(6*sqrt(3))]
[1 0; 1/2 sqrt(3)/2]'*[1/3, 0]
[1 0; 1/2 sqrt(3)/2]'*[0, 1/3]

pa = [2/(3*sqrt(3)), 1/(3*sqrt(3))] ∈ triangular
pb = [1/(3*sqrt(3)), 2/(3*sqrt(3))] ∈ triangular
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

pa|>euclidean
pb|>euclidean

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

unitcell = Subset(pa, pb)
visualize(unitcell)

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

    pa|>euclidean

    unitcell = Subset(pa, pb, pc, pd, pe, pf)


crystal = Crystal(unitcell, [24, 24])
# collect(crystal|>getspace|>getbasisvectors)
# new1 = [2/3, -1/3] ∈ triangular
# new2 = [-1/3, 2/3] ∈ triangular
# new3 = [0, 0] ∈ triangular
# new4 = [1/3, 1/3] ∈ triangular
# parallelogramref = Subset(new1,new2,new3,new4)
# new1.pos-new1.pos
# new1.pos
# [sqrt(3)/2 -1/2; 0. 1.]'*new1.pos
c3
inv([2/3 -1/3; -1/3 2/3]')

res = gethexagonalregion(rot = inv([2/3 -1/3; -1/3 2/3]'),crystal=crystal, center=new3, metricspace=crystal|>getspace)
visualize(res)

reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1 = members(modes)

t_a = -0.3 + 0im
t_b = -0.7 + 0im
tₙ = -1 + 0im

onsite = [
    (m1, m1) => t_b,
    (m0, m0) => t_a
]

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:r => [-1, 0] ∈ triangular)) => tₙ,
    (m0, m1 |> setattr(:r => [0, 1] ∈ triangular)) => tₙ]

# bonds::FockMap = bondmap([onsite..., nearestneighbor...])
bonds::FockMap = bondmap([nearestneighbor...])

energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstateprojector = groundstates|>crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector

H = crystalfockmap(energyspectrum)

function groupdictwifvalue(dict)
    sortedtupledata = sort([(value,key) for (key,value) in dict],by=first)
    refvalue = round(sortedtupledata[1][1],digits=10)
    result = []
    subresult = Subset(sortedtupledata[1][2])+Subset(sortedtupledata[2][2])
    for pair in sortedtupledata
        if round(pair[1],digits=10)==refvalue
            subresult = subresult+Subset(pair[2])
        else
            append!(result,[refvalue,subresult])
            refvalue = round(pair[1],digits=10)
            subresult = Subset(pair[2])
        end
    end
    return result
end

# [norm(euclidean(triangular*pt)) for pt in getsphericalregionwifrotsymm(crystal=crystal, radius=1, metricspace=triangular)]

crystalfock = correlations|>getoutspace

rgscale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...")
@info("Generating rgblocking transformation...")
rgblocker = @time rgscale * crystalfock
@info("Performing rgblocking on correlations...")
rgblockedcorrelations = @time rgblocker * correlations * rgblocker'
rgblockedcrystalfock = rgblockedcorrelations|>getoutspace
rgblockedcrystal::Crystal = rgblockedcrystalfock|>getcrystal
rgblockedspace::RealSpace = rgblockedcrystal|>getspace

hexagonalregion = (getsphericalregion(crystal=rgblockedcrystal, radius=0.9, metricspace=rgblockedspace|>orthospace))
visualize(hexagonalregion)
hexagonalregionfock = quantize(hexagonalregion,1)
localcorrelations = regioncorrelations(rgblockedcorrelations,quantize(hexagonalregion,1))
localspectrum = localcorrelations|>eigspech
localspectrum|>visualize
ref = groupdictwifvalue(localspectrum|>geteigenvalues)
for data in ref[1:20]
    println(data)
end

localiso = localcorrelations|>modeselectionbycount(3)
localwannierseedslists(localiso)






# @info "Computing frozen restricter..."
# distillregion::Region = getsphericalregion(crystal=rgblockedcrystal, radius=1, metricspace=rgblockedspace|>orthospace)
# visualize(distillregion)
# frozenseedingfock::RegionFock = quantize(distillregion, 1)
#     frozenrestrict = fourier(blockedcrystalfock, frozenseedingfock)

ratio = 2
rgscale2 = Scale([ratio 0; 0 ratio], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...")
@info("Generating rgblocking transformation...")
rgblocker2 = @time rgscale2 * rgblockedcrystalfock
@info("Performing rgblocking on correlations...")
rgblockedcorrelations2 = @time rgblocker2 * rgblockedcorrelations * rgblocker2'
rgblockedcrystalfock2 = rgblockedcorrelations2|>getoutspace
rgblockedcrystal2::Crystal = rgblockedcrystalfock2|>getcrystal
rgblockedspace2::RealSpace = rgblockedcrystal2|>getspace
rgblockedcrystal2|>getunitcell
distillregion2::Region = getsphericalregion(crystal=rgblockedcrystal2, radius=1, metricspace=rgblockedspace2|>orthospace)
visualize(distillregion2)
108/6
(rgblockedspace2|>orthospace).rep
(rgblockedspace2).rep

visualize(distillregion2)
[pt for pt in rgblockedcrystal2|>getunitcell  if (((pt.pos)[1]<1/2) & ((pt.pos)[2]<1/2))]

center1 = [0,0] ∈ rgblockedspace2
hexagonalregion1 = gethexagonalregion(crystal=rgblockedcrystal2,center = center1,metricspace=rgblockedspace2,ratio=ratio)
visualize(hexagonalregion1)

center2 = [1,1/2] ∈ rgblockedspace2
hexagonalregion2 = gethexagonalregion(crystal=rgblockedcrystal2,center = center2,metricspace=rgblockedspace2,ratio=ratio)
visualize(hexagonalregion2)

center3 = [1/2,1] ∈ rgblockedspace2
hexagonalregion3 = gethexagonalregion(crystal=rgblockedcrystal2,center = center3,metricspace=rgblockedspace2,ratio=ratio)
visualize(hexagonalregion3)

wholeregion = (hexagonalregion1+hexagonalregion2+hexagonalregion3)
rgblockedcrystal2|>getunitcell
intersect(wholeregion,rgblockedcrystal2|>getunitcell)==rgblockedcrystal2|>getunitcell

visualize(regioncorrelations(rgblockedcorrelations2,quantize(hexagonalregion,1))|>eigspech)
visualize(regioncorrelations(rgblockedcorrelations2,quantize(shifthexagonalregion1,1))|>eigspech)
visualize(regioncorrelations(rgblockedcorrelations2,quantize(shifthexagonalregion2,1))|>eigspech)

localiso = regioncorrelations(rgblockedcorrelations2,quantize(hexagonalregion,1))|>modeselectionbycount(18)

localwannierseedslists(localiso)