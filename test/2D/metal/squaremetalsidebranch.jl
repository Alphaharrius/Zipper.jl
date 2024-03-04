using LinearAlgebra
using Zipper, Plots

plotlyjs()

setmaxthreads(Threads.nthreads())

square = euclidean(RealSpace, 2)
point = [1/2, 1/2] ∈ square
spatialsnappingcalibration([point])

kspace = convert(MomentumSpace, square)

c4 = pointgrouptransform([0 -1; 1 0], localspace=square)
c2 = pointgrouptransform([-1 0; 0 -1], localspace=square)
m135 = pointgrouptransform([0 -1; -1 0], localspace=square)
m45 = pointgrouptransform([0 1; 1 0], localspace=square)

xxminusyy = BasisFunction(:xx=>1, :yy=>-1, dimension=2)
seteigenfunction(m135, xxminusyy)
seteigenfunction(m45, xxminusyy)

unitcell = Subset(point)
crystal = Crystal(unitcell, [128, 128])
reciprocalhashcalibration(crystal.sizes)

fiodir("/Users/alphaharrius/ZERData/squaremetal/rg1")

@info "Loading the global distiller..."
globaldistiller = fioload("globaldistiller")
globaldistiller|>crystalspectrum|>visualize

@info "Distilling frozen bands..."
frozenbands = groupbands(globaldistiller|>crystalspectrum, bandgrouping=[12], frombottom=false)[1]
frozenprojector = frozenbands|>crystalprojector
frozencorrelations = 1 - frozenprojector

localcrystal = frozenbands|>getcrystal
localspace = localcrystal|>getspace
localcrystalfock = frozencorrelations|>getoutspace

normalvector = (1, 1) ∈ localspace
extendedrestrict = extendedcrystalrestrict(
    crystal=localcrystal, normalvector=normalvector, stripradius=0.5)
restrict = extendedrestrict * localcrystalfock
stripcorrelations = restrict * frozencorrelations * restrict'
stripspectrum = stripcorrelations|>crystalspectrum
stripspectrum|>linespectrum|>visualize

scaledcrystal = stripspectrum|>getcrystal
scaledspace = scaledcrystal|>getspace

wannierregion = getcrosssection(crystal=localcrystal, normalvector=normalvector*0.874, radius=0.5, minbottomheight=0.25) .- normalvector*0.5
wannierregion|>visualize
regionfock = quantize(wannierregion, 1)
restrict = fourier(localcrystalfock, regionfock) / (scaledcrystal|>vol|>sqrt)
localcorrelations = restrict' * frozencorrelations * restrict
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[8])[1]
localstates = m45 * localstates
localstates = m135 * localstates
localstates = localstates[3] + localstates[4] + localstates[6]
localstates|>normalize
visualize(localstates, markersize=5, logscale=0.5)

fullstates = []
for g in pointgroupelements(c4|>recenter([0.25, 0.25]∈localspace))
    push!(fullstates, g .* localstates)
end
fullstates = sum(fullstates)
visualize(fullstates|>FockMap|>RegionState, markersize=5, logscale=0.7)
fullseeds = fullstates|>FockMap

transform = fourier(localcrystalfock, fullseeds|>getoutspace) / (localcrystal|>vol|>sqrt)
crystalseeds = transform * fullseeds
crystalseeds = Dict(k=>crystalseeds[k, :] for k in localcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh if v < 0.2)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

wannierfilledisometry = wannierprojection(
    crystalisometries=frozenbands|>geteigenvectors, crystal=localcrystal, crystalseeds=crystalseeds)

wannierfock = wannierfilledisometry|>getinspace|>unitcellfock|>RegionFock
visualregion = getsphericalregion(crystal=localcrystal, radius=2, metricspace=localspace)
visualfock = quantize(visualregion, 1)
leftrestrict = fourier(wannierfilledisometry|>getoutspace, visualfock) / (localcrystal|>vol|>sqrt)
rightrestrict = fourier(wannierfilledisometry|>getinspace, wannierfock)
filledlocalstates = leftrestrict' * wannierfilledisometry * rightrestrict
visualize(filledlocalstates|>RegionState|>normalize, markersize=7, logscale=0.7)

wannierfilledisometry' * frozencorrelations * wannierfilledisometry |>crystalspectrum|>visualize
wannierfilledprojector = wannierfilledisometry * wannierfilledisometry'

projecteddistiller = wannierfilledprojector * globaldistiller * wannierfilledprojector
projectedspectrum = projecteddistiller|>crystalspectrum
projectedspectrum|>visualize

fiodir("/Users/alphaharrius/ZERData/squaremetal/rg1/branch")
fiosave(projecteddistiller, name="metalbanddistiller")
fiosave(wannierfilledisometry, name="wanniermetalisometry")

distillation(projectedspectrum, :projected=>(v->v>0.5))

fiodir("/Users/alphaharrius/ZERData/squaremetal")
inputcorrelations = fioload("inputcorrelations")
fiodir("/Users/alphaharrius/ZERData/squaremetal/rg1")
block = fioload("block")
blockedcorrelations = block * inputcorrelations * block'

metalcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
metalspectrum = metalcorrelations|>crystalspectrum
metalspectrum = metalspectrum|>roundingpurification
metalspectrum|>visualize
metalcorrelations = metalspectrum|>CrystalFockMap
fiosave(metalcorrelations, name="metalcorrelations")
fiosave(wannierfilledisometry, name="wanniermetalisometry")
wannierfilledisometry|>getinspace|>unitcellfock|>RegionFock|>getregion|>visualize

fiodir("/Users/alphaharrius/ZERData/squaremetal/rg1/branch")
metalcorrelations = fioload("metalcorrelations")
metalcorrelations|>crystalspectrum|>visualize
transform = c4 * getoutspace(metalcorrelations)
