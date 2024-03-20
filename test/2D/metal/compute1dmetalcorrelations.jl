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

wannierregion = getsphericalregion(crystal=localcrystal, radius=0.25, metricspace=localcrystal|>getspace)
wannierregion = wannierregion .+ getspace(wannierregion)*(0.25, 0.25)
regionfock = quantize(wannierregion, 1)
restrict = fourier(localcrystalfock, regionfock) / (localcrystal|>vol|>sqrt)
localcorrelations = restrict' * frozencorrelations * restrict
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3, 1])[1]
localstates = m135 * localstates
localstates = m45 * localstates

localstates = sum(g.*localstates for g in pointgroupelements(c4|>recenter(localspace*(0.5, 0.5))))
visualize(localstates|>FockMap|>RegionState, markersize=5, logscale=0.5)
localseeds = localstates|>FockMap
transform = fourier(localcrystalfock, localseeds|>getoutspace) / (localcrystal|>vol|>sqrt)
crystalseeds = transform * localseeds
crystalseeds = Dict(k=>crystalseeds[k, :] for k in localcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh if v < 0.2)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

wannierfilledisometry = wannierprojection(
    crystalisometries=frozenbands|>geteigenvectors, crystal=localcrystal, crystalseeds=crystalseeds)

wannierfock = wannierfilledisometry|>getinspace|>unitcellfock|>RegionFock
visualregion = getsphericalregion(crystal=localcrystal, radius=2, metricspace=localspace) .+ localspace*(0.75, 0.75)
visualfock = quantize(visualregion, 1)
leftrestrict = fourier(wannierfilledisometry|>getoutspace, visualfock) / (localcrystal|>vol|>sqrt)
rightrestrict = fourier(wannierfilledisometry|>getinspace, wannierfock)
filledlocalstates = leftrestrict' * wannierfilledisometry * rightrestrict
visualize(filledlocalstates|>RegionState|>normalize, markersize=7, logscale=0.5)

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

fiodir("/Users/alphaharrius/ZERData/squaremetal/rg1/branch")
fiosave(metalcorrelations, name="metalcorrelations")
fiosave(wannierfilledisometry, name="metalisometry")
