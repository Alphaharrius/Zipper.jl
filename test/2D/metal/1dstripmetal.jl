using LinearAlgebra
using Zipper, Plots
plotlyjs()

usecrystaldensemap()
setmaxthreads(Threads.nthreads())

square = euclidean(RealSpace, 2)
point = square*(1/2, 1/2)
spatialsnappingcalibration([point])

kspace = convert(MomentumSpace, square)

c4 = pointgrouptransform([0 -1; 1 0], localspace=square)
c2 = pointgrouptransform([-1 0; 0 -1], localspace=square)
m135 = pointgrouptransform([0 -1; -1 0], localspace=square)
m45 = pointgrouptransform([0 1; 1 0], localspace=square)

unitcell = Subset(point)
bc = BoundaryCondition(square, 128, 128)
crystal = Crystal(unitcell, bc)
reciprocalhashcalibration(bc.bounds)

@info "Loading data..."
fiodir("/Users/alphaharrius/ZERData/squaremetal/RG2")
metalcorrelations = fioload("stripmetalcorrelations")
metalspectrum = metalcorrelations|>crystalspectrum
metalspectrum|>visualize
metalcrystal = metalspectrum|>getcrystal
metalspace = metalcrystal|>getspace

rotation = Scale([1 -1; 1 1], metalspace)
rotate = rotation*(metalcorrelations|>getoutspace)
rotatedcorrelations = rotate * metalcorrelations * rotate'
rotatedcorrelations|>crystalspectrum|>visualize

rotatedcorrelations = couriercorrelations

blocking = Scale([2 0; 0 1], rotatedcorrelations|>getoutspace|>getcrystal|>getspace)
block = blocking*(rotatedcorrelations|>getoutspace)
blockedcorrelations = block * rotatedcorrelations * block'
blockedcorrelations|>crystalspectrum|>visualize

blockedcrystal = block|>getoutspace|>getcrystal

blockedcrystal|>getunitcell|>visualize

blockedcrystal|>brillouinzone|>visualize

region = blockedcrystal|>getunitcell
region = region + (region.+metalspace*(-1, 1)/2)
visualize(region, blockedcrystal|>getunitcell)
regionfock = getregionfock(blockedcorrelations|>getoutspace, region)
transform = fourier(blockedcorrelations|>getoutspace, regionfock)
localcorrelations = transform'*blockedcorrelations*transform/(blockedcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[4, 24-4-4, 4])
filledstates = localstates[1]
emptystates = localstates[3]

filledseeds = filledstates|>FockMap
filledisometry = transform * filledseeds
filledprojector = filledisometry .* filledisometry'
filledprojector|>crystalspectrum|>visualize

emptyseeds = emptystates|>FockMap
emptyisometry = transform * emptyseeds
emptyprojector = emptyisometry .* emptyisometry'
emptyprojector|>crystalspectrum|>visualize

globaldistiller = emptyprojector - filledprojector
distillerspectrum = globaldistiller|>crystalspectrum
distillerspectrum|>visualize

bands = groupbands(globaldistiller|>crystalspectrum, :empty=>(v -> v > 1e-5), :filled=>(v -> v < -1e-5))
courierbands = bands[:others]
courierprojector = courierbands|>crystalprojector
testcorrelations = courierprojector * blockedcorrelations * courierprojector'
testcorrelations|>crystalspectrum|>visualize
testspectrum = roundingpurification(testcorrelations|>crystalspectrum, tolerance=0.4)
testspectrum|>visualize
testcorrelations = testspectrum|>crystalfockmap

couriercorrelations = 1 - courierprojector

blockedspace = blockedcrystal|>getspace
blockedcrystal|>getunitcell|>visualize
region = Subset(blockedspace*(0.125, 0), blockedspace*(0.375, 0))
regions = [region, region.+blockedspace*(-0.03125, 0.25)]
regions = [regions..., (g.+blockedspace*(0, 0.5) for g in regions)...]
regions = [regions..., (g.+blockedspace*(0.5, 0) for g in regions)...]
visualize(regions...)
regionfocks = [getregionfock(couriercorrelations|>getoutspace, g) for g in regions]
transforms = [fourier(couriercorrelations|>getoutspace, fock) for fock in regionfocks]
localcorrelations = [t'*couriercorrelations*t for t in transforms]
localstates = sum(getregionstates(localcorrelations=corr, grouping=[1])[1] for corr in localcorrelations)
visualize(localstates|>normalize, markersize=5, logscale=0.9)

localseeds = localstates|>FockMap
transform = fourier(couriercorrelations|>getoutspace, localseeds|>getoutspace)
crystalseeds = transform * localseeds
crystalseeds = Dict(k=>crystalseeds[k, :] for k in couriercorrelations|>getoutspace|>getcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing courier states..."
courierisometry = wannierprojection(
    crystalisometries=courierbands|>geteigenvectors, crystal=couriercorrelations|>getoutspace|>getcrystal, crystalseeds=crystalseeds)

righttransform = fourier(courierisometry|>getinspace, courierisometry|>getinspace|>unitcellfock|>RegionFock)
courierlocalstates = transform'*courierisometry*righttransform|>RegionState
visualize(courierlocalstates|>normalize, markersize=5, logscale=0.9)

couriercorrelations = courierisometry'*blockedcorrelations*courierisometry
couriercorrelations = roundingpurification(couriercorrelations|>crystalspectrum, tolerance=0.4)|>crystalfockmap
couriercorrelations|>crystalspectrum|>visualize
couriercorrelations|>getoutspace|>getcrystal|>getunitcell|>visualize

# Strip distillation approach
metalcrystal = metalspectrum|>getcrystal
metalspace = metalcrystal|>getspace
normalvector = metalspace*(1, 1)
contractbasis = metalspace*(0, getboundsize(metalcrystal)|>first)
@info "Contract basis: $contractbasis"
scaledspace = affinespace(normalvector, contractbasis)

extendedscale = Scale(scaledspace|>rep, metalspace)
extendedrestrict = ExtendedRestrict(extendedscale, normalvector, 1.5)
restrict = extendedrestrict * (metalcorrelations|>getoutspace)
stripcorrelations = restrict * metalcorrelations * restrict'
stripspectrum = stripcorrelations|>crystalspectrum
stripspectrum|>linespectrum|>visualize

frozenthreshold = 2e-4
stripfrozenbands = groupbands(stripspectrum, :frozen=>(v -> v < frozenthreshold || v > 1-frozenthreshold))[:frozen]
stripfrozenbands
stripfrozenprojector = stripfrozenbands|>crystalprojector
stripfrozencorrelations = 1 - stripfrozenprojector

stripunitcell = stripspectrum|>getcrystal|>getunitcell
truncateregion = sum(stripunitcell.+normalvector*n for n in 0:5)
visualize(truncateregion, stripunitcell)
truncatedcorrelations = truncatetoregion(stripfrozencorrelations, truncateregion)
truncatedspectrum = truncatedcorrelations|>crystalspectrum
truncatedspectrum|>linespectrum|>visualize

metallicbands = groupbands(truncatedspectrum, bandgrouping=[2])[1]
metallicprojector = metallicbands|>crystalprojector
metalliccorrelations = 1 - metallicprojector

region = stripunitcell
region|>visualize
regionfock = getregionfock(metalliccorrelations|>getoutspace, region)
transform = fourier(metalliccorrelations|>getoutspace, regionfock)
localcorrelations = transform'*metalliccorrelations*transform/(stripspectrum|>getcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates1 = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]
localstates = m45*localstates1
visualize(localstates|>normalize, markersize=5, logscale=0.5)

region = stripunitcell .+ normalvector
regionfock = getregionfock(metalliccorrelations|>getoutspace, region)
transform = fourier(metalliccorrelations|>getoutspace, regionfock)
localcorrelations = transform'*metalliccorrelations*transform/(stripspectrum|>getcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates2 = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]
localstates += m45*localstates
visualize(localstates|>normalize, markersize=5, logscale=0.5)

localseeds = localstates|>FockMap

normalvector = metalspace*(1, 1)*2
contractbasis = metalspace*(0, getboundsize(metalcrystal)|>first)
@info "Contract basis: $contractbasis"
scaledspace = affinespace(normalvector, contractbasis)

extendedscale = Scale(scaledspace|>rep, metalspace)
extendedrestrict = ExtendedRestrict(extendedscale, normalvector, 1.5)
restrict = extendedrestrict * (metalcorrelations|>getoutspace)
translatecrystalfock = restrict|>getoutspace

remapper = spatialremapper(localseeds|>getoutspace, translatecrystalfock)
localseeds = remapper * localseeds
transform = fourier(translatecrystalfock, localseeds|>getoutspace)
crystalseeds = transform * localseeds
crystalseeds = Dict(k=>crystalseeds[k, :] for k in translatecrystalfock|>getcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing metal bands..."
metallicisometry = wannierprojection(
    crystalisometries=metallicbands|>geteigenvectors, crystal=metalliccorrelations|>getoutspace|>getcrystal, crystalseeds=crystalseeds)
rightfock = metallicisometry|>getinspace|>unitcellfock|>RegionFock
rightrestrict = fourier(metallicisometry|>getinspace, rightfock)
metallocalstates = transform'*metallicisometry*rightrestrict|>RegionState|>normalize
visualize(metallocalstates|>normalize, markersize=10, logscale=0.9)

scale = Scale([1 -1; 1 1], metalspace)
rotate = scale*(metalcorrelations|>getoutspace)
rotate|>getoutspace|>getcrystal|>getunitcell|>visualize
remapper = spatialremapper(metallocalstates|>getoutspace, rotate|>getoutspace)
metallocalisometry = remapper*FockMap(metallocalstates)
transform = fourier(rotate|>getoutspace, metallocalisometry|>getoutspace)
crystalisometry = transform * metallocalisometry
stripmetalprojector = crystalisometry .* crystalisometry'
stripmetalspectrum = stripmetalprojector|>crystalspectrum


rotatedcorrelations = rotate * metalcorrelations * rotate'
stripmetalprojector*rotatedcorrelations*stripmetalprojector'|>crystalspectrum|>visualize

projectedspectrum = metallicisometry'*stripcorrelations*metallicisometry|>crystalspectrum
projectedspectrum|>linespectrum|>visualize

visualize(getsphericalregion(crystal=stripspectrum|>getcrystal, radius=16, metricspace=scaledspace|>euclidean), stripspectrum|>getcrystal|>getunitcell)

bands = groupbands(stripspectrum, bandgrouping=[2, getbandcount(stripspectrum)-4, 2])

filledbands = bands[1]

region1 = getcrosssection(crystal=metalspectrum|>getcrystal, normalvector=inv(c4)*normalvector*1.75, radius=1.25, minbottomheight=0.25)
region1 = region1 .- inv(c4)*normalvector
visualize(region1, stripspectrum|>getcrystal|>getunitcell)

region2|>getcenter|>latticeoff

regionfock = getregionfock(stripcorrelations|>getoutspace, region1)
transform = fourier(stripcorrelations|>getoutspace, regionfock)
localcorrelations = transform'*stripcorrelations*transform/(stripspectrum|>getcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates1 = getregionstates(localcorrelations=localcorrelations, grouping=[2, 20, 2])
regionfock = getregionfock(stripcorrelations|>getoutspace, region2)
transform = fourier(stripcorrelations|>getoutspace, regionfock)
localcorrelations = transform'*stripcorrelations*transform/(stripspectrum|>getcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates2 = getregionstates(localcorrelations=localcorrelations, grouping=[2, 6])

localfilled = localstates1[1]
visualize(localfilled|>normalize, markersize=5, logscale=0.5)

filledseeds = localfilled|>FockMap
visualize(filledseeds|>RegionState, markersize=5, logscale=0.5)
transform = fourier(stripcorrelations|>getoutspace, filledseeds|>getoutspace)/(stripspectrum|>getcrystal|>vol)
filledcrystalseeds = transform * filledseeds
filledcrystalseeds = Dict(k=>filledcrystalseeds[k, :] for k in stripspectrum|>getcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in filledcrystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing strip filled states..."
wannierfilledisometry = wannierprojection(
    crystalisometries=filledbands|>geteigenvectors, crystal=stripspectrum|>getcrystal, crystalseeds=filledcrystalseeds)
wanniercrystal = wannierfilledisometry|>getinspace|>getcrystal
rightrestrict = fourier(wannierfilledisometry|>getinspace, wannierfilledisometry|>getinspace|>unitcellfock|>RegionFock)
wannierlocalisometry = transform' * wannierfilledisometry * rightrestrict
stripfilledlocalstates = wannierlocalisometry|>RegionState|>normalize
visualize(stripfilledlocalstates|>normalize, markersize=5, logscale=0.9)

scale = Scale([1 -1; 1 1], metalspace)
scaledoutspace = scale*(metalcorrelations|>getoutspace)|>getoutspace

remapper = spatialremapper(wannierlocalisometry|>getoutspace, scaledoutspace)
filledisometry = remapper * wannierlocalisometry
visualize(filledisometry|>RegionState|>normalize, markersize=5, logscale=0.9)
transform = fourier(scaledoutspace, filledisometry|>getoutspace)
wannierfilledisometry = transform * filledisometry
wannierfilledprojector = wannierfilledisometry .* wannierfilledisometry'
wannierfilledprojector|>crystalspectrum|>visualize

emptybands = bands[3]

localempty = localstates1[3]
visualize(localempty|>normalize, markersize=5, logscale=0.5)

emptyseeds = localempty|>FockMap
visualize(emptyseeds|>RegionState, markersize=5, logscale=0.5)
transform = fourier(stripcorrelations|>getoutspace, emptyseeds|>getoutspace)/(stripspectrum|>getcrystal|>vol)
emptycrystalseeds = transform * emptyseeds
emptycrystalseeds = Dict(k=>emptycrystalseeds[k, :] for k in stripspectrum|>getcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in emptycrystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing strip empty states..."
wannieremptyisometry = wannierprojection(
    crystalisometries=emptybands|>geteigenvectors, crystal=stripspectrum|>getcrystal, crystalseeds=emptycrystalseeds)
wanniercrystal = wannieremptyisometry|>getinspace|>getcrystal
rightrestrict = fourier(wannieremptyisometry|>getinspace, wannieremptyisometry|>getinspace|>unitcellfock|>RegionFock)
wannierlocalisometry = transform' * wannieremptyisometry * rightrestrict
stripemptylocalstates = wannierlocalisometry|>RegionState|>normalize
visualize(stripemptylocalstates|>normalize, markersize=5, logscale=0.9)

scale = Scale([1 -1; 1 1], metalspace)
scaledoutspace = scale*(metalcorrelations|>getoutspace)|>getoutspace

remapper = spatialremapper(wannierlocalisometry|>getoutspace, scaledoutspace)
emptyisometry = remapper * wannierlocalisometry
visualize(emptyisometry|>RegionState|>normalize, markersize=5, logscale=0.9)
transform = fourier(scaledoutspace, emptyisometry|>getoutspace)
wannieremptyisometry = transform * emptyisometry
wannieremptyprojector = wannieremptyisometry .* wannieremptyisometry'
wannieremptyprojector|>crystalspectrum|>visualize

globaldistiller = wannierfilledprojector
distillerspectrum = globaldistiller|>crystalspectrum
distillerspectrum|>visualize

distilledbands = groupbands(distillerspectrum, :courier=>(v -> v < 1e-5 && v > -1e-5))

courierband = distilledbands[:courier]
courierprojector = courierband|>crystalprojector
couriercorrelations = 1 - courierprojector

frozenprojector = distilledbands[:others]|>crystalprojector

rotate = scale*(metalcorrelations|>getoutspace)
rotatedcorrelations = rotate * metalcorrelations * rotate'
courierprojector*rotatedcorrelations*courierprojector'|>crystalspectrum|>visualize
rotatedcorrelations|>crystalspectrum|>visualize

couriercrystal = couriercorrelations|>getoutspace|>getcrystal
region = couriercrystal|>getunitcell
regionfock = getregionfock(couriercorrelations|>getoutspace, region)
transform = fourier(couriercorrelations|>getoutspace, regionfock)
localcorrelations = transform'*couriercorrelations*transform/(couriercrystal|>vol)
localcorrelations|>eigspech|>visualize
localstate = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]
visualize(localstate|>normalize, markersize=5, logscale=0.9)

blockcenter = block|>getoutspace|>getcrystal|>getunitcell|>getcenter
localstates = sum(g.*localstate for g in c4|>recenter(blockcenter)|>pointgroupelements)
visualize(localstates|>normalize, markersize=5, logscale=0.9)

localseeds = localstates|>FockMap
visualize(localseeds|>RegionState, markersize=5, logscale=0.9)
remapper = spatialremapper(localseeds|>getoutspace, block|>getoutspace)
localseeds = remapper * localseeds
transform = fourier(block|>getoutspace, localseeds|>getoutspace)/(block|>getoutspace|>getcrystal|>vol)
blockseeds = transform * localseeds
blockseeds = Dict(k=>blockseeds[k, :] for k in block|>getoutspace|>getcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in blockseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing courier states..."
wannierisometry = wannierprojection(
    crystalisometries=courierband|>geteigenvectors, crystal=couriercorrelations|>getoutspace|>getcrystal, crystalseeds=blockseeds)
