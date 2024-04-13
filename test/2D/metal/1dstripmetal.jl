using LinearAlgebra
using Zipper, Plots
plotlyjs()

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
fiodir("/Users/alphaharrius/ZERData/squaremetal/RG1")
metalcorrelations = fioload("stripmetalcorrelations")
metalspectrum = metalcorrelations|>crystalspectrum
metalspectrum|>visualize

metalcrystal = metalspectrum|>getcrystal
metalspace = metalcrystal|>getspace
normalvector = metalspace*(-1, 1)
contractbasis = metalspace*(0, getboundsize(metalcrystal)|>first)
@info "Contract basis: $contractbasis"
scaledspace = affinespace(normalvector, contractbasis)

extendedscale = Scale(scaledspace|>rep, metalspace)
extendedrestrict = ExtendedRestrict(extendedscale, normalvector, 1.5)
restrict = extendedrestrict * (metalcorrelations|>getoutspace)
stripcorrelations = restrict * metalcorrelations * restrict'
stripspectrum = stripcorrelations|>crystalspectrum
stripspectrum|>linespectrum|>visualize

bands = groupbands(stripspectrum, bandgrouping=[2, getbandcount(stripspectrum)-4, 2])

filledbands = bands[1]

function _getregionfock(crystalfock::CrystalFock, region::Region)
    crystalspace = crystalfock|>getcrystal|>getspace
    transformed = (crystalspace*r for r in region)
    offsets = (r-basispoint(r) for r in transformed)
    modes = (m+r for (m, r) in Iterators.product(crystalfock|>unitcellfock|>RegionFock, offsets))

    return transformed

    physregion = Subset(r|>euclidean for r in region)
    return RegionFock(m for m in modes if m|>getpos|>euclidean ∈ physregion)
end

region = getcrosssection(crystal=metalspectrum|>getcrystal, normalvector=normalvector*1.5, radius=0.5, minbottomheight=0.1)
region = region .- metalspace*(1, 1)*0.5
sspace = stripspectrum|>getcrystal|>getspace
[b-basispoint(b) for b in regionfock][1]|>getspace|>rep
visualize(stripspectrum|>getcrystal|>getunitcell, region)
regionfock = _getregionfock(stripcorrelations|>getoutspace, region)
transform = fourier(stripcorrelations|>getoutspace, regionfock)
localcorrelations = transform'*stripcorrelations*transform/(stripspectrum|>getcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[2, 17-2-2, 2])

localfilled = localstates[1]
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
rightrestrict = (
    fourier(wannierfilledisometry|>getinspace, wannierfilledisometry|>getinspace|>unitcellfock|>RegionFock) / (wanniercrystal|>vol|>sqrt))
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

localempty = localstates[3]
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
rightrestrict = (
    fourier(wannieremptyisometry|>getinspace, wannieremptyisometry|>getinspace|>unitcellfock|>RegionFock) / (wanniercrystal|>vol|>sqrt))
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

globaldistiller = wannieremptyprojector - wannierfilledprojector
distillerspectrum = globaldistiller|>crystalspectrum
distillerspectrum|>visualize

distilledbands = groupbands(distillerspectrum, :courier=>(v -> v < 1e-5 && v > -1e-5))

courierband = distilledbands[:courier]
courierprojector = courierband|>crystalprojector
couriercorrelations = 1 - courierprojector

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
