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

unitcell = Subset(point)
crystal = Crystal(unitcell, [128, 128])
reciprocalhashcalibration(crystal.sizes)

fiodir("/Users/alphaharrius/ZERData/squaremetal/RG1")
globaldistiller = fioload("rg1globaldistiller")
globaldistillerspectrum = globaldistiller|>crystalspectrum
globaldistillerspectrum|>visualize

frozenbands = distillation(globaldistillerspectrum, :frozen=>(v -> v > 0.5))[:frozen]
frozenbands|>visualize
frozenprojector = frozenbands|>crystalprojector
frozencorrelations = idmap(frozenprojector|>getoutspace) - frozenprojector

crystalfock = globaldistiller|>getoutspace
localcrystal = globaldistillerspectrum|>getcrystal
localspace = localcrystal|>getspace
normalvector = (1, 1) ∈ localspace

extendedrestrict = extendedcrystalrestrict(
        crystal=localcrystal, normalvector=normalvector, stripradius=0.5)
restrict = extendedrestrict * crystalfock
stripcorrelations = restrict * frozencorrelations * restrict'
stripspectrum = stripcorrelations|>crystalspectrum
stripspectrum|>linespectrum|>visualize

stripfrozensstates = distillation(stripspectrum, :filled=>(v -> v < 0.5))[:filled]

fiodir("/Users/alphaharrius/ZERData/squaremetal")
physcorrelations = fioload("inputC")
filledprojector = idmap(physcorrelations|>getoutspace) - physcorrelations

stripfrozenprojector = stripfrozensstates|>crystalprojector
stripfrozenprojector = restrict' * stripfrozenprojector * restrict
fiodir("/Users/alphaharrius/ZERData/squaremetal/RG1")
block = fioload("rg1block")
filledprojector = block * filledprojector * block'
stripfilledprojector = filledprojector * stripfrozenprojector * filledprojector'

stripfilledprojector = restrict * stripfilledprojector * restrict'
stripfilledcorrelations = idmap(stripfilledprojector|>getoutspace) - stripfilledprojector
stripfilledspectrum = stripfilledcorrelations|>crystalspectrum
stripfilledspectrum|>linespectrum|>visualize

@info "Computing strip truncation restricted Fourier transform..."
stripunitcell = getcrosssection(crystal=localcrystal, normalvector=normalvector*3, radius=0.5)
striphomefock = stripfilledcorrelations|>getoutspace|>unitcellfock
basistohomemodes = ((mode|>getattr(:b)|>basispoint)=>mode for mode in striphomefock)
conversionmappings = Dict(
    mode=>(mode|>setattr(:r=>getattr(mode, :b)-b)|>setattr(:b=>b)) for (b, mode) in basistohomemodes)
actualstriphomefock = conversionmappings|>values|>RegionFock
truncationregionfock = RegionFock(
    mode|>setattr(:r=>getattr(mode, :r)+normalvector*n) for mode in actualstriphomefock for n in 0:2)
homemappings = Dict(tomode=>frommode|>removeattr(:r) for (tomode, frommode) in conversionmappings)
truncate = (
    fourier(stripfilledcorrelations|>getoutspace, truncationregionfock, homemappings) / sqrt(stripfilledcorrelations|>getoutspace|>getcrystal|>vol))

@info "Performing truncation..."
truncatedregioncorrelations = truncate' * stripfilledcorrelations * truncate
offsets = Subset(normalvector * n for n in 0:2)
remapper = spatialremapper(truncationregionfock; offsets=offsets, unitcell=stripunitcell)
truncatedregioncorrelations = remapper * truncatedregioncorrelations * remapper'

truncationregionindices = Iterators.product(
    truncatedregioncorrelations|>getoutspace, truncatedregioncorrelations|>getoutspace)
truncationregionbonds = Dict(
    (getattr(tomode, :r) - getattr(frommode, :r), frommode|>removeattr(:r), tomode|>removeattr(:r))=>(frommode, tomode) 
    for (frommode, tomode) in truncationregionindices)
pruningindices = (index for (_, index) in truncationregionbonds)
prunedcorrelations = extractindices(truncatedregioncorrelations, pruningindices)
prunedcorrelations = remapper' * prunedcorrelations * remapper
stripfouriers = (truncate[subspace, :] for (_, subspace) in truncate|>getoutspace|>crystalsubspaces)
stripcrystal = truncate|>getoutspace|>getcrystal
truncatedstripfilledcorrelations = crystaldirectsum(
    (transform * prunedcorrelations * transform' for transform in stripfouriers), 
    outcrystal=stripcrystal, incrystal=stripcrystal)
truncatedstripfilledcorrelationspectrum = truncatedstripfilledcorrelations|>crystalspectrum
truncatedstripfilledcorrelationspectrum|>linespectrum|>visualize

stripfilledstates = groundstatespectrum(truncatedstripfilledcorrelationspectrum, perunitcellfillings=6)
stripfilledstates|>linespectrum|>visualize

stripfilledprojector = stripfilledstates|>crystalprojector
stripfilledcorrelations = idmap(stripfilledprojector|>getoutspace) - stripfilledprojector

scaledcrystal = stripfilledstates|>getcrystal
scaledspace = scaledcrystal|>getspace

wannierregion1 = getcrosssection(crystal=localcrystal, normalvector=normalvector*0.875, radius=0.5, minbottomheight=0.15)
regionfock = quantize(wannierregion1, 1)
remapper = spatialremapper(regionfock, offsets=scaledspace|>getorigin|>Subset, unitcell=scaledcrystal|>getunitcell)
wannierregionfock = remapper|>getoutspace
restrict = fourier(stripfilledcorrelations|>getoutspace, wannierregionfock) / (scaledcrystal|>vol|>sqrt)
localcorrelations = restrict' * stripfilledcorrelations * restrict
localcorrelations = remapper' * localcorrelations * remapper
localcorrelations|>eigspech|>visualize

region1states = getregionstates(localcorrelations=localcorrelations, grouping=[3])[1]
region1states = m45 * region1states
region1states = m135 * region1states
region1states = remapper*FockMap(region1states)|>RegionState
visualize(region1states, markersizemultiplier=20, markersizescaling=0.3)

wannierregion2 = wannierregion1 .- normalvector*0.5
regionfock = quantize(wannierregion2, 1)
offsets = Subset(-(scaledspace*normalvector), scaledspace|>getorigin)
remapper = spatialremapper(regionfock, offsets=offsets, unitcell=scaledcrystal|>getunitcell)
wannierregionfock = remapper|>getoutspace
restrict = fourier(stripfilledcorrelations|>getoutspace, wannierregionfock) / (scaledcrystal|>vol|>sqrt)
localcorrelations = restrict' * stripfilledcorrelations * restrict
localcorrelations = remapper' * localcorrelations * remapper
localcorrelations|>eigspech|>visualize

region2states = getregionstates(localcorrelations=localcorrelations, grouping=[3])[1]
region2states = m45 * region2states
region2states = m135 * region2states
region2states = remapper*FockMap(region2states)|>RegionState
visualize(region2states, markersizemultiplier=20, markersizescaling=0.3)

localseeds = FockMap(region1states + region2states)
visualize(localseeds|>RegionState, markersizemultiplier=20, markersizescaling=0.3)

transform = fourier(stripfilledcorrelations|>getoutspace, localseeds|>getoutspace) / (scaledcrystal|>vol|>sqrt)
crystalseeds = transform * FockMap(localseeds)
crystalseeds = Dict(k=>crystalseeds[k, :] for k in scaledcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh if v < 0.2)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing strip filled states..."
wannierfilledisometry = wannierprojection(
    crystalisometries=stripfilledstates|>geteigenvectors, crystal=scaledcrystal, crystalseeds=crystalseeds)
wanniercrystal = wannierfilledisometry|>getinspace|>getcrystal
rightrestrict = (
    fourier(wannierfilledisometry|>getinspace, wannierfilledisometry|>getinspace|>unitcellfock|>RegionFock) / (wanniercrystal|>vol|>sqrt))
wannierlocalisometry = transform' * wannierfilledisometry * rightrestrict
stripfilledlocalstates = wannierlocalisometry|>RegionState
visualize(stripfilledlocalstates, markersizemultiplier=20, markersizescaling=0.3)

wannierlocalisometry|>getinspace|>unitcellfock|>showmodes

stripfilledprojector = wannierfilledisometry * wannierfilledisometry'

stripspectrum|>linespectrum|>visualize
distillation(stripspectrum, :frozen=>(v -> v < 0.5))[:frozen]
stripfrozenprojector = idmap(stripcorrelations|>getoutspace) - stripcorrelations
stripfrozenprojector|>crystalspectrum|>linespectrum|>visualize
trialspectrum = stripfilledprojector * stripfrozenprojector * stripfilledprojector'|>crystalspectrum
distillation(trialspectrum, :frozen=>(v -> v > 0.5))[:frozen]

blockedcorrelations = block * physcorrelations * block'

wannierfilledisometry' * wannierfilledisometry |>FockMap|>visualize

blockedcorrelations|>crystalspectrum|>visualize
distillation(blockedcorrelations|>crystalspectrum, :frozen=>(v -> v > 0.5))[:frozen]

blockedprojector = idmap(blockedcorrelations|>getoutspace) - blockedcorrelations
blockedfrozenprojector = frozenprojector * blockedprojector * frozenprojector'
blockedfrozencorrelations = idmap(blockedfrozenprojector|>getoutspace) - blockedfrozenprojector
distillation(blockedfrozencorrelations|>crystalspectrum, :frozen=>(v -> v > 0.5))[:frozen]
purifiedfrozencorrelations = blockedfrozencorrelations|>crystalspectrum|>roundingpurification|>CrystalFockMap
purifiedfrozencorrelations|>crystalspectrum|>visualize

wannierfilledprojector = wannierfilledisometry * wannierfilledisometry'

restrict = extendedrestrict * crystalfock
globalfilledisometry = wannierfilledprojector * restrict

purifiedfrozenprojector = idmap(purifiedfrozencorrelations|>getoutspace) - purifiedfrozencorrelations

projectedstripfilledprojector = globalfilledisometry * purifiedfrozenprojector * globalfilledisometry'
projectedstripfilledprojector |>crystalspectrum|>linespectrum|>visualize

projectedstripfilledcorrelations = idmap(projectedstripfilledprojector|>getoutspace) - projectedstripfilledprojector

projectedstripfilledcorrelations|>crystalspectrum|>linespectrum|>visualize
distillation(projectedstripfilledcorrelations|>crystalspectrum, :frozen=>(v -> v < 0.5))[:frozen]
