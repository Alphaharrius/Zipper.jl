using LinearAlgebra
using Zipper, Plots

plotlyjs()

setmaxthreads(Threads.nthreads())

@info "Preparing environment..."
square = euclidean(RealSpace, 2)
point = [1/2, 1/2] ∈ square
spatialsnappingcalibration([point])

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

@info "Loading data..."
fiodir("/Users/alphaharrius/ZERData/squaremetal/rg1/branch")
correlations = fioload("metalcorrelations")

crystalfock = correlations|>getoutspace
crystal = crystalfock|>getcrystal
space = crystal|>getspace

distillfock = crystalfock|>unitcellfock|>RegionFock
transform = fourier(crystalfock, distillfock)
localcorrelations = transform'*correlations*transform/(crystal|>vol|>sqrt)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1, 10, 1])

localfilledisometry = localstates[1]|>FockMap
crystalfilledisometry = transform * localfilledisometry
filledprojector = crystalfilledisometry .* crystalfilledisometry'

localemptyisometry = localstates[3]|>FockMap
crystalemptyisometry = transform * localemptyisometry
emptyprojector = crystalemptyisometry .* crystalemptyisometry'

globaldistiller = emptyprojector - filledprojector

distillregion = getsphericalregion(crystal=crystal, radius=0.5, metricspace=space)
distillfock = getregionfock(crystalfock, distillregion)
transform = fourier(crystalfock, distillfock)
localcorrelations = transform'*correlations*transform/(crystal|>vol|>sqrt)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1, 10, 1])

localfilledisometry = localstates[1]|>FockMap
crystalfilledisometry = transform * localfilledisometry
filledprojector = crystalfilledisometry .* crystalfilledisometry'

localemptyisometry = localstates[3]|>FockMap
crystalemptyisometry = transform * localemptyisometry
emptyprojector = crystalemptyisometry .* crystalemptyisometry'

globaldistiller += (emptyprojector - filledprojector)

distillerspectrum = globaldistiller|>crystalspectrum
distillerspectrum|>visualize

distilled = groupbands(distillerspectrum, :courier=>(v -> v > -1e-5 && v < 1e-5))

frozenbands = distilled[:others]
frozenprojector = frozenbands|>crystalprojector
frozencorrelations = 1 - frozenprojector

regionfock = getregionfock(frozencorrelations|>getoutspace, crystal|>getunitcell)
transform = fourier(frozencorrelations|>getoutspace, regionfock)
localcorrelations = transform'*frozencorrelations*transform/(crystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[2])[1]

regionfock = getregionfock(frozencorrelations|>getoutspace, getsphericalregion(crystal=crystal, radius=0.5, metricspace=space))
transform = fourier(frozencorrelations|>getoutspace, regionfock)
localcorrelations = transform'*frozencorrelations*transform/(crystal|>vol)
localcorrelations|>eigspech|>visualize
localstates += getregionstates(localcorrelations=localcorrelations, grouping=[2])[1]

localseeds = localstates|>FockMap
visualize(localseeds|>RegionState, markersize=5, logscale=0.5)
transform = fourier(frozencorrelations|>getoutspace, localseeds|>getoutspace)
crystalseeds = transform * localseeds
crystalseeds = Dict(k=>crystalseeds[k, :] for k in crystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

wannierfrozenisometry = wannierprojection(
    crystalisometries=frozenbands|>geteigenvectors, crystal=crystal, crystalseeds=crystalseeds)

visualregion = getsphericalregion(crystal=crystal, radius=2, metricspace=space)
visualfock = quantize(visualregion, 3)
lefttransform = fourier(wannierfrozenisometry|>getoutspace, visualfock) / (crystal|>vol|>sqrt)
righttransform = fourier(wannierfrozenisometry|>getinspace, wannierfrozenisometry|>getinspace|>unitcellfock|>RegionFock)
frozenlocalstates = lefttransform' * wannierfrozenisometry * righttransform
visualize(frozenlocalstates|>RegionState, markersize=1, logscale=0.7)

courierbands = distilled[:courier]
courierprojector = courierbands|>crystalprojector
couriercorrelations = 1 - courierprojector
regionfock = getregionfock(couriercorrelations|>getoutspace, space*(0.25, 0.25)|>Subset)
transform = fourier(couriercorrelations|>getoutspace, regionfock)
localcorrelations = transform'*couriercorrelations*transform/(crystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1, 1, 1])
localstates = localstates[1] + localstates[3]
localstates = sum(g.*localstates for g in c4|>recenter(space*(0.5, 0.5))|>pointgroupelements)
visualize(localstates, markersize=5, logscale=0.5)

localseeds = localstates|>FockMap
visualize(localseeds|>RegionState, markersize=5, logscale=0.5)
transform = fourier(couriercorrelations|>getoutspace, localseeds|>getoutspace)
crystalseeds = transform * localseeds
crystalseeds = Dict(k=>crystalseeds[k, :] for k in crystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

wanniercourierisometry = wannierprojection(
    crystalisometries=courierbands|>geteigenvectors, crystal=crystal, crystalseeds=crystalseeds)

visualregion = getsphericalregion(crystal=crystal, radius=1.5, metricspace=space) .+ space*(0.5, 0.5)
visualfock = quantize(visualregion, 2)
lefttransform = fourier(wanniercourierisometry|>getoutspace, visualfock) / (crystal|>vol|>sqrt)
righttransform = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
courierlocalstates = lefttransform' * wanniercourierisometry * righttransform
visualize(courierlocalstates|>RegionState|>normalize, markersize=5, logscale=0.5)

fiodir("/Users/alphaharrius/ZERData/squaremetal/rg1/branch")
metalcorrelations = fioload("metalcorrelations")

frozencorrelations = wannierfrozenisometry' * metalcorrelations * wannierfrozenisometry
frozencorrelations|>crystalspectrum|>visualize

couriercorrelations = wanniercourierisometry' * metalcorrelations * wanniercourierisometry
couriercorrelations|>crystalspectrum
couriercorrelations|>crystalspectrum|>visualize

@info "Performing blocking..."
couriercrystal = couriercorrelations|>getoutspace|>getcrystal
courierspace = couriercrystal|>getspace

scale = Scale([2 0; 0 2], courierspace)
block = scale * (couriercorrelations|>getoutspace)
blockedcorrelations = block * couriercorrelations * block'
blockedcrystalfock = blockedcorrelations|>getoutspace
blockedcrystal = blockedcrystalfock|>getcrystal

@info "Performing extended restrictions..."
normalvector = (1, 1) ∈ (blockedcrystal|>getspace)
extendedrestrict = extendedcrystalrestrict(
    crystal=blockedcrystal, normalvector=normalvector, stripradius=0.5)

restrict = extendedrestrict * blockedcrystalfock
stripcorrelations = restrict * blockedcorrelations * restrict'
stripoutspace = stripcorrelations|>getoutspace
stripcrystal = stripoutspace|>getcrystal
stripspace = stripcrystal|>getspace
stripspectrum = stripcorrelations|>crystalspectrum
stripspectrum|>linespectrum|>visualize

function truncatetoregion(crystalmap::CrystalFockMap, region::Region)
    crystalfock = crystalmap|>getoutspace
    crystal = crystalfock|>getcrystal
    truncatefock = getregionfock(crystalfock, region)
    transform = fourier(crystalfock, truncatefock)

    localcorrelations = transform' * crystalmap * transform / (crystal|>vol)
    modecombinations = Iterators.product(localcorrelations|>getoutspace, localcorrelations|>getoutspace)

    getconn(from::Mode, to::Mode) = (
        getattr(to, :r) - getattr(from, :r), from|>removeattr(:r), to|>removeattr(:r))
    # Using a Dict to filter out unique (:r, :b) -> (:r', :b') connections and get the corresponding indices 
    # that we wish to kept.
    modebonds = Dict(
        getconn(frommode, tomode)=>(frommode, tomode) for (frommode, tomode) in modecombinations)
    pruningindices = (index for (_, index) in modebonds)
    crystalpruned = transform * extractindices(localcorrelations, pruningindices)
    return crystalpruned .* crystalpruned'
end

truncationregion = sum((stripspectrum|>getcrystal|>getunitcell).+normalvector*n for n in 0:2)
truncatedstripcorrelations = truncatetoregion(stripcorrelations, truncationregion)
truncatedstripspectrum = truncatedstripcorrelations|>crystalspectrum
truncatedstripspectrum|>linespectrum|>visualize

seedingbands = groupbands(truncatedstripspectrum, bandgrouping=[4])[1]
seedingbands|>linespectrum|>visualize

seedingprojector = seedingbands|>crystalprojector
seedingcorrelations = 1 - seedingprojector

wannierregion = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.875, radius=0.5, minbottomheight=0.15)
wannierfock = getregionfock(seedingcorrelations|>getoutspace, wannierregion)
restrict = fourier(seedingcorrelations|>getoutspace, wannierfock)
localcorrelations = restrict'*seedingcorrelations*restrict/(stripcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates1 = getregionstates(localcorrelations=localcorrelations, grouping=[4])[1]
localstates1 = m45 * localstates1

function snap2unitcell(regionfock::RegionFock)
    spatials = ((mode|>getpos, mode) for mode in regionfock)
    spatials = ((r, r|>basispoint, mode) for (r, mode) in spatials)
    spatials = ((r-b, b, mode) for (r, b, mode) in spatials)
    return (m|>setattr(:r=>r, :b=>b) for (r, b, m) in spatials)|>mapmodes(m->m)|>RegionFock
end

localstates1|>getoutspace|>snap2unitcell|>showmodes

wannierregion = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.875, radius=0.5, minbottomheight=0.15) .- normalvector*0.5
wannierfock = getregionfock(seedingcorrelations|>getoutspace, wannierregion)
restrict = fourier(seedingcorrelations|>getoutspace, wannierfock)
localcorrelations = restrict'*seedingcorrelations*restrict/(stripcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates2 = getregionstates(localcorrelations=localcorrelations, grouping=[2])[1]

localstates = localstates1 + localstates2
visualize(localstates1, markersize=5, logscale=0.5)
localseeds = localstates|>FockMap
transform = fourier(seedingcorrelations|>getoutspace, localseeds|>getoutspace)
crystalseeds = transform * localseeds
crystalseeds = Dict(k=>crystalseeds[k, :] for k in stripcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

wanniercourierisometry = wannierprojection(
    crystalisometries=seedingbands|>geteigenvectors, crystal=stripcrystal, crystalseeds=crystalseeds)

striprestrict = extendedrestrict * blockedcrystalfock

function blockdiag(fockmap::CrystalFockMap)
    @assert hassamespan(fockmap|>getoutspace, fockmap|>getinspace)
    blocks = Dict((ok, ik)=>block for ((ok, ik), block) in fockmap.blocks if ok == ik)
    crystal = fockmap|>getoutspace|>getcrystal
    return CrystalFockMap(crystal, crystal, blocks)
end

courierprojector = striprestrict' * wanniercourierisometry * wanniercourierisometry' * striprestrict
courierprojector = courierprojector|>blockdiag

courierprojector * blockedcorrelations * courierprojector |>crystalspectrum|>visualize

1-blockedcorrelations|>crystalspectrum|>visualize
FockMap(courierprojector)[1:24, 1:24]|>visualize
A = FockMap(courierprojector) * FockMap(1 - blockedcorrelations) * FockMap(courierprojector)
A = CrystalFockMap(A|>getoutspace|>getcrystal, A|>getoutspace|>getcrystal, Dict((k, k)=>A[subspace, subspace] for (k, subspace) in A|>getoutspace|>crystalsubspaces))
A|>crystalspectrum|>visualize

righttransform = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
courierlocalstates = transform' * wanniercourierisometry * righttransform
visualize(courierlocalstates|>RegionState|>normalize, markersize=5, logscale=0.9)

translateregion = getregionfock(couriercorrelations|>getoutspace, courierlocalstates|>getoutspace|>getregion)|>getregion
remapper = spatialremapper(courierlocalstates|>getoutspace, translateregion)
courierlocalstates = remapper * courierlocalstates
transform = fourier(couriercorrelations|>getoutspace, courierlocalstates|>getoutspace) / (couriercrystal|>vol|>sqrt)
courierisometry = transform * courierlocalstates
courierprojector = courierisometry .* courierisometry'
courierprojector|>crystalspectrum|>visualize

courierspectrum = couriercorrelations|>crystalspectrum|>roundingpurification
purifiedcorrelations = courierspectrum|>CrystalFockMap

courierprojector * couriercorrelations * courierprojector |>crystalspectrum
