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
distillfock = fockspaceunion(g*distillfock|>getoutspace for g in pointgroupelements(c4))|>RegionFock
transform = fourier(crystalfock, distillfock)
localcorrelations = transform'*correlations*transform/(crystal|>vol|>sqrt)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[2, 44, 2])

localfilledisometry = localstates[1]|>FockMap
crystalfilledisometry = transform * localfilledisometry
filledprojector = crystalfilledisometry .* crystalfilledisometry'
filledprojector|>crystalspectrum|>visualize

localemptyisometry = localstates[3]|>FockMap
crystalemptyisometry = transform * localemptyisometry
emptyprojector = crystalemptyisometry .* crystalemptyisometry'
emptyprojector|>crystalspectrum|>visualize

globaldistiller = emptyprojector - filledprojector
distillerspectrum = globaldistiller|>crystalspectrum
distillerspectrum|>visualize

distilled = distillation(distillerspectrum, :courier=>(v -> v > -1e-5 && v < 1e-5))

frozenbands = distilled[:others]
frozenprojector = frozenbands|>crystalprojector
frozencorrelations = 1 - frozenprojector
regionfock = crystalfock|>unitcellfock|>RegionFock
transform = fourier(frozencorrelations|>getoutspace, regionfock)
localcorrelations = transform'*frozencorrelations*transform/(crystal|>vol)
localcorrelations|>eigspech|>visualize
localstate = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]
localstates = sum(g.*localstate for g in c4|>pointgroupelements)

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

    wannierfrozenisometry|>getinspace|>unitcellfock|>showmodes

visualregion = getsphericalregion(crystal=crystal, radius=4, metricspace=space)
visualfock = quantize(visualregion, 3)
lefttransform = fourier(wannierfrozenisometry|>getoutspace, visualfock) / (crystal|>vol|>sqrt)
righttransform = fourier(wannierfrozenisometry|>getinspace, wannierfrozenisometry|>getinspace|>unitcellfock|>RegionFock)
frozenlocalstates = lefttransform' * wannierfrozenisometry * righttransform
visualize(frozenlocalstates|>RegionState, markersize=1, logscale=0.7)

courierbands = distilled[:courier]
courierprojector = courierbands|>crystalprojector
couriercorrelations = 1 - courierprojector
regionfock = crystalfock|>unitcellfock|>RegionFock
transform = fourier(couriercorrelations|>getoutspace, regionfock)
localcorrelations = transform'*couriercorrelations*transform/(crystal|>vol)
localcorrelations|>eigspech|>visualize
localstate = getregionstates(localcorrelations=localcorrelations, grouping=[2])[1]
localstates = sum(g.*localstate for g in c4|>pointgroupelements)
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

visualregion = getsphericalregion(crystal=crystal, radius=2, metricspace=space)
visualfock = quantize(visualregion, 2)
lefttransform = fourier(wanniercourierisometry|>getoutspace, visualfock) / (crystal|>vol|>sqrt)
righttransform = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
courierlocalstates = lefttransform' * wanniercourierisometry * righttransform
visualize(courierlocalstates|>RegionState, markersize=1, logscale=0.7)

fiodir("/Users/alphaharrius/ZERData/squaremetal/rg1/branch")
metalcorrelations = fioload("metalcorrelations")

frozencorrelations = wannierfrozenisometry' * metalcorrelations * wannierfrozenisometry
frozencorrelations|>crystalspectrum|>visualize

couriercorrelations = wanniercourierisometry' * metalcorrelations * wanniercourierisometry
couriercorrelations|>crystalspectrum|>visualize

@info "Performing blocking..."
couriercrystal = couriercorrelations|>getoutspace|>getcrystal
courierspace = couriercrystal|>getspace

scale = Scale([2 0; 0 2], courierspace)
block = scale * (couriercorrelations|>getoutspace)
blockedcorrelations = block * couriercorrelations * block'
blockedcrystalfock = blockedcorrelations|>getoutspace
blockedcrystal = blockedcrystalfock|>getcrystal

function getregionfock(crystalfock::CrystalFock, region::Region)::RegionFock
    unitcelloffsets = Subset(getattr(m, :b) - basispoint(getattr(m, :b)) for m in crystalfock|>unitcellfock)
    if length(unitcelloffsets) != 1 || first(unitcelloffsets) != unitcelloffsets|>getspace|>getorigin
        error("Unusual unit cell detected!")
    end

    offsets = (r-basispoint(r) for r in region)
    return RegionFock(m+r for (m, r) in Iterators.product(crystalfock|>unitcellfock|>RegionFock, offsets))
end

function Base.:*(restrict::ExtendedRestrict, crystalfock::CrystalFock)
    crystal::Crystal = crystalfock|>getcrystal
    scaledcrystal::Crystal = (restrict.scale) * crystal
    scaledspace::RealSpace = scaledcrystal|>getspace
    scaledkspace::MomentumSpace = convert(MomentumSpace, scaledspace)
    # This is the unit cell of the 1-D embedded crystal defining the strip region in the blocked crystal metric space.
    unscaledextendedunitcell::Region = getcrosssection(crystal=crystal, normalvector=restrict.normalvector, radius=restrict.stripradius)
    unscaledextendedhomefock = getregionfock(crystalfock, unscaledextendedunitcell)

    bz::Subset{Momentum} = crystal|>brillouinzone
    momentummappings::Base.Generator = (basispoint(scaledkspace * k) => k for k in bz)

    restrictedfourier = fourier(crystalfock, unscaledextendedhomefock)'
    blocks::Dict = Dict()

    stripunitcell::Region = Subset(scaledspace * r for r in unscaledextendedunitcell)
    stripcrystal::Crystal = Crystal(stripunitcell, scaledcrystal|>size)
    volumeratio::Real = vol(crystal) / vol(stripcrystal)

    scaledksubspaces::Dict{Momentum, FockSpace} = Dict()
    for (scaledk, k) in momentummappings
        kfourier::FockMap = restrictedfourier[:, k] / sqrt(volumeratio)
        if !haskey(scaledksubspaces, scaledk)
            scaledksubspaces[scaledk] = FockSpace(
                setattr(mode, :k=>scaledk, :b=>(scaledspace * getpos(mode)))|>removeattr(:r) for mode in kfourier|>getoutspace)
        end
        blocks[(scaledk, k)] = FockMap(scaledksubspaces[scaledk], kfourier|>getinspace, kfourier|>rep)
    end

    return CrystalFockMap(stripcrystal, crystal, blocks)
end

@info "Performing extended restrictions..."
normalvector = (1, 1) ∈ (blockedcrystal|>getspace)
extendedrestrict = extendedcrystalrestrict(
    crystal=blockedcrystal, normalvector=normalvector, stripradius=0.5)

blockedcrystal|>getunitcell|>visualize

getcrosssection(crystal=blockedcrystal, normalvector=extendedrestrict.normalvector, radius=extendedrestrict.stripradius)|>visualize

restrict = extendedrestrict * blockedcrystalfock
stripcorrelations = restrict * blockedcorrelations * restrict'
stripcorrelations|>getoutspace|>unitcellfock|>RegionFock|>getregion|>visualize
stripspectrum = stripcorrelations|>crystalspectrum
stripspectrum|>linespectrum|>visualize

stripcorrelations|>getoutspace|>unitcellfock|>showmodes

