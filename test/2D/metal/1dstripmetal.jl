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
stripcorrelations = fioload("stripmetalcorrelations")
stripspectrum = roundingpurification(stripcorrelations|>crystalspectrum, tolerance=0.4)
stripspectrum|>visualize
stripcorrelations = stripspectrum|>CrystalFockMap
stripcorrelations|>crystalspectrum|>visualize

@info "Rotating correlations..."
scale = Scale([1 -1; 1 1], stripspectrum|>getcrystal|>getspace)
rotate = scale*getoutspace(stripcorrelations)
rotate|>getoutspace|>getcrystal|>getunitcell|>visualize
rotatedcorrelations = rotate*stripcorrelations*rotate'
rotatedspectrum = rotatedcorrelations|>crystalspectrum
rotatedspectrum|>visualize

normalvector = rotatedspectrum|>getcrystal|>getspace|>getbasisvectors|>first
unitcell = rotatedspectrum|>getcrystal|>getunitcell
distillregion = unitcell + (unitcell.+normalvector)
regionfock = getregionfock(rotatedcorrelations|>getoutspace, distillregion)
restrict = fourier(rotatedcorrelations|>getoutspace, regionfock)
localcorrelations = restrict'*rotatedcorrelations*restrict/(rotatedspectrum|>getcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1, 48-2, 1])
localfrozen = (localstates[1]+localstates[3])|>FockMap
transform = fourier(rotatedcorrelations|>getoutspace, localfrozen|>getoutspace)
crystalfrozen = transform * localfrozen
frozenprojector = crystalfrozen .* crystalfrozen'
frozenspectrum = frozenprojector|>crystalspectrum
frozenspectrum|>visualize

bands = groupbands(frozenspectrum, :frozen=>(v->v>1e-5))
courierbands = bands[:others]
courierprojector = courierbands|>crystalprojector

courierprojector*rotatedcorrelations*courierprojector|>crystalspectrum|>visualize

visualize(stripspectrum|>getcrystal|>brillouinzone, rotatedspectrum|>getcrystal|>brillouinzone, scale*(rotatedspectrum|>getcrystal)|>brillouinzone)
visualize(rotatedspectrum|>getcrystal|>brillouinzone, scale*(rotatedspectrum|>getcrystal)|>brillouinzone)

@info "Blocking correlations..."
scale = Scale([2 0; 0 1], rotatedspectrum|>getcrystal|>getspace)
block = scale*getoutspace(rotatedcorrelations)
blockedcorrelations = block*rotatedcorrelations*block'
blockedspectrum = blockedcorrelations|>crystalspectrum
blockedspectrum|>visualize

blockedcorrelations|>getoutspace|>unitcellfock

normalvector = blockedspectrum|>getcrystal|>getspace|>getbasisvectors|>first
unitcell = blockedspectrum|>getcrystal|>getunitcell
unitcell|>visualize
distillregion = unitcell + (unitcell.+normalvector)
distillregion|>visualize
regionfock = getregionfock(blockedcorrelations|>getoutspace, distillregion)
restrict = fourier(blockedcorrelations|>getoutspace, regionfock)
localcorrelations = restrict'*blockedcorrelations*restrict/(blockedspectrum|>getcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[2, 96-4, 2])

localfrozen = (localstates[1]+localstates[3])|>FockMap
transform = fourier(blockedcorrelations|>getoutspace, localfrozen|>getoutspace)
crystalfrozen = transform * localfrozen
frozenprojector = crystalfrozen .* crystalfrozen'
frozenspectrum = frozenprojector|>crystalspectrum
frozenspectrum|>visualize

bands = groupbands(frozenspectrum, :frozen=>(v->v>0.1))
courierbands = bands[:others]
courierprojector = courierbands|>crystalprojector

courierprojector*blockedcorrelations*courierprojector|>crystalspectrum|>visualize
