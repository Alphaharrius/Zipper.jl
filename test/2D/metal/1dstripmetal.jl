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
normalvector = metalspace*(-2, 2)
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

region = getcrosssection(crystal=metalspectrum|>getcrystal, normalvector=normalvector/2, radius=1.5)
region|>visualize
regionfock = getregionfock(stripcorrelations|>getoutspace, region)
transform = fourier(stripcorrelations|>getoutspace, regionfock)
localcorrelations = transform'*stripcorrelations*transform/(stripspectrum|>getcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]
localstates = m135*localstates
visualize(localstates|>normalize, markersize=5, logscale=0.5)

localseeds = localstates

region = region .+ normalvector/2
region|>visualize
regionfock = getregionfock(stripcorrelations|>getoutspace, region)
transform = fourier(stripcorrelations|>getoutspace, regionfock)
localcorrelations = transform'*stripcorrelations*transform/(stripspectrum|>getcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]
localstates = m135*localstates
visualize(localstates|>normalize, markersize=5, logscale=0.5)

localseeds += localstates

localseeds = localseeds|>FockMap
visualize(localseeds|>RegionState, markersize=5, logscale=0.5)
transform = fourier(stripcorrelations|>getoutspace, localseeds|>getoutspace)/(stripspectrum|>getcrystal|>vol)
crystalseeds = transform * localseeds
crystalseeds = Dict(k=>crystalseeds[k, :] for k in stripspectrum|>getcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing strip filled states..."
wannierfilledisometry = wannierprojection(
    crystalisometries=filledbands|>geteigenvectors, crystal=stripspectrum|>getcrystal, crystalseeds=crystalseeds)
wanniercrystal = wannierfilledisometry|>getinspace|>getcrystal
rightrestrict = (
    fourier(wannierfilledisometry|>getinspace, wannierfilledisometry|>getinspace|>unitcellfock|>RegionFock) / (wanniercrystal|>vol|>sqrt))
wannierlocalisometry = transform' * wannierfilledisometry * rightrestrict
stripfilledlocalstates = wannierlocalisometry|>RegionState|>normalize
visualize(stripfilledlocalstates|>normalize, markersize=5, logscale=0.9)

remapper = spatialremapper(wannierlocalisometry|>getoutspace, metalcorrelations|>getoutspace)
filledisometry = remapper * wannierlocalisometry
transform = fourier(metalcorrelations|>getoutspace, filledisometry|>getoutspace)
wanniercrystalisometry = transform * filledisometry
wanniermetalprojector = wanniercrystalisometry .* wanniercrystalisometry'
wanniermetalprojector|>crystalspectrum|>visualize

groupbands(wanniermetalprojector|>crystalspectrum, :filled=>(v->v<1e-5))
