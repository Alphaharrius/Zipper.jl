using LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Zipper, Plots

setmaxthreads(Threads.nthreads())

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [48, 48])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1 = members(modes)

t_a = 1
t_b = -1
tₙ = -1 + 0im

onsite = [
    (m1, m1) => t_b,
    (m0, m0) => t_a
]

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:r => [-1, 0] ∈ triangular)) => tₙ,
    (m0, m1 |> setattr(:r => [0, 1] ∈ triangular)) => tₙ]

bonds::FockMap = bondmap([onsite..., nearestneighbor...])

energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstateprojector = groundstates|>crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector

@info("Starting RG...")
crystalfock = correlations|>getoutspace

scale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
@info("Performing blocking...")
@info("Generating blocking transformation...")
blocker = @time scale * crystalfock
@info("Performing blocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'
blockedcrystalfock = blockedcorrelations|>getoutspace
blockedcrystal::Crystal = blockedcrystalfock|>getcrystal
blockedspace::RealSpace = blockedcrystal|>getspace

@info "Computing frozen restricter..."
distillregion::Region = getsphericalregion(crystal=blockedcrystal, radius=0.9, metricspace=blockedspace|>orthospace)
visualize(distillregion)
# refrot = inv([2/3 -1/3; -1/3 2/3]')
# center = [0,0] ∈ blockedspace
# distillregion::Region = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=center, metricspace=blockedspace)
frozenseedingfock::RegionFock = quantize(distillregion, 1)
frozenrestrict = fourier(blockedcrystalfock, frozenseedingfock)

@info "Computing local correlations..."
localrestrict = fourier(blockedcrystalfock, frozenseedingfock) / (blockedcrystal|>vol|>sqrt)
localcorrelations = localrestrict' * blockedcorrelations * localrestrict
localspectrum = localcorrelations|>eigspech
visualize(localspectrum)
filledfock = FockSpace(m for (m, v) in localspectrum|>geteigenvalues if v < 0.000007)
emptyfock = FockSpace(m for (m, v) in localspectrum|>geteigenvalues if v > (1-0.000007))

@info "Computing frozen projectors..."
localfilledisometry = geteigenvectors(localspectrum)[:, filledfock]
crystalfilledisometry = frozenrestrict * localfilledisometry
filledprojector = crystalfilledisometry .* crystalfilledisometry'

localemptyisometry = geteigenvectors(localspectrum)[:, emptyfock]
crystalemptyisometry = frozenrestrict * localemptyisometry
emptyprojector = crystalemptyisometry .* crystalemptyisometry'

@info "Computing global distiller..."
globaldistiller = emptyprojector - filledprojector

globaldistiller|>eigspech|>visualize

@info "Distilling..."
distilled = groupbands(globaldistiller|>crystalspectrum, :filled=>(v->v < -1e-7), :empty=>(v->v > 1e-7))

@info "Searching seed for filled bands..."
filledbands = distilled[:filled]
filledprojector = crystalprojector(filledbands)
filledcorrelations = idmap(filledprojector|>getoutspace) - filledprojector
filledlocalcorrelations = localrestrict' * filledcorrelations * localrestrict
filledlocalcorrelations|>eigspech|>visualize
filledseeds = getregionstates(localcorrelations=filledlocalcorrelations, grouping=[3])[1]
filledseeds = c3 * filledseeds

@info "Searching seed for empty bands..."
emptybands = distilled[:empty]
emptyprojector = crystalprojector(emptybands)
emptycorrelations = idmap(emptyprojector|>getoutspace) - emptyprojector
emptylocalcorrelations = localrestrict' * emptycorrelations * localrestrict
emptylocalcorrelations|>eigspech|>visualize
emptyseeds = getregionstates(localcorrelations=emptylocalcorrelations, grouping=[3])[1]
emptyseeds = c3 * emptyseeds

@info "Wannierizing filled bands..."
crystalfilledrseeds = crystalisometries(localisometry=filledseeds|>FockMap, crystalfock=blockedcrystalfock, addinspacemomentuminfo=true)
wannierfilledisometry = wannierprojection(
    crystalisometries=distilled[:filled]|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledrseeds|>Dict)

visualregion = getsphericalregion(crystal=blockedcrystal, radius=3, metricspace=blockedspace|>orthospace)
visualfock = quantize(visualregion, 1)

@info "Computing filled correlations..."
filledcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
filledcorrelations|>eigspech|>visualize