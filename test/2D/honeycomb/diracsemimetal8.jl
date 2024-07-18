using LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Zipper, Plots

setmaxthreads(Threads.nthreads())

function focktraceL1norm(fockmap,systemsize)
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end

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

t_a = 0
t_b = 0
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

scale = Scale([8 0; 0 8], crystalfock|>getcrystal|>getspace)
@info("Performing blocking...")
@info("Generating blocking transformation...")
blocker = @time scale * crystalfock
@info("Performing blocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'
blockedcrystalfock = blockedcorrelations|>getoutspace
blockedcrystal::Crystal = blockedcrystalfock|>getcrystal
blockedspace::RealSpace = blockedcrystal|>getspace
blockedunitcell = blockedcrystal|>getunitcell

@info "Computing frozen restricter..."
distillregion = blockedunitcell+c3*blockedunitcell+(c3)^2*blockedunitcell
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
sortedevalpair = sort([(m,v) for (m, v) in localspectrum|>geteigenvalues],by=x->x[2])
filledfock = FockSpace(m for (m,v) in sortedevalpair[1:48])
emptyfock = FockSpace(m for (m,v) in sortedevalpair[end-47:end])

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
centerfilledseeds = getregionstates(localcorrelations=filledlocalcorrelations, grouping=[48])[1]
# filledseeds = c6 * filledseeds

@info "Searching seed for empty bands..."
emptybands = distilled[:empty]
emptyprojector = crystalprojector(emptybands)
emptycorrelations = idmap(emptyprojector|>getoutspace) - emptyprojector
emptylocalcorrelations = localrestrict' * emptycorrelations * localrestrict
# emptylocalcorrelations|>eigspech|>visualize
centeremptyseeds = getregionstates(localcorrelations=emptylocalcorrelations, grouping=[48])[1]
# emptyseeds = c6 * emptyseeds

# @info "Searching seed at site A and B..."
# siteAseedingregion = intersect(blockedunitcell,(c6)^5*blockedunitcell)
# # visualize(siteAseedingregion)
# siteAseedingfock::RegionFock = quantize(siteAseedingregion|>Subset, 1)
# siteArestrict = fourier(blockedcrystalfock, siteAseedingfock) / (blockedcrystal|>vol|>sqrt)
# siteAlocalfilledcorrelations = siteArestrict' * filledcorrelations * siteArestrict
# # siteAlocalfilledcorrelations|>eigspech|>visualize
# siteAlocalemptycorrelations = siteArestrict' * emptycorrelations * siteArestrict
# # siteAlocalemptycorrelations|>eigspech|>visualize

# siteAfilledseed = getregionstates(localcorrelations=siteAlocalfilledcorrelations, grouping=[3])[1]
# siteAemptyseed = getregionstates(localcorrelations=siteAlocalemptycorrelations, grouping=[3])[1]
# # siteAfilledseed = c3 * siteAfilledseed
# # siteAemptyseed = c3 * siteAemptyseed

# siteBseedingregion = intersect(blockedunitcell,(c6)*blockedunitcell)
# siteBseedingfock::RegionFock = quantize(siteBseedingregion|>Subset, 1)
# siteBrestrict = fourier(blockedcrystalfock, siteBseedingfock) / (blockedcrystal|>vol|>sqrt)
# siteBlocalfilledcorrelations = siteBrestrict' * filledcorrelations * siteBrestrict
# # siteBlocalfilledcorrelations|>eigspech|>visualize
# siteBlocalemptycorrelations = siteBrestrict' * emptycorrelations * siteBrestrict
# # siteBlocalemptycorrelations|>eigspech|>visualize

# siteBfilledseed = getregionstates(localcorrelations=siteBlocalfilledcorrelations , grouping=[3])[1]
# siteBemptyseed = getregionstates(localcorrelations=siteBlocalemptycorrelations, grouping=[3])[1]
# # siteBfilledseed = c3 * siteBfilledseed
# # siteBemptyseed = c3 * siteBemptyseed

filledseeds = centerfilledseeds
emptyseeds = centeremptyseeds

@info "Wannierizing filled bands..."
crystalfilledrseeds = crystalisometries(localisometry=filledseeds|>FockMap, crystalfock=blockedcrystalfock, addinspacemomentuminfo=true)
wannierfilledisometry = wannierprojection(
    crystalisometries=distilled[:filled]|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledrseeds|>Dict)

visualregion = getsphericalregion(crystal=blockedcrystal, radius=3, metricspace=blockedspace|>orthospace)
visualfock = quantize(visualregion, 1)

@info "Computing local filled states..."
leftrestrict = fourier(wannierfilledisometry|>getoutspace, visualfock) / (blockedcrystal|>vol|>sqrt)
rightrestrict = fourier(wannierfilledisometry|>getinspace, wannierfilledisometry|>getinspace|>unitcellfock|>RegionFock)
wannierfilledstates = leftrestrict' * wannierfilledisometry * rightrestrict
# wannierfilledstates|>visualize

@info "Wannierizing empty bands..."
crystalemptyrseeds = crystalisometries(localisometry=emptyseeds|>FockMap, crystalfock=blockedcrystalfock, addinspacemomentuminfo=true)
wannieremptyisometry = wannierprojection(
        crystalisometries=distilled[:empty]|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyrseeds|>Dict)

@info "Computing local empty states..."
leftrestrict = fourier(wannieremptyisometry|>getoutspace, visualfock) / (blockedcrystal|>vol|>sqrt)
rightrestrict = fourier(wannieremptyisometry|>getinspace, wannieremptyisometry|>getinspace|>unitcellfock|>RegionFock)
wannieremptystates = leftrestrict' * wannieremptyisometry * rightrestrict

@info "Computing filled correlations..."
filledcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
filledcorrelations|>eigspech|>visualize

@info "Computing empty correlations..."
emptycorrelations = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry
emptycorrelations|>eigspech|>visualize


@info "Searching seed for courier bands..."
courierbands = distilled[:others]
courierprojector = crystalprojector(courierbands)
couriercorrelations = idmap(courierprojector|>getoutspace) - courierprojector
courierlocalcorrelations = localrestrict' * couriercorrelations * localrestrict
courierlocalcorrelations|>eigspech|>visualize
siteAseedingregion = intersect(blockedunitcell,(c6)^5*blockedunitcell)
siteAseedingfock::RegionFock = quantize(siteAseedingregion|>Subset, 1)
siteArestrict = fourier(blockedcrystalfock, siteAseedingfock) / (blockedcrystal|>vol|>sqrt)
siteBseedingregion = intersect(blockedunitcell,(c6)*blockedunitcell)
siteBseedingfock::RegionFock = quantize(siteBseedingregion|>Subset, 1)
siteBrestrict = fourier(blockedcrystalfock, siteBseedingfock) / (blockedcrystal|>vol|>sqrt)
siteAlocalcouriercorrelations = siteArestrict' * couriercorrelations * siteArestrict
siteAlocalcouriercorrelations|>eigspech|>visualize
siteBrestrict = fourier(blockedcrystalfock, siteBseedingfock) / (blockedcrystal|>vol|>sqrt)
siteBlocalcouriercorrelations = siteBrestrict' * couriercorrelations * siteBrestrict
siteBlocalcouriercorrelations|>eigspech|>visualize

siteAcourierseed = getregionstates(localcorrelations=siteAlocalcouriercorrelations, grouping=[16])[1]
siteBcourierseed = getregionstates(localcorrelations=siteBlocalcouriercorrelations, grouping=[16])[1]
courierseeds = siteAcourierseed+siteBcourierseed
crystalcourierseeds = crystalisometries(localisometry=courierseeds|>FockMap, crystalfock=blockedcrystalfock, addinspacemomentuminfo=true)
wanniercourierisometry = wannierprojection(
        crystalisometries=distilled[:others]|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds|>Dict)


@info "Computing local courier states..."
leftrestrict = fourier(wanniercourierisometry|>getoutspace, visualfock) / (blockedcrystal|>vol|>sqrt)
rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
wanniercourierstate = leftrestrict' * wanniercourierisometry * rightrestrict

(wanniercourierstate|>RegionState)[6]|>visualize

couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
couriercorrelations|>eigspech|>visualize
couriercorrelations|>getinspace|>getcrystal|>getunitcell|>visualize