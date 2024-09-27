using LinearAlgebra
using Zipper, Plots

plotlyjs()
# usecrystaldensemap()
setmaxthreads(Threads.nthreads())

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
# spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

unitcell = Subset(pa, pb)
visualize(unitcell)
crystal = Crystal(unitcell, [32, 32])
space = crystal|>getspace
# reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1 = members(modes)

t_a = 0 + 0im
t_b = -0 + 0im
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
sort([(mode,abs(eval)) for (mode,eval) in energyspectrum|>geteigenvalues],by=x->x[2])

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstateprojector = groundstates|>crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector

crystalfock = correlations|>getoutspace
localrestrict = fourier(crystalfock, quantize(unitcell,1)) / (crystal|>vol|>sqrt)
localcorrelations = localrestrict'*correlations*localrestrict
[(localcorrelations[mode,mode]|>rep)[1,1]|>abs for mode in quantize(unitcell,1)]


function zer(correlations)
    @info("Starting RG...")
    crystalfock = correlations|>getoutspace

    scale = Scale([2 0; 0 2], crystalfock|>getcrystal|>getspace)
    @info("Performing blocking...")
    @info("Generating blocking transformation...")
    block = @time scale * crystalfock
    @info("Performing blocking on correlations...")
    blockedcorrelations = @time block * correlations * block'
    blockedcrystalfock = blockedcorrelations|>getoutspace
    blockedcrystal::Crystal = blockedcrystalfock|>getcrystal
    blockedspace::RealSpace = blockedcrystal|>getspace

    @info "Computing frozen restricter..."
    distillregion::Region = getsphericalregion(crystal=blockedcrystal, radius=1, metricspace=blockedspace|>orthospace)
    frozenseedingfock::RegionFock = quantize(distillregion, 1)
    frozenrestrict = fourier(blockedcrystalfock, frozenseedingfock)

    @info "Computing local correlations..."
    localrestrict = fourier(blockedcrystalfock, frozenseedingfock) / (blockedcrystal|>vol|>sqrt)
    localcorrelations = localrestrict' * blockedcorrelations * localrestrict
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3, 18, 3])

    @info "Computing frozen projectors..."
    localfilledisometry = localstates[1]|>FockMap
    crystalfilledisometry = frozenrestrict * localfilledisometry
    filledprojector = crystalfilledisometry .* crystalfilledisometry'

    localemptyisometry = localstates[3]|>FockMap
    crystalemptyisometry = frozenrestrict * localemptyisometry
    emptyprojector = crystalemptyisometry .* crystalemptyisometry'

    @info "Computing global distiller..."
    globaldistiller = emptyprojector - filledprojector

    @info "Distilling..."
    distilled = distillation(globaldistiller|>crystalspectrum, :filled=>(v->v < -1e-3), :empty=>(v->v > 1e-3))

    @info "Searching seed for filled bands..."
    filledbands = distilled[:filled]
    filledprojector = crystalprojector(filledbands)
    filledcorrelations = idmap(filledprojector|>getoutspace) - filledprojector
    filledlocalcorrelations = localrestrict' * filledcorrelations * localrestrict
    filledlocalcorrelations|>eigspech|>visualize
    filledseeds = getregionstates(localcorrelations=filledlocalcorrelations, grouping=[3])[1]
    filledseeds = c6 * filledseeds

    @info "Searching seed for empty bands..."
    emptybands = distilled[:empty]
    emptyprojector = crystalprojector(emptybands)
    emptycorrelations = idmap(emptyprojector|>getoutspace) - emptyprojector
    emptylocalcorrelations = localrestrict' * emptycorrelations * localrestrict
    emptylocalcorrelations|>eigspech|>visualize
    emptyseeds = getregionstates(localcorrelations=emptylocalcorrelations, grouping=[3])[1]
    emptyseeds = c6 * emptyseeds

    @info "Searching seed at site A and B..."
    siteAseedingcenter::Offset = [1/3, 2/3] ∈ blockedspace
    siteAseedingfock::RegionFock = quantize(siteAseedingcenter|>Subset, 1)
    siteArestrict = fourier(blockedcrystalfock, siteAseedingfock) / (blockedcrystal|>vol|>sqrt)
    siteAlocalcorrelations = siteArestrict' * blockedcorrelations * siteArestrict
    siteAseed = getregionstates(localcorrelations=siteAlocalcorrelations, grouping=[1])[1]
    siteAseed = c3 * siteAseed

    siteBseedingcenter::Offset = [2/3, 1/3] ∈ blockedspace
    siteBseedingfock::RegionFock = quantize(siteBseedingcenter|>Subset, 1)
    siteBrestrict = fourier(blockedcrystalfock, siteBseedingfock) / (blockedcrystal|>vol|>sqrt)
    siteBlocalcorrelations = siteBrestrict' * blockedcorrelations * siteBrestrict
    siteBseed = getregionstates(localcorrelations=siteBlocalcorrelations, grouping=[1])[1]
    siteBseed = c3 * siteBseed

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

    @info "Wannierizing empty bands..."
    crystalemptyrseeds = crystalisometries(localisometry=emptyseeds|>FockMap, crystalfock=blockedcrystalfock, addinspacemomentuminfo=true)
    wannieremptyisometry = wannierprojection(
        crystalisometries=distilled[:empty]|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyrseeds|>Dict)

    @info "Computing local empty states..."
    leftrestrict = fourier(wannieremptyisometry|>getoutspace, visualfock) / (blockedcrystal|>vol|>sqrt)
    rightrestrict = fourier(wannieremptyisometry|>getinspace, wannieremptyisometry|>getinspace|>unitcellfock|>RegionFock)
    wannieremptystates = leftrestrict' * wannieremptyisometry * rightrestrict

    @info "Wannierizing courier bands..."
    courierseeds = siteAseed + siteBseed
    crystalcourierseeds = crystalisometries(localisometry=courierseeds|>FockMap, crystalfock=blockedcrystalfock, addinspacemomentuminfo=true)
    wanniercourierisometry = wannierprojection(
        crystalisometries=distilled[:others]|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds|>Dict)

    @info "Computing local courier states..."
    leftrestrict = fourier(wanniercourierisometry|>getoutspace, visualfock) / (blockedcrystal|>vol|>sqrt)
    rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
    wanniercourierstate = leftrestrict' * wanniercourierisometry * rightrestrict

    @info "Computing filled correlations..."
    filledcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
    @info "Computing empty correlations..."
    emptycorrelations = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry
    @info "Computing courier correlations..."
    couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
    display(couriercorrelations|>eigspech|>visualize)

    couriercorrelationspectrum = couriercorrelations|>crystalspectrum
    @info "Purifiying courier correlations..."
    purifiedcorrelations = roundingpurification(couriercorrelationspectrum)|>crystalfockmap

    return Dict(
        :block=>block,
        :globaldistiller=>globaldistiller,
        :couriercorrelations=>purifiedcorrelations,
        :filledcorrelations=>filledcorrelations,
        :emptycorrelations=>emptycorrelations,
        :wanniercourierisometry=>wanniercourierisometry,
        :wannierfilledisometry=>wannierfilledisometry,
        :wannieremptyisometry=>wannieremptyisometry,
        :wannierfilledstates=>wannierfilledstates|>RegionState,
        :wannieremptystates=>wannieremptystates|>RegionState,
        :wanniercourierstates=>wanniercourierstate|>RegionState,
        :rawcouriercorrelations=>couriercorrelations)
end

region::Region = getsphericalregion(crystal=crystal, radius=2, metricspace=space|>orthospace)
seedingfock::RegionFock = quantize(region, 1)
rscorrelations = regioncorrelations(correlations,seedingfock)
rscorrelations|>eigspec|>visualize

rg1 = @time zer(correlations)

rg1correlations = rg1[:couriercorrelations]
rg1correlations|>getoutspace|>getcrystal
rgedcrystal = rg1correlations|>getoutspace|>getcrystal
rgedspace = rgedcrystal|>getspace

rg1region::Region = getsphericalregion(crystal=rgedcrystal, radius=2, metricspace=rgedspace|>orthospace)
rg1seedingfock::RegionFock = quantize(rg1region, 1)
rg1rscorrelations = regioncorrelations(rg1correlations,rg1seedingfock)
rg1rscorrelations|>eigspec|>visualize

rg2 = @time zer(rg1correlations)

rg2correlations = rg2[:couriercorrelations]
rg2correlations|>getoutspace|>getcrystal
rg2edcrystal = rg2correlations|>getoutspace|>getcrystal
rg2edspace = rg2edcrystal|>getspace

rg2region::Region = getsphericalregion(crystal=rg2edcrystal, radius=2, metricspace=rg2edspace|>orthospace)
rg2seedingfock::RegionFock = quantize(rg2region, 1)
rg2rscorrelations = regioncorrelations(rg2correlations,rg2seedingfock)
rg2rscorrelations|>eigspec|>visualize

rg3 = @time zer(rg2correlations)

rg3correlations = rg3[:couriercorrelations]
rg3correlations|>getoutspace|>getcrystal
rg3edcrystal = rg3correlations|>getoutspace|>getcrystal
rg3edspace = rg3edcrystal|>getspace

rg3region::Region = getsphericalregion(crystal=rg3edcrystal, radius=2, metricspace=rg3edspace|>orthospace)
rg3seedingfock::RegionFock = quantize(rg3region, 1)
rg3rscorrelations = regioncorrelations(rg3correlations,rg3seedingfock)
rg3rscorrelations|>eigspec|>visualize

rg4 = @time zer(rg3correlations)
