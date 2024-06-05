using LinearAlgebra
using Zipper, Plots
plotlyjs()

usecrystaldensemap()
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
bc = BoundaryCondition(triangular, 192, 192)
crystal = Crystal(unitcell, bc)
reciprocalhashcalibration(bc.bounds)

fiodir("/Users/alphaharrius/ZERData/chern192x192/RG2")
courieriso = fioload("wanniercourierisometry")
corr = fioload("couriercorrelations")
blk = fioload("block")
courierzipper = blk' * courieriso
courierzipper * corr * courierzipper'
fillediso = fioload("wannierfilledisometry")
filledzipper = blk' * fillediso
emptyiso = fioload("wannieremptyisometry")
emptyzipper = blk' * emptyiso
region = buildregion(crystal, 10, 10)
region = region + c3*region + c3*c3*region
region|>visualize
regionfock = getregionfock(blk|>getinspace, region)
transform = fourier(blk|>getinspace, regionfock)
rt = fourier(emptyzipper|>getinspace, emptyzipper|>getinspace|>unitcellfock|>RegionFock)
localmap = transform' * emptyzipper * rt
fiosave(localmap, name="emptylocalstates")
fiosave(fioload("filledlocalstates"), name="filledlocalstatesjson")
visualize(localmap|>RegionState|>normalize, markersize=5, logscale=0.5)

modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1 = members(modes)

t_a = 0 + 0im
t_b = 0 + 0im
tₙ = -1 + 0im
tₕ = 0.4im

onsite = [
    (m1, m1) => t_b,
    (m0, m0) => t_a
]

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:r => [-1, 0] ∈ triangular)) => tₙ,
    (m0, m1 |> setattr(:r => [0, 1] ∈ triangular)) => tₙ]

haldane = [
    (m0, setattr(m0, :r => Point([1, 1], triangular))) => tₕ,
    (m0, setattr(m0, :r => Point([-1, 0], triangular))) => tₕ,
    (m0, setattr(m0, :r => Point([0, -1], triangular))) => tₕ,
    (m1, setattr(m1, :r => Point([1, 1], triangular))) => -tₕ,
    (m1, setattr(m1, :r => Point([-1, 0], triangular))) => -tₕ,
    (m1, setattr(m1, :r => Point([0, -1], triangular))) => -tₕ]

bonds::FockMap = bondmap([onsite..., nearestneighbor..., haldane...])

energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
energyspectrum|>visualize
groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstateprojector = groundstates|>crystalprojector
correlations = 1 - groundstateprojector

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
        :rawcouriercorrelations=>couriercorrelations,
        :courierzipper=>block'*wanniercourierisometry,
        :filledzipper=>block'*wannierfilledisometry,
        :emptyzipper=>block'*wannieremptyisometry)
end

rg1 = @time zer(correlations)
rg1[:block]'*rg1[:wanniercourierisometry]
rg6[:globaldistiller]|>crystalspectrum|>visualize
visualize(rg1[:wanniercourierstates]|>normalize, markersize=5, logscale=0.5)

rg2 = @time zer(rg1[:couriercorrelations])
rg2[:globaldistiller]|>crystalspectrum|>visualize
rg3 = @time zer(rg2[:couriercorrelations])
rg4 = @time zer(rg3[:couriercorrelations])
rg5 = @time zer(rg4[:couriercorrelations])
rg6 = @time zer(rg5[:couriercorrelations])
rg7 = @time zer(rg6[:couriercorrelations])

fiodir("/Users/alphaharrius/ZERData/chern192x192/wannierisometries")
fiodir("/Users/alphaharrius/ZERData/chern192x192/wannierisometries/RG1")
fiosave(rg1[:wanniercourierisometry]|>CrystalFockMap, name="wanniercourierisometry")
fiosave(rg1[:wannierfilledisometry]|>CrystalFockMap, name="wannierfilledisometry")
fiosave(rg1[:wannieremptyisometry]|>CrystalFockMap, name="wannieremptyisometry")
fiodir("/Users/alphaharrius/ZERData/chern192x192/wannierisometries/RG2")
fiosave(rg2[:wanniercourierisometry]|>CrystalFockMap, name="wanniercourierisometry")
fiosave(rg2[:wannierfilledisometry]|>CrystalFockMap, name="wannierfilledisometry")
fiosave(rg2[:wannieremptyisometry]|>CrystalFockMap, name="wannieremptyisometry")
fiodir("/Users/alphaharrius/ZERData/chern192x192/wannierisometries/RG3")
fiosave(rg3[:wanniercourierisometry]|>CrystalFockMap, name="wanniercourierisometry")
fiosave(rg3[:wannierfilledisometry]|>CrystalFockMap, name="wannierfilledisometry")
fiosave(rg3[:wannieremptyisometry]|>CrystalFockMap, name="wannieremptyisometry")
fiodir("/Users/alphaharrius/ZERData/chern192x192/wannierisometries/RG4")
fiosave(rg4[:wanniercourierisometry]|>CrystalFockMap, name="wanniercourierisometry")
fiosave(rg4[:wannierfilledisometry]|>CrystalFockMap, name="wannierfilledisometry")
fiosave(rg4[:wannieremptyisometry]|>CrystalFockMap, name="wannieremptyisometry")
fiodir("/Users/alphaharrius/ZERData/chern192x192/wannierisometries/RG5")
fiosave(rg5[:wanniercourierisometry]|>CrystalFockMap, name="wanniercourierisometry")
fiosave(rg5[:wannierfilledisometry]|>CrystalFockMap, name="wannierfilledisometry")
fiosave(rg5[:wannieremptyisometry]|>CrystalFockMap, name="wannieremptyisometry")
fiodir("/Users/alphaharrius/ZERData/chern192x192/wannierisometries/RG6")
fiosave(rg6[:wanniercourierisometry]|>CrystalFockMap, name="wanniercourierisometry")
fiosave(rg6[:wannierfilledisometry]|>CrystalFockMap, name="wannierfilledisometry")
fiosave(rg6[:wannieremptyisometry]|>CrystalFockMap, name="wannieremptyisometry")

region = buildregion(crystal, 5, 5)
region = region + c3*region + c3*c3*region
region|>visualize

regionfock = getregionfock(rg1[:courierzipper]|>getoutspace, region)
lt = fourier(rg1[:courierzipper]|>getoutspace, regionfock)
rt = fourier(rg1[:courierzipper]|>getinspace, rg1[:courierzipper]|>getinspace|>unitcellfock|>RegionFock)
rg1couriercrystalstate = FockMap(rg1[:courierzipper]|>CrystalFockMap) * FockMap(rt)
rg1courierlocalstate = FockMap(lt') * rg1couriercrystalstate
visualize(rg1courierlocalstate|>RegionState|>normalize, markersize=5, logscale=0.5)
lt = fourier(rg1[:filledzipper]|>getoutspace, regionfock)
rt = fourier(rg1[:filledzipper]|>getinspace, rg1[:filledzipper]|>getinspace|>unitcellfock|>RegionFock)
rg1filledcrystalstate = FockMap(rg1[:filledzipper]|>CrystalFockMap) * FockMap(rt)
rg1filledlocalstate = FockMap(lt') * rg1filledcrystalstate
visualize(rg1filledlocalstate|>RegionState|>normalize, markersize=5, logscale=0.5)
lt = fourier(rg1[:emptyzipper]|>getoutspace, regionfock)
rt = fourier(rg1[:emptyzipper]|>getinspace, rg1[:emptyzipper]|>getinspace|>unitcellfock|>RegionFock)
rg1emptycrystalstate = FockMap(rg1[:emptyzipper]|>CrystalFockMap) * FockMap(rt)
rg1emptylocalstate = FockMap(lt') * rg1emptycrystalstate
visualize(rg1emptylocalstate|>RegionState|>normalize, markersize=5, logscale=0.5)

region|>visualize
rg2courierzipper = FockMap(rg1[:courierzipper]|>CrystalFockMap)*FockMap(rg2[:courierzipper]|>CrystalFockMap)
rg2filledzipper = FockMap(rg1[:courierzipper]|>CrystalFockMap)*FockMap(rg2[:filledzipper]|>CrystalFockMap)
rg2emptyzipper = FockMap(rg1[:courierzipper]|>CrystalFockMap)*FockMap(rg2[:emptyzipper]|>CrystalFockMap)
region = buildregion(crystal, 8, 8)
c3r = c3|>recenter(triangular*(3,3))
c6r = c6|>recenter(triangular*(3,3))
region = region + c3r*region + c3r*c3r*region
region = region + c6r*region
region|>visualize
regionfock = getregionfock(rg2courierzipper|>getoutspace, region)
lt = fourier(rg2courierzipper|>getoutspace, regionfock)
rt = fourier(rg2courierzipper|>getinspace, rg2courierzipper|>getinspace|>unitcellfock|>RegionFock)
rg2couriercrystalstate = FockMap(rg2courierzipper) * FockMap(rt)
rg2courierlocalstate = FockMap(lt') * rg2couriercrystalstate
visualize(rg2courierlocalstate|>RegionState|>normalize, markersize=5, logscale=0.5)
region = buildregion(crystal, 5, 5)
region = region + c3*region + c3*c3*region
lt = fourier(rg2filledzipper|>getoutspace, regionfock)
rt = fourier(rg2filledzipper|>getinspace, rg2filledzipper|>getinspace|>unitcellfock|>RegionFock)
rg2filledcrystalstate = FockMap(rg2filledzipper) * FockMap(rt)
rg2filledlocalstate = FockMap(lt') * rg2filledcrystalstate
visualize(rg2filledlocalstate|>RegionState|>normalize, markersize=5, logscale=0.5)
lt = fourier(rg2emptyzipper|>getoutspace, regionfock)
rt = fourier(rg2emptyzipper|>getinspace, rg2emptyzipper|>getinspace|>unitcellfock|>RegionFock)
rg2emptycrystalstate = FockMap(rg2emptyzipper) * FockMap(rt)
rg2emptylocalstate = FockMap(lt') * rg2emptycrystalstate
visualize(rg2emptylocalstate|>RegionState|>normalize, markersize=5, logscale=0.5)

rg3courierzipper = rg2courierzipper*FockMap(rg3[:courierzipper]|>CrystalFockMap)
rg3filledzipper = rg2courierzipper*FockMap(rg3[:filledzipper]|>CrystalFockMap)
rg3emptyzipper = rg2courierzipper*FockMap(rg3[:emptyzipper]|>CrystalFockMap)
region = buildregion(crystal, 9, 9)
c3r = c3|>recenter(triangular*(3,3))
c6r = c6|>recenter(triangular*(3,3))
region = region + c3r*region + c3r*c3r*region
region = region + c6r*region
region|>visualize
regionfock = getregionfock(rg3courierzipper|>getoutspace, region)
lt = fourier(rg3courierzipper|>getoutspace, regionfock)
rt = fourier(rg3courierzipper|>getinspace, rg3courierzipper|>getinspace|>unitcellfock|>RegionFock)
rg3couriercrystalstate = FockMap(rg3courierzipper) * FockMap(rt)
rg3courierlocalstate = FockMap(lt') * rg3couriercrystalstate
visualize(rg3courierlocalstate|>RegionState|>normalize, markersize=5, logscale=0.5)
region = buildregion(crystal, 6, 6)
region = region + c3*region + c3*c3*region
regionfock = getregionfock(rg3filledzipper|>getoutspace, region)
lt = fourier(rg3filledzipper|>getoutspace, regionfock)
rt = fourier(rg3filledzipper|>getinspace, rg3filledzipper|>getinspace|>unitcellfock|>RegionFock)
rg3filledcrystalstate = FockMap(rg3filledzipper) * FockMap(rt)
rg3filledlocalstate = FockMap(lt') * rg3filledcrystalstate
visualize(rg3filledlocalstate|>RegionState|>normalize, markersize=5, logscale=0.5)
lt = fourier(rg3emptyzipper|>getoutspace, regionfock)
rt = fourier(rg3emptyzipper|>getinspace, rg3emptyzipper|>getinspace|>unitcellfock|>RegionFock)
rg3emptycrystalstate = FockMap(rg3emptyzipper) * FockMap(rt)
rg3emptylocalstate = FockMap(lt') * rg3emptycrystalstate
visualize(rg3emptylocalstate|>RegionState|>normalize, markersize=5, logscale=0.5)

rg4courierzipper = rg3courierzipper*FockMap(rg4[:courierzipper]|>CrystalFockMap)
rg4filledzipper = rg3courierzipper*FockMap(rg4[:filledzipper]|>CrystalFockMap)
rg4emptyzipper = rg3courierzipper*FockMap(rg4[:emptyzipper]|>CrystalFockMap)
region = buildregion(crystal, 12, 12)
c3r = c3|>recenter(triangular*(3,3))
c6r = c6|>recenter(triangular*(3,3))
region = region + c3r*region + c3r*c3r*region
region = region + c6r*region
region|>visualize
regionfock = getregionfock(rg4courierzipper|>getoutspace, region)
lt = fourier(rg4courierzipper|>getoutspace, regionfock)
rt = fourier(rg4courierzipper|>getinspace, rg4courierzipper|>getinspace|>unitcellfock|>RegionFock)
rg4couriercrystalstate = FockMap(rg4courierzipper) * FockMap(rt)
rg4courierlocalstate = FockMap(lt') * rg4couriercrystalstate
visualize(rg4courierlocalstate|>RegionState|>normalize, markersize=5, logscale=0.5)
region = buildregion(crystal, 8, 8)
region = region + c3*region + c3*c3*region
region|>visualize
regionfock = getregionfock(rg4filledzipper|>getoutspace, region)
lt = fourier(rg4filledzipper|>getoutspace, regionfock)
rt = fourier(rg4filledzipper|>getinspace, rg4filledzipper|>getinspace|>unitcellfock|>RegionFock)
rg4filledcrystalstate = FockMap(rg4filledzipper) * FockMap(rt)
rg4filledlocalstate = FockMap(lt') * rg4filledcrystalstate
visualize(rg4filledlocalstate|>RegionState|>normalize, markersize=5, logscale=0.5)
lt = fourier(rg4emptyzipper|>getoutspace, regionfock)
rt = fourier(rg4emptyzipper|>getinspace, rg4emptyzipper|>getinspace|>unitcellfock|>RegionFock)
rg4emptycrystalstate = FockMap(rg4emptyzipper) * FockMap(rt)
rg4emptylocalstate = FockMap(lt') * rg4emptycrystalstate
visualize(rg4emptylocalstate|>RegionState|>normalize, markersize=5, logscale=0.5)

fiodir("/Users/alphaharrius/ZERData/chern192x192json/rglocalstates")
fiosave(rg1courierlocalstate, name="rg1courierlocalstate")
fiosave(rg1filledlocalstate, name="rg1filledlocalstate")
fiosave(rg1emptylocalstate, name="rg1emptylocalstate")
fiosave(rg2courierlocalstate, name="rg2courierlocalstate")
fiosave(rg2filledlocalstate, name="rg2filledlocalstate")
fiosave(rg2emptylocalstate, name="rg2emptylocalstate")
fiosave(rg3courierlocalstate, name="rg3courierlocalstate")
fiosave(rg3filledlocalstate, name="rg3filledlocalstate")
fiosave(rg3emptylocalstate, name="rg3emptylocalstate")
fiosave(rg4courierlocalstate, name="rg4courierlocalstate")
fiosave(rg4filledlocalstate, name="rg4filledlocalstate")
fiosave(rg4emptylocalstate, name="rg4emptylocalstate")

rg5courierzipper = rg4courierzipper*FockMap(rg5[:courierzipper]|>CrystalFockMap)
rg5filledzipper = rg4courierzipper*FockMap(rg5[:filledzipper]|>CrystalFockMap)
rg5emptyzipper = rg4courierzipper*FockMap(rg5[:emptyzipper]|>CrystalFockMap)
region = buildregion(crystal, 15, 15)
c3r = c3|>recenter(triangular*(3,3))
c6r = c6|>recenter(triangular*(3,3))
region = region + c3r*region + c3r*c3r*region
region = region + c6r*region
region|>visualize
regionfock = getregionfock(rg5courierzipper|>getoutspace, region)
lt = fourier(rg5courierzipper|>getoutspace, regionfock)
rt = fourier(rg5courierzipper|>getinspace, rg5courierzipper|>getinspace|>unitcellfock|>RegionFock)
rg5couriercrystalstate = FockMap(rg5courierzipper) * FockMap(rt)
rg5courierlocalstate = FockMap(lt') * rg5couriercrystalstate
visualize(rg5courierlocalstate|>RegionState|>normalize, markersize=5, logscale=0.5)

using JSON
function Zipper.fiosave(object; name::String)
    storageobject = object
    type = typeof(object)
    # We would like to see if there are any storage type registered for the parametric type first, 
    # for example NormalFock{Region} vs NormalFock, we will flavor the former first.
    if haskey(Zipper.STORAGE_TYPES, type)
        storageobject = convert(Zipper.STORAGE_TYPES[type.name.wrapper], object)
    elseif haskey(Zipper.STORAGE_TYPES, type.name.wrapper)
        storageobject = convert(Zipper.STORAGE_TYPES[type], object)
    end
    filepath = joinpath(fiodir(), "$name.dat")
    # Set the target name for the current thread.
    Zipper.FIO_STATE.threadtargetname[Threads.threadid()] = name
    watchprogress(desc="fiosave lowering ($name)")
    jsonstring = JSON.json(storageobject)
    unwatchprogress()
    # Reset the target name for the current thread.
    Zipper.FIO_STATE.threadtargetname[Threads.threadid()] = ""
    # lzwcompressed::Vector = lzwcompress(jsonstring)
    open(filepath, "w") do io
        write(io, jsonstring)
        @info "Saved $type object to $filepath"
    end
    return filepath
end

fiodir("/Users/alphaharrius/ZERData/chern192x192json")
fiosave(energyspectrum|>crystalfockmap, name="hamiltonian")
fiosave(correlations, name="correlations")

function _momentumoccupations(correlations)
    kcorrs = correlations|>crystalsubmaps|>Dict
    crystal::Crystal = correlations|>getoutspace|>getcrystal
    center::Offset = crystal|>getspace|>getorigin
    tracecrystal::Crystal = Crystal(center|>Subset, crystal|>getbc)
    mode::Mode = Mode(:b=>center)

    function tracing(k::Momentum, corr::FockMap)
        space::FockSpace = mode|>setattr(:k=>k)|>FockSpace
        return (k, k)=>FockMap(space, space, [corr|>tr][:, :])/dimension(corr|>getoutspace)
    end

    blocks = paralleltasks(
        name="momentumoccupations",
        tasks=(()->tracing(k, corr) for (k, corr) in kcorrs),
        count=length(kcorrs))|>parallel|>Dict
    
    return crystalfockmap(tracecrystal, tracecrystal, blocks)
end

rg1filledzipper = rg1[:block]'*rg1[:wannierfilledisometry]
rg1emptyzipper = rg1[:block]'*rg1[:wannieremptyisometry]
rg1courierzipper = rg1[:block]'*rg1[:wanniercourierisometry]
mo = _momentumoccupations(rg1filledzipper * rg1filledzipper') + _momentumoccupations(rg1emptyzipper * rg1emptyzipper')

mo|>crystalspectrum|>getbands|>first|>visualize

import CairoMakie
CairoMakie.save("/Users/alphaharrius/ZERData/momentumoccupationsRG1.svg", p)

rg2filledzipper = rg2[:block]'*rg2[:wannierfilledisometry]
rg2emptyzipper = rg2[:block]'*rg2[:wannieremptyisometry]
rg2courierzipper = rg2[:block]'*rg2[:wanniercourierisometry]
fmo = _momentumoccupations(rg1courierzipper * rg2filledzipper * rg2filledzipper' * rg1courierzipper')
emo = _momentumoccupations(rg1courierzipper * rg2emptyzipper * rg2emptyzipper' * rg1courierzipper')

p = fmo+emo|>crystalspectrum|>getbands|>first|>visualize
CairoMakie.save("/Users/alphaharrius/ZERData/momentumoccupationsRG2.svg", p)

rg3filledzipper = rg3[:block]'*rg3[:wannierfilledisometry]
rg3emptyzipper = rg3[:block]'*rg3[:wannieremptyisometry]
rg3courierzipper = rg3[:block]'*rg3[:wanniercourierisometry]
fmo = _momentumoccupations(rg1courierzipper * rg2courierzipper * rg3filledzipper * rg3filledzipper' * rg2courierzipper' * rg1courierzipper')
emo = _momentumoccupations(rg1courierzipper * rg2courierzipper * rg3emptyzipper * rg3emptyzipper' * rg2courierzipper' * rg1courierzipper')

p = fmo+emo|>crystalspectrum|>getbands|>first|>visualize
CairoMakie.save("/Users/alphaharrius/ZERData/momentumoccupationsRG3.svg", p)

rg4filledzipper = rg4[:block]'*rg4[:wannierfilledisometry]
rg4emptyzipper = rg4[:block]'*rg4[:wannieremptyisometry]
rg4courierzipper = rg4[:block]'*rg4[:wanniercourierisometry]

rg4fi = rg1courierzipper * rg2courierzipper * rg3courierzipper * rg4filledzipper
rg4fi*rg4fi'
fmo = _momentumoccupations(rg1courierzipper * rg2courierzipper * rg3courierzipper * rg4filledzipper * rg4filledzipper' * rg3courierzipper' * rg2courierzipper' * rg1courierzipper')
emo = _momentumoccupations(rg1courierzipper * rg2courierzipper * rg3courierzipper * rg4emptyzipper * rg4emptyzipper' * rg3courierzipper' * rg2courierzipper' * rg1courierzipper')

p = fmo+emo|>crystalspectrum|>getbands|>first|>visualize
CairoMakie.save("/Users/alphaharrius/ZERData/momentumoccupationsRG4.svg", p)

rg5filledzipper = rg5[:block]'*rg5[:wannierfilledisometry]
rg5emptyzipper = rg5[:block]'*rg5[:wannieremptyisometry]
rg5courierzipper = rg5[:block]'*rg5[:wanniercourierisometry]
rg5fi = rg1courierzipper * rg2courierzipper * rg3courierzipper * rg4courierzipper * rg5filledzipper
rg5fi = rg5fi|>CrystalFockMap|>FockMap
rg5fc = rg5fi*rg5fi'
blocks = Dict((k, k)=>rg5fc[subspace, subspace] for (k, subspace) in rg5fc|>getoutspace|>crystalsubspaces)
rg5fcd = CrystalFockMap(rg5fc|>getoutspace|>getcrystal, rg5fc|>getinspace|>getcrystal, blocks)
fmo = _momentumoccupations(rg5fcd)
rg5ei = rg1courierzipper * rg2courierzipper * rg3courierzipper * rg4courierzipper * rg5emptyzipper
rg5ei = rg5ei|>CrystalFockMap|>FockMap
rg5ec = rg5ei*rg5ei'
blocks = Dict((k, k)=>rg5ec[subspace, subspace] for (k, subspace) in rg5ec|>getoutspace|>crystalsubspaces)
rg5ecd = CrystalFockMap(rg5ec|>getoutspace|>getcrystal, rg5ec|>getinspace|>getcrystal, blocks)
emo = _momentumoccupations(rg5ecd)
corei = rg1courierzipper * rg2courierzipper * rg3courierzipper * rg4courierzipper * rg5courierzipper
corei = corei|>CrystalFockMap|>FockMap
corec = corei*corei'
blocks = Dict((k, k)=>corec[subspace, subspace] for (k, subspace) in corec|>getoutspace|>crystalsubspaces)
corecd = CrystalFockMap(corec|>getoutspace|>getcrystal, corec|>getinspace|>getcrystal, blocks)
cmo = _momentumoccupations(corecd)

p = fmo+emo|>crystalspectrum|>getbands|>first|>visualize
import CairoMakie
CairoMakie.save("/Users/alphaharrius/ZERData/momentumoccupationsRG5.svg", p)
p = cmo|>crystalspectrum|>getbands|>first|>visualize
CairoMakie.save("/Users/alphaharrius/ZERData/momentumoccupationsCore.svg", p)

p = mo|>crystalspectrum|>getbands|>first|>visualize
import CairoMakie
CairoMakie.save("/Users/alphaharrius/ZERData/momentumoccupationsRG4.svg", p)

function savedata(rgdata, rgname)
    fiodir("/Users/alphaharrius/ZERData/chern192x192json/$rgname")
    fiosave(rgdata[:block], name="block")
    fiosave(rgdata[:couriercorrelations], name="couriercorrelations")
    fiosave(rgdata[:filledcorrelations], name="filledcorrelations")
    fiosave(rgdata[:emptycorrelations], name="emptycorrelations")
    fiosave(rgdata[:wannierfilledisometry], name="wannierfilledisometry")
    fiosave(rgdata[:wannieremptyisometry], name="wannieremptyisometry")
    fiosave(rgdata[:wanniercourierisometry], name="wanniercourierisometry")
    fiosave(rgdata[:rawcouriercorrelations], name="rawcouriercorrelations")
    fiosave(rgdata[:wannierfilledstates]|>FockMap, name="wannierfilledstates")
    fiosave(rgdata[:wannieremptystates]|>FockMap, name="wannieremptystates")
    fiosave(rgdata[:wanniercourierstates]|>FockMap, name="wanniercourierstates")
end

savedata(rg1, "rg1")
savedata(rg2, "rg2")
savedata(rg3, "rg3")
savedata(rg4, "rg4")
savedata(rg5, "rg5")
savedata(rg6, "rg6")
savedata(rg7, "rg7")

visualize(rg6[:wanniercourierstates]|>normalize, markersize=5, logscale=0.5)
