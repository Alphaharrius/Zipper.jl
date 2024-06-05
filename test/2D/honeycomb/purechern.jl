using LinearAlgebra
using Zipper, Plots
plotlyjs()

x = -1:0.01:1
f = 0.5 .* (1 .+ x.^2)
plot(x, f)
plot!(legend=false)

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

fiodir("/Users/alphaharrius/ZERData/chern192x192")
hamiltonian = fioload("hamiltonian")
correlations = fioload("correlations")

function purezer(correlations)
    @info("Starting RG...")
    crystalfock = correlations|>getoutspace
    sortedhomefock = NormalFock(sort(crystalfock|>unitcellfock|>collect, by=mode->mode|>getattr(:b)|>vec|>first))
    crystalfock = CrystalFock(crystalfock|>getcrystal, crystalfock.korderings, sortedhomefock)

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
    distilled = groupbands(globaldistiller|>crystalspectrum, :filled=>(v->v < -1e-3), :empty=>(v->v > 1e-3))

    blockedspectrum = blockedcorrelations|>crystalspectrum
    blocked = groupbands(blockedspectrum, :filled=>(v->v < 0.5), :empty=>(v->v > 0.5))

    filledbands = distilled[:filled]
    filledisometries = filledbands|>geteigenvectors

    emptybands = distilled[:empty]
    emptyisometries = emptybands|>geteigenvectors

    function svdpurification(k, ψ, Φ)
        ϕ = Φ[k]
        p = ϕ'ψ
        u, s, vd = p|>svd
        @assert all([v for (_, v) in s].>0.9)
        m = u*vd
        return k=>ϕ*m
    end

    purifiedfilled = paralleltasks(
        name="svdpurification",
        tasks=(()->svdpurification(k, ψ, blocked[:filled]|>geteigenvectors) for (k, ψ) in filledisometries),
        count=length(filledisometries))|>parallel|>Dict

    purifiedempty = paralleltasks(
        name="svdpurification",
        tasks=(()->svdpurification(k, ψ, blocked[:empty]|>geteigenvectors) for (k, ψ) in emptyisometries),
        count=length(emptyisometries))|>parallel|>Dict

    purefilledbands = CrystalSpectrum{2}(filledbands|>getcrystal, filledbands|>geteigenmodes, filledbands|>geteigenvalues, purifiedfilled, filledbands.bandcount)
    pureemptybands = CrystalSpectrum{2}(emptybands|>getcrystal, emptybands|>geteigenmodes, emptybands|>geteigenvalues, purifiedempty, emptybands.bandcount)

    puredistiller = (pureemptybands|>crystalprojector) - (purefilledbands|>crystalprojector)

    purebands = groupbands(puredistiller|>crystalspectrum, :filled=>(v -> v < -1e-3), :empty=>(v -> v > 1e-3))

    @info "Searching seed for filled bands..."
    filledprojector = crystalprojector(filledbands)
    filledcorrelations = idmap(filledprojector|>getoutspace) - filledprojector
    filledlocalcorrelations = localrestrict' * filledcorrelations * localrestrict
    filledseeds = getregionstates(localcorrelations=filledlocalcorrelations, grouping=[3])[1]
    filledseeds = c6 * filledseeds

    @info "Searching seed for empty bands..."
    emptyprojector = crystalprojector(emptybands)
    emptycorrelations = idmap(emptyprojector|>getoutspace) - emptyprojector
    emptylocalcorrelations = localrestrict' * emptycorrelations * localrestrict
    emptyseeds = getregionstates(localcorrelations=emptylocalcorrelations, grouping=[3])[1]
    emptyseeds = c6 * emptyseeds

    @info "Wannierizing filled bands..."
    crystalfilledrseeds = crystalisometries(localisometry=filledseeds|>FockMap, crystalfock=blockedcrystalfock, addinspacemomentuminfo=true)
    wannierfilledisometry = wannierprojection(
        crystalisometries=purebands[:filled]|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledrseeds|>Dict)

    @info "Wannierizing empty bands..."
    crystalemptyrseeds = crystalisometries(localisometry=emptyseeds|>FockMap, crystalfock=blockedcrystalfock, addinspacemomentuminfo=true)
    wannieremptyisometry = wannierprojection(
        crystalisometries=purebands[:empty]|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyrseeds|>Dict)

    visualregion = getsphericalregion(crystal=blockedcrystal, radius=3, metricspace=blockedspace|>orthospace)
    visualfock = quantize(visualregion, 1)

    @info "Computing local filled states..."
    leftrestrict = fourier(wannierfilledisometry|>getoutspace, visualfock) / (blockedcrystal|>vol|>sqrt)
    rightrestrict = fourier(wannierfilledisometry|>getinspace, wannierfilledisometry|>getinspace|>unitcellfock|>RegionFock)
    wannierfilledstates = leftrestrict' * wannierfilledisometry * rightrestrict

    @info "Computing local empty states..."
    leftrestrict = fourier(wannieremptyisometry|>getoutspace, visualfock) / (blockedcrystal|>vol|>sqrt)
    rightrestrict = fourier(wannieremptyisometry|>getinspace, wannieremptyisometry|>getinspace|>unitcellfock|>RegionFock)
    wannieremptystates = leftrestrict' * wannieremptyisometry * rightrestrict

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

    @info "Wannierizing courier bands..."
    courierseeds = siteAseed + siteBseed
    crystalcourierseeds = crystalisometries(localisometry=courierseeds|>FockMap, crystalfock=blockedcrystalfock, addinspacemomentuminfo=true)
    wanniercourierisometry = wannierprojection(
        crystalisometries=purebands[:others]|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds|>Dict)

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

    return Dict(
        :filledstates=>wannierfilledstates,
        :emptystates=>wannieremptystates,
        :courierstate=>wanniercourierstate,
        :filledcorrelations=>filledcorrelations,
        :emptycorrelations=>emptycorrelations,
        :couriercorrelations=>couriercorrelations,
        :initialdistiller=>globaldistiller,
        :puredistiller=>puredistiller,
        :wannierfilledisometry=>wannierfilledisometry,
        :wannieremptyisometry=>wannieremptyisometry,
        :wanniercourierisometry=>wanniercourierisometry)
end

rg1 = purezer(correlations)
rg2 = purezer(rg1[:couriercorrelations])
rg3 = purezer(rg2[:couriercorrelations])
rg4 = purezer(rg3[:couriercorrelations])
rg5 = purezer(rg4[:couriercorrelations])
rg6 = purezer(rg5[:couriercorrelations])

rg6region = buildregion(rg6[:couriercorrelations]|>getoutspace|>getcrystal, 3, 3)
rg6regionfock = getregionfock(rg6[:couriercorrelations]|>getoutspace, rg6region)
restrict = fourier(rg6[:couriercorrelations]|>getoutspace, rg6regionfock) / (rg6[:couriercorrelations]|>getoutspace|>getcrystal|>vol|>sqrt)
rg6localcorrelations = restrict' * rg6[:couriercorrelations] * restrict
rg6localcorrelations|>visualize
rg6localcorrelations|>eigspech|>visualize

rg6[:couriercorrelations]|>crystalspectrum|>visualize
visualize(rg3[:courierstate]|>RegionState|>normalize, markersize=5, logscale=0.5)
visualize(rg3[:filledstates]|>RegionState|>normalize, markersize=5, logscale=0.5)
visualize(rg3[:emptystates]|>RegionState|>normalize, markersize=5, logscale=0.5)

fiodir("/Users/alphaharrius/ZERData/chern192x192/adriandata")
fiosave(rg6localcorrelations, name="corelocalcorrelations")
fiodir("/Users/alphaharrius/ZERData/chern192x192/adriandata/RG1")
fiosave(rg1[:wanniercourierisometry]|>CrystalFockMap, name="wanniercourierisometry")
fiosave(rg1[:wannierfilledisometry]|>CrystalFockMap, name="wannierfilledisometry")
fiosave(rg1[:wannieremptyisometry]|>CrystalFockMap, name="wannieremptyisometry")
fiosave(rg1[:courierstate], name="courierstate")
fiosave(rg1[:filledstates], name="filledstates")
fiosave(rg1[:emptystates], name="emptystates")
fiodir("/Users/alphaharrius/ZERData/chern192x192/adriandata/RG2")
fiosave(rg2[:wanniercourierisometry]|>CrystalFockMap, name="wanniercourierisometry")
fiosave(rg2[:wannierfilledisometry]|>CrystalFockMap, name="wannierfilledisometry")
fiosave(rg2[:wannieremptyisometry]|>CrystalFockMap, name="wannieremptyisometry")
fiosave(rg2[:courierstate], name="courierstate")
fiosave(rg2[:filledstates], name="filledstates")
fiosave(rg2[:emptystates], name="emptystates")
fiodir("/Users/alphaharrius/ZERData/chern192x192/adriandata/RG3")
fiosave(rg3[:wanniercourierisometry]|>CrystalFockMap, name="wanniercourierisometry")
fiosave(rg3[:wannierfilledisometry]|>CrystalFockMap, name="wannierfilledisometry")
fiosave(rg3[:wannieremptyisometry]|>CrystalFockMap, name="wannieremptyisometry")
fiosave(rg3[:courierstate], name="courierstate")
fiosave(rg3[:filledstates], name="filledstates")
fiosave(rg3[:emptystates], name="emptystates")
fiodir("/Users/alphaharrius/ZERData/chern192x192/adriandata/RG4")
fiosave(rg4[:wanniercourierisometry]|>CrystalFockMap, name="wanniercourierisometry")
fiosave(rg4[:wannierfilledisometry]|>CrystalFockMap, name="wannierfilledisometry")
fiosave(rg4[:wannieremptyisometry]|>CrystalFockMap, name="wannieremptyisometry")
fiosave(rg4[:courierstate], name="courierstate")
fiosave(rg4[:filledstates], name="filledstates")
fiosave(rg4[:emptystates], name="emptystates")
fiodir("/Users/alphaharrius/ZERData/chern192x192/adriandata/RG5")
fiosave(rg5[:wanniercourierisometry]|>CrystalFockMap, name="wanniercourierisometry")
fiosave(rg5[:wannierfilledisometry]|>CrystalFockMap, name="wannierfilledisometry")
fiosave(rg5[:wannieremptyisometry]|>CrystalFockMap, name="wannieremptyisometry")
fiosave(rg5[:courierstate], name="courierstate")
fiosave(rg5[:filledstates], name="filledstates")
fiosave(rg5[:emptystates], name="emptystates")
fiodir("/Users/alphaharrius/ZERData/chern192x192/adriandata/RG6")
fiosave(rg6[:wanniercourierisometry]|>CrystalFockMap, name="wanniercourierisometry")
fiosave(rg6[:wannierfilledisometry]|>CrystalFockMap, name="wannierfilledisometry")
fiosave(rg6[:wannieremptyisometry]|>CrystalFockMap, name="wannieremptyisometry")
fiosave(rg6[:courierstate], name="courierstate")
fiosave(rg6[:filledstates], name="filledstates")
fiosave(rg6[:emptystates], name="emptystates")
