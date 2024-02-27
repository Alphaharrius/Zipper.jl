using LinearAlgebra
using Zipper

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

m = quantize(unitcell, 1)|>first

tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m, m|>setattr(:r=>[1, 0]∈square))=>tₙ,
    (m, m|>setattr(:r=>[0, 1]∈square))=>tₙ])

@info "Computing energy spectrum..."
energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
@info "Resolving ground states..."
groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=0.5)
@info "Computing ground state correlations..."
groundstateprojector = groundstates |> crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector
@info "Performing initial blocking..."
initscale = Scale([2 0; 0 2], square)
initblock = initscale * (correlations|>getoutspace)
inputcorrelations = initblock * correlations * initblock'

function zer(correlations)
    @info "Performing blocking..."
    scale = Scale([2 0; 0 2], square)
    block = scale * (correlations|>getoutspace)
    blockedcorrelations = block * correlations * block'
    blockedcrystalfock = blockedcorrelations|>getoutspace
    blockedcrystal = blockedcrystalfock|>getcrystal

    @info "Performing extended restrictions..."
    normalvector = (1, 1) ∈ (blockedcrystal|>getspace)
    extendedrestrict = extendedcrystalrestrict(
        crystal=blockedcrystal, normalvector=normalvector, stripradius=0.5)
    restrict = extendedrestrict * blockedcrystalfock
    stripcorrelations = restrict * blockedcorrelations * restrict'
    stripspectrum = stripcorrelations|>crystalspectrum

    @info "Computing strip frozen correlations..."
    stripfrozenstates = distillation(stripspectrum, :frozen=>(v -> v < 0.003 || v > 0.99))[:frozen]
    stripfrozenprojector = stripfrozenstates|>crystalprojector
    stripfrozencorrelations = idmap(stripfrozenprojector|>getoutspace) - stripfrozenprojector

    @info "Computing strip truncation restricted Fourier transform..."
    stripunitcell = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*3, radius=0.5)
    striphomefock = stripfrozencorrelations|>getoutspace|>unitcellfock
    basistohomemodes = ((mode|>getattr(:b)|>basispoint)=>mode for mode in striphomefock)
    conversionmappings = Dict(
        mode=>(mode|>setattr(:r=>getattr(mode, :b)-b)|>setattr(:b=>b)) for (b, mode) in basistohomemodes)
    actualstriphomefock = conversionmappings|>values|>RegionFock
    truncationregionfock = RegionFock(
        mode|>setattr(:r=>getattr(mode, :r)+normalvector*n) for mode in actualstriphomefock for n in 0:2)
    homemappings = Dict(tomode=>frommode|>removeattr(:r) for (tomode, frommode) in conversionmappings)
    truncate = (
        fourier(stripfrozencorrelations|>getoutspace, truncationregionfock, homemappings) / sqrt(stripfrozencorrelations|>getoutspace|>getcrystal|>vol))

    @info "Performing truncation..."
    truncatedregioncorrelations = truncate' * stripfrozencorrelations * truncate
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
    truncatedstripfrozencorrelations = crystaldirectsum(
        (transform * prunedcorrelations * transform' for transform in stripfouriers), 
        outcrystal=stripcrystal, incrystal=stripcrystal)
    truncatedstripfrozencorrelationspectrum = truncatedstripfrozencorrelations|>crystalspectrum

    @info "Extracting strip filled states..."
    stripfilledstates = groundstatespectrum(truncatedstripfrozencorrelationspectrum, perunitcellfillings=8)
    stripfilledprojector = stripfilledstates|>crystalprojector
    stripfilledcorrelations = idmap(stripfilledprojector|>getoutspace) - stripfilledprojector

    @info "Searching strip filled seeds..."
    scaledcrystal = stripfilledstates|>getcrystal
    scaledspace = scaledcrystal|>getspace

    wannierregion1 = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.875, radius=0.5, minbottomheight=0.15)
    regionfock = quantize(wannierregion1, 1)
    remapper = spatialremapper(regionfock, offsets=scaledspace|>getorigin|>Subset, unitcell=scaledcrystal|>getunitcell)
    wannierregionfock = remapper|>getoutspace
    restrict = fourier(stripfilledcorrelations|>getoutspace, wannierregionfock) / (scaledcrystal|>vol|>sqrt)
    localcorrelations = restrict' * stripfilledcorrelations * restrict
    localcorrelations = remapper' * localcorrelations * remapper
    region1localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3])[1]
    region1localstates = m135 * region1localstates
    region1localstates = remapper*FockMap(region1localstates)|>RegionState

    wannierregion2 = wannierregion1 .- normalvector*0.5
    regionfock = quantize(wannierregion2, 1)
    offsets = Subset(-(scaledspace*normalvector), scaledspace|>getorigin)
    remapper = spatialremapper(regionfock, offsets=offsets, unitcell=scaledcrystal|>getunitcell)
    wannierregionfock = remapper|>getoutspace
    restrict = fourier(stripfilledcorrelations|>getoutspace, wannierregionfock) / (scaledcrystal|>vol|>sqrt)
    localcorrelations = restrict' * stripfilledcorrelations * restrict
    localcorrelations = remapper' * localcorrelations * remapper
    region2localstates = getregionstates(localcorrelations=localcorrelations, grouping=[2, 5])[2]
    region2localstates = m135 * region2localstates
    region2localstates = remapper*FockMap(region2localstates)|>RegionState

    localseeds = region1localstates+region2localstates|>FockMap
    transform = fourier(stripfilledcorrelations|>getoutspace, localseeds|>getoutspace) / (scaledcrystal|>vol|>sqrt)
    crystalseeds = transform * localseeds
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

    @info "Computing strip filled projector..."
    remapper = spatialremapper(
        wannierlocalisometry|>getoutspace, getsphericalregion(crystal=blockedcrystal, radius=2, metricspace=blockedcrystal|>getspace))
    wannierlocalisometry = remapper * wannierlocalisometry
    transform = fourier(blockedcrystalfock, wannierlocalisometry|>getoutspace)
    wanniercrystalisometry = transform * wannierlocalisometry
    wanniermetalprojector = wanniercrystalisometry .* wanniercrystalisometry'

    @info "Computing global distiller..."
    transform = c4 * blockedcrystalfock
    globaldistiller = wanniermetalprojector + transform * wanniermetalprojector * transform'
    globaldistillerspectrum = globaldistiller|>crystalspectrum

    @info "Distilling courier bands..."
    courierstates = groundstatespectrum(globaldistillerspectrum, perunitcellfillings=4)
    courierprojector = courierstates|>crystalprojector
    couriercorrelations = idmap(courierprojector|>getoutspace) - courierprojector

    @info "Searching courier seeds..."
    wannierregion = getsphericalregion(crystal=blockedcrystal, radius=0.25, metricspace=blockedcrystal|>getspace)
    wannierregion = wannierregion .+ ([0.25, 0.25]∈getspace(wannierregion))
    regionfock = quantize(wannierregion, 1)
    restrict = fourier(blockedcrystalfock, regionfock) / (blockedcrystal|>vol|>sqrt)
    localcorrelations = restrict' * couriercorrelations * restrict
    regionseeds = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]

    regionseedmap = regionseeds|>FockMap
    for element in pointgroupelements(c4|>recenter([0.5, 0.5]∈getspace(wannierregion)))
        lefttransform = element * RegionFock(regionseedmap|>getoutspace)
        righttransform = element * RegionFock(regionseedmap|>getinspace)
        regionseedmap += lefttransform * regionseedmap * righttransform'
    end
    courierseeds = RegionFock(regionseedmap|>getoutspace) * regionseedmap * RegionFock(regionseedmap|>getinspace)

    transform = fourier(blockedcrystalfock, courierseeds|>getoutspace) / (blockedcrystal|>vol|>sqrt)
    crystalcourierseeds = transform * courierseeds
    crystalcourierseeds = Dict(k=>crystalcourierseeds[k, :] for k in blockedcrystal|>brillouinzone)
    pseudoidens = (u' * u for (_, u) in crystalcourierseeds)
    lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
    @info "Courier seeds linear-dependence metric: $lineardepmetric"

    @info "Wannierizing courier bands..."
    wanniercourierisometry = wannierprojection(
        crystalisometries=courierstates|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds)
    visualregion = getsphericalregion(crystal=blockedcrystal, radius=2, metricspace=blockedcrystal|>getspace)
    visualfock = quantize(visualregion, 1)
    leftrestrict = fourier(blockedcrystalfock, visualfock) / (blockedcrystal|>vol|>sqrt)
    rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
    courierlocalstates = leftrestrict' * wanniercourierisometry * rightrestrict

    @info "Computing courier correlations..."
    couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
    courierspectrum = couriercorrelations|>crystalspectrum
    @info "Purifying courier states..."
    purifiedspectrum = roundingpurification(courierspectrum)
    purifiedcorrelations = purifiedspectrum|>CrystalFockMap

    return Dict(
        :couriercorrelations=>couriercorrelations,
        :courierisometry=>wanniercourierisometry,
        :correlations=>purifiedcorrelations,
        :courierlocalstates=>courierlocalstates|>RegionState,
        :stripfilledlocalstates=>stripfilledlocalstates,
        :globaldistiller=>globaldistiller,
        :block=>block)
end

rg1 = zer(inputcorrelations)
rg2 = zer(rg1[:correlations])
rg3 = zer(rg2[:correlations])
rg4 = zer(rg3[:correlations])

fiodir("/Users/alphaharrius/ZERData/squaremetal")
fiosave(energyspectrum|>CrystalFockMap, name="H")
fiosave(correlations, name="C")
fiosave(initblock, name="initblock")
fiosave(inputcorrelations, name="inputC")

fiodir("/Users/alphaharrius/ZERData/squaremetal/RG1")
fiosave(rg1[:couriercorrelations], name="rg1couriercorrelations")
fiosave(rg1[:courierisometry], name="rg1courierisometry")
fiosave(rg1[:correlations], name="rg1correlations")
fiosave(rg1[:globaldistiller], name="rg1globaldistiller")
fiosave(rg1[:block], name="rg1block")

fiodir("/Users/alphaharrius/ZERData/squaremetal/RG2")
fiosave(rg2[:couriercorrelations], name="rg2couriercorrelations")
fiosave(rg2[:courierisometry], name="rg2courierisometry")
fiosave(rg2[:correlations], name="rg2correlations")
fiosave(rg2[:globaldistiller], name="rg2globaldistiller")
fiosave(rg2[:block], name="rg2block")

fiodir("/Users/alphaharrius/ZERData/squaremetal/RG3")
fiosave(rg3[:couriercorrelations], name="rg3couriercorrelations")
fiosave(rg3[:courierisometry], name="rg3courierisometry")
fiosave(rg3[:correlations], name="rg3correlations")
fiosave(rg3[:globaldistiller], name="rg3globaldistiller")
fiosave(rg3[:block], name="rg3block")

fiodir("/Users/alphaharrius/ZERData/squaremetal/RG4")
fiosave(rg4[:couriercorrelations], name="rg4couriercorrelations")
fiosave(rg4[:courierisometry], name="rg4courierisometry")
fiosave(rg4[:correlations], name="rg4correlations")
fiosave(rg4[:globaldistiller], name="rg4globaldistiller")
fiosave(rg4[:block], name="rg4block")
