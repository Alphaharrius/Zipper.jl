using Revise, LinearAlgebra, AbstractAlgebra, SmithNormalForm
using Zipper

setmaxthreads(Threads.nthreads())

square = euclidean(RealSpace, 2)
point = [1/2, 1/2] ∈ square
spatialsnappingcalibration([point])

kspace = convert(MomentumSpace, square)

c4 = pointgrouptransform([0 -1; 1 0])

unitcell = Subset(point)
crystal = Crystal(unitcell, [128, 128])
reciprocalhashcalibration(crystal.sizes)

m = quantize(unitcell, 1)|>first

tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m, m |> setattr(:r => [1, 0] ∈ square)) => tₙ,
    (m, m |> setattr(:r => [0, 1] ∈ square)) => tₙ])

@info("Computing energy spectrum...")
energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)

@info("Computing ground state correlations...")
@time begin
    groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=0.5)
    groundstateprojector = groundstates |> crystalprojector
    correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector
end

crystalfock = correlations|>getoutspace

scale = Scale([4 0; 0 4], square)

@info("Performing blocking...")
@time begin
    blocker = scale * (correlations|>getoutspace)
    blockedcorrelations = blocker * correlations * blocker'
    blockedcrystalfock = blockedcorrelations|>getoutspace
    blockedcrystal = blockedcrystalfock|>getcrystal
end

normalvector = (1, 1) ∈ (blockedcrystal|>getspace)

extendedrestrict = extendedcrystalrestrict(crystal=blockedcrystal, normalvector=normalvector, stripradius=0.5)

@info("Performing strip restriction...")
@time begin
    restriction = extendedrestrict * blockedcrystalfock
    stripcorrelations = restriction * blockedcorrelations * restriction'
    stripcorrelationspectrum = stripcorrelations|>crystalspectrum
end
stripcorrelationspectrum|>linespectrum|>visualize

@info("Computing strip frozen correlations...")
@time begin
    stripfrozenstates = distillation(stripcorrelationspectrum, :frozen => v -> v < 0.003 || v > 0.99)[:frozen]
    stripfrozenprojector = stripfrozenstates|>crystalprojector
    stripfrozencorrelations = idmap(stripfrozenprojector|>getoutspace) - stripfrozenprojector
    stripfrozencorrelationspectrum = stripfrozencorrelations|>crystalspectrum
end
stripfrozencorrelationspectrum|>linespectrum|>visualize

@info("Computing strip truncation restricted Fourier transform...")
@time begin
    stripunitcell = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*5, radius=0.5)
    striphomefock = stripfrozencorrelations|>getoutspace|>unitcellfock
    basistohomemodes = ((mode|>getattr(:b)|>basispoint)=>mode for mode in striphomefock)
    conversionmappings = Dict(mode=>(mode|>setattr(:r=>getattr(mode, :b)-b)|>setattr(:b=>b)) for (b, mode) in basistohomemodes)
    actualstriphomefock = conversionmappings|>values|>RegionFock
    truncationregionfock = RegionFock(mode|>setattr(:r=>getattr(mode, :r)+normalvector*n) for mode in actualstriphomefock for n in 0:2)
    homemappings = Dict(tomode=>frommode|>removeattr(:r) for (tomode, frommode) in conversionmappings)
    truncator = fourier(stripfrozencorrelations|>getoutspace, truncationregionfock, homemappings) / sqrt(stripfrozencorrelations|>getoutspace|>getcrystal|>vol)
end

@info("Performing truncation on strip frozen correlations...")
@time begin
    truncationregioncorrelations = truncator' * stripfrozencorrelations * truncator
    offsets = Subset(normalvector * n for n in 0:2)
    remapper = spatialremapper(truncationregionfock; offsets=offsets, unitcell=stripunitcell)
    truncationregioncorrelations = remapper * truncationregioncorrelations * remapper'

    truncationregionindices = Iterators.product(truncationregioncorrelations|>getoutspace, truncationregioncorrelations|>getoutspace)
    truncationregionbonds = Dict((getattr(tomode, :r) - getattr(frommode, :r), frommode|>removeattr(:r), tomode|>removeattr(:r)) => (frommode, tomode) for (frommode, tomode) in truncationregionindices)
    pruningindices = (index for (_, index) in truncationregionbonds)
    prunedcorrelations = extractindices(truncationregioncorrelations, pruningindices)
    prunedcorrelations = remapper' * prunedcorrelations * remapper
    stripfouriers = (truncator[subspace, :] for (_, subspace) in truncator|>getoutspace|>crystalsubspaces)
    stripcrystal = truncator|>getoutspace|>getcrystal
    truncatedstripfrozencorrelations = crystaldirectsum((transform * prunedcorrelations * transform' for transform in stripfouriers), outcrystal=stripcrystal, incrystal=stripcrystal)
    truncatedstripfrozencorrelationspectrum = truncatedstripfrozencorrelations|>crystalspectrum
end
truncatedstripfrozencorrelationspectrum|>linespectrum|>visualize

@info("Retrieving strip quasi-metallic states...")
@time begin
    quasistripmetalstate = groundstatespectrum(truncatedstripfrozencorrelationspectrum, perunitcellfillings=8)
    quasistripmetalprojector = quasistripmetalstate|>crystalprojector
    quasistripmetalcorrelations = idmap(quasistripmetalprojector|>getoutspace) - quasistripmetalprojector
end
quasistripmetalstate|>linespectrum|>visualize

c2 = pointgrouptransform([-1 0; 0 -1], localspace=square)
m135 = pointgrouptransform([0 -1; -1 0], localspace=square)
m45 = pointgrouptransform([0 1; 1 0], localspace=square)

wannierregion = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.875, radius=0.5, minbottomheight=0.15)
wannierregionfock = quantize(wannierregion, 1)
scaledcrystal = quasistripmetalstate|>getcrystal
scaledspace = scaledcrystal|>getspace
remapper = spatialremapper(wannierregionfock, offsets=scaledspace|>getorigin|>Subset, unitcell=quasistripmetalstate|>getcrystal|>getunitcell)
wannierregionfock = remapper|>getoutspace
wannierrestrict = fourier(quasistripmetalcorrelations|>getoutspace, wannierregionfock) / (scaledcrystal|>vol|>sqrt)
localcorrelations = wannierrestrict' * quasistripmetalcorrelations * wannierrestrict
localcorrelations|>eigspech|>visualize
localstatesR00 = getregionstates(localcorrelations=localcorrelations, grouping=[3, 5])|>collect

wannierregion = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.875, radius=0.5, minbottomheight=0.15) - normalvector*0.5
wannierregionfock = quantize(wannierregion, 1)
offsets = Subset(-(scaledspace*normalvector), scaledspace|>getorigin)
remapper = spatialremapper(wannierregionfock, offsets=offsets, unitcell=quasistripmetalstate|>getcrystal|>getunitcell)
wannierregionfock = remapper|>getoutspace
wannierrestrict = fourier(quasistripmetalcorrelations|>getoutspace, wannierregionfock) / (scaledcrystal|>vol|>sqrt)
localcorrelations = wannierrestrict' * quasistripmetalcorrelations * wannierrestrict
localcorrelations|>eigspech|>visualize
localstatesR05 = getregionstates(localcorrelations=localcorrelations, grouping=[2, 5])|>collect

wannierseedstate = localstatesR00[1] + localstatesR05[2]
wannierseeds = wannierseedstate|>FockMap

seedstransform = fourier(quasistripmetalcorrelations|>getoutspace, wannierseeds|>getoutspace) / (scaledcrystal|>vol|>sqrt)
crystalseeds = seedstransform * wannierseeds
pseudoidentities = (crystalseeds[subspace, :]' * crystalseeds[subspace, :] for (_, subspace) in crystalseeds|>getoutspace|>crystalsubspaces)
(v for id in pseudoidentities for (_, v) in id|>eigvalsh)|>minimum

crystalseeds = Dict(k=>crystalseeds[subspace, :] for (k, subspace) in crystalseeds|>getoutspace|>crystalsubspaces)
wanniermetalisometries = wannierprojection(crystalisometries=quasistripmetalstate|>geteigenvectors, crystal=scaledcrystal, crystalseeds=crystalseeds)
wanniermetalisometries' * stripcorrelations * wanniermetalisometries |>crystalspectrum|>linespectrum|>visualize

leftrestrict = fourier(wanniermetalisometries|>getoutspace, wannierseeds|>getoutspace) / (scaledcrystal|>vol|>sqrt)
wanniercrystal = wanniermetalisometries|>getinspace|>getcrystal
rightrestrict = fourier(wanniermetalisometries|>getinspace, wanniermetalisometries|>getinspace|>unitcellfock|>RegionFock) / (wanniercrystal|>vol|>sqrt)
wannierlocalisometries = leftrestrict' * FockMap(wanniermetalisometries) * rightrestrict
visualize(wannierlocalisometries|>RegionState, markersizemultiplier=20, markersizescaling=0.1)

remapper = spatialremapper(wannierlocalisometries|>getoutspace, getsphericalregion(crystal=blockedcrystal, radius=2, metricspace=blockedcrystal|>getspace))
wannierlocalisometries = remapper * wannierlocalisometries
wanniermetalprojector = crystalprojector(localisometry=wannierlocalisometries, crystalfock=blockedcrystalfock)
wanniermetalprojector|>crystalspectrum|>visualize

c4T = c4 * blockedcrystalfock
c4wanniermetalprojector = c4T * wanniermetalprojector * c4T'

globaldistiller = wanniermetalprojector + c4wanniermetalprojector
globaldistillerspectrum = globaldistiller|>crystalspectrum
globaldistillerspectrum|>visualize
courierspectrum = groundstatespectrum(globaldistillerspectrum, perunitcellfillings=4)
courierspectrum|>visualize
courierprojector = courierspectrum|>crystalprojector
couriercorrelations = idmap(courierprojector|>getoutspace) - courierprojector
couriercorrelationspectrum = couriercorrelations|>crystalspectrum

wannierregion = couriercorrelations|>getoutspace|>getcrystal|>getunitcell
wannierregion1 = wannierregion - ([0.25, 0.25]∈getspace(wannierregion))
wannierregion2 = wannierregion + ([0.25, 0.25]∈getspace(wannierregion))
wannierregion3 = wannierregion + ([0.25, -0.25]∈getspace(wannierregion))
wannierregion4 = wannierregion + ([-0.25, 0.25]∈getspace(wannierregion))
visualize(wannierregion, wannierregion1, wannierregion2, wannierregion3, wannierregion4)

restrictfock = quantize(wannierregion1, 1)
restrict = fourier(couriercorrelations|>getoutspace, restrictfock) / (couriercorrelations|>getoutspace|>getcrystal|>vol|>sqrt)
localcorrelations = restrict' * couriercorrelations * restrict
localstates1 = getregionstates(localcorrelations=localcorrelations, grouping=[1, 2, 4])|>collect

restrictfock = quantize(wannierregion2, 1)
restrict = fourier(couriercorrelations|>getoutspace, restrictfock) / (couriercorrelations|>getoutspace|>getcrystal|>vol|>sqrt)
localcorrelations = restrict' * couriercorrelations * restrict
localstates2 = getregionstates(localcorrelations=localcorrelations, grouping=[1, 2, 4])|>collect

restrictfock = quantize(wannierregion3, 1)
restrict = fourier(couriercorrelations|>getoutspace, restrictfock) / (couriercorrelations|>getoutspace|>getcrystal|>vol|>sqrt)
localcorrelations = restrict' * couriercorrelations * restrict
localstates3 = getregionstates(localcorrelations=localcorrelations, grouping=[1, 2, 4])|>collect

restrictfock = quantize(wannierregion4, 1)
restrict = fourier(couriercorrelations|>getoutspace, restrictfock) / (couriercorrelations|>getoutspace|>getcrystal|>vol|>sqrt)
localcorrelations = restrict' * couriercorrelations * restrict
localstates4 = getregionstates(localcorrelations=localcorrelations, grouping=[1, 2, 4])|>collect

wannierregionorigin = wannierregion
wannierregionorigin|>visualize
restrictfock = quantize(wannierregionorigin, 1)
restrict = fourier(couriercorrelations|>getoutspace, restrictfock) / (couriercorrelations|>getoutspace|>getcrystal|>vol|>sqrt)
localcorrelations = restrict' * couriercorrelations * restrict
localcorrelations|>eigspech|>visualize
localstatesorigin = getregionstates(localcorrelations=localcorrelations, grouping=[2, 2])|>collect
visualize(localstatesorigin[1], markersizemultiplier=20, markersizescaling=0.3)

seeds = localstates1[1] + localstates2[1] + localstates3[1] + localstates4[1]
visualize(seeds|>FockMap|>RegionState, markersizemultiplier=20, markersizescaling=0.3)
wannierseeds = seeds|>FockMap

seedstransform = fourier(couriercorrelations|>getoutspace, seeds|>FockMap|>getoutspace) / sqrt(couriercorrelations|>getoutspace|>getcrystal|>vol)
crystalseeds = seedstransform * wannierseeds
pseudoidentities = (crystalseeds[subspace, :]' * crystalseeds[subspace, :] for (_, subspace) in crystalseeds|>getoutspace|>crystalsubspaces)
(v for id in pseudoidentities for (_, v) in id|>eigvalsh)|>minimum

crystalseeds = Dict(k=>crystalseeds[subspace, :] for (k, subspace) in crystalseeds|>getoutspace|>crystalsubspaces)
wanniercourierisometries = wannierprojection(crystalisometries=courierspectrum|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalseeds)
visualregion = getsphericalregion(crystal=blockedcrystal, radius=2, metricspace=blockedcrystal|>getspace) + ([0.5, 0.5]∈getspace(blockedcrystal))
visualfock = quantize(visualregion, 1)
leftfourier = fourier(wanniercourierisometries|>getoutspace, visualfock) / sqrt(wanniercourierisometries|>getoutspace|>getcrystal|>vol)
rightfourier = fourier(wanniercourierisometries|>getinspace, wanniercourierisometries|>getinspace|>unitcellfock|>RegionFock)
wannierlocalstates = leftfourier' * FockMap(wanniercourierisometries) * rightfourier
visualize(wannierlocalstates|>RegionState, markersizemultiplier=20, markersizescaling=0.3)

wanniercourierisometries' * blockedcorrelations * wanniercourierisometries |>crystalspectrum|>visualize

# TODO: This multiplication is wrong using CrystalFockMap
c4transform = c4 * blockedcrystalfock |>FockMap
fullwanniermetalprojector = FockMap(wanniermetalprojector) + c4transform * FockMap(wanniermetalprojector) * c4transform'
fullwanniermetalprojector|>crystalspectrum|>visualize

courierstates = groundstatespectrum(fullwanniermetalprojector|>crystalspectrum, perunitcellfillings=4)
courierstates|>visualize

courierprojector = courierstates|>crystalprojector
courierprojector * blockedcorrelations * courierprojector' |>crystalspectrum|>visualize

# ===================================================

wannierregion = getcrosssection(crystal=blockedcrystal, normalvector=normalvector * 1.5, radius=0.5, minbottomheight=0.15)
wannierregion|>visualize
scaledspace = stripcorrelations|>getoutspace|>getcrystal|>getspace
scaledunitcell = Subset(scaledspace * r for r in stripunitcell)
wannierregionfock = quantize(Subset(scaledspace * r for r in wannierregion), 1)
remapper = spatialremapper(wannierregionfock, offsets=Subset(-(scaledspace * normalvector), scaledspace|>getorigin, scaledspace * normalvector), unitcell=scaledunitcell)
wannierregionfock = remapper|>getoutspace

using Base.Iterators
function getregionstates(; localcorrelations::FockMap{RegionFock, <:Zipper.FockSpace}, grouping::Vector{<:Integer})
    localspectrum::EigenSpectrum = localcorrelations|>eigspech
    current = 1
    states = []
    for group in grouping
        currentstate = geteigenvectors(localspectrum)[:, current:current+group-1]
        push!(states, currentstate|>RegionState)
        current += group
    end
    return states
end

@info("Searching for Wannier seeds in strip metallic states at r=0...")
wannierregionfockR00 = wannierregionfock - normalvector*0.75
wannierregionfockR00|>getregion|>visualize

regionrestrictor = fourier(quasistripmetalcorrelations|>getoutspace, wannierregionfockR00, homemappings) / sqrt(quasistripmetalcorrelations|>getoutspace|>getcrystal|>vol)
localcorrelations = regionrestrictor' * quasistripmetalcorrelations * regionrestrictor
localcorrelations|>eigspech|>visualize

r00states = getregionstates(localcorrelations=localcorrelations, grouping=[2, 1, 1, 1, 1, 1, 1])|>collect
visualize(r00states[1], markersizemultiplier=20, markersizescaling=0.1)

visualize(m135 * r00states[7], markersizemultiplier=20, markersizescaling=0.1)
visualize(c2 * r00states[3], markersizemultiplier=20, markersizescaling=0.1)
visualize(m45 * r00states[2], markersizemultiplier=20, markersizescaling=0.1)

@info("Searching for Wannier seeds in strip metallic states at r=0.5...")
wannierregionfockR05 = wannierregionfock - normalvector*0.25
wannierregionfockR05|>getregion|>visualize
regionrestrictor = fourier(quasistripmetalcorrelations|>getoutspace, wannierregionfockR05, homemappings) / sqrt(quasistripmetalcorrelations|>getoutspace|>getcrystal|>vol)
localcorrelations = regionrestrictor' * quasistripmetalcorrelations * regionrestrictor
localcorrelations|>eigspech|>visualize
r05states = getregionstates(localcorrelations=localcorrelations, grouping=[4, 2, 2, 2])|>collect
visualize(m135 * r05states[1], markersizemultiplier=20, markersizescaling=0.1)
visualize(m45 * r05states[1], markersizemultiplier=20, markersizescaling=0.1)
visualize(c2 * r05states[1], markersizemultiplier=20, markersizescaling=0.1)

seeds = m135 * (r00states[2] + r00states[4] + r00states[5] + r00states[7]) + m135 * r05states[1]
visualize(seeds|>FockMap|>RegionState, markersizemultiplier=20, markersizescaling=0.1)

seedtransform = fourier(quasistripmetalcorrelations|>getoutspace, seeds|>FockMap|>getoutspace, homemappings) / sqrt(quasistripmetalcorrelations|>getoutspace|>getcrystal|>vol)
seedtransform|>visualize
seedtransform' * seedtransform |>visualize

crystalseeds = seedtransform * (seeds|>FockMap)
pseudoidentities = (crystalseeds[subspace, :]' * crystalseeds[subspace, :] for (_, subspace) in crystalseeds|>getoutspace|>crystalsubspaces)
(v for id in pseudoidentities for (_, v) in id|>eigvalsh)|>minimum

[v for id in pseudoidentities for (_, v) in id|>eigvalsh]|>sort

(crystalseeds' * crystalseeds)|>eigvalsh

crystalseeds = Dict(k=>crystalseeds[subspace, :] for (k, subspace) in crystalseeds|>getoutspace|>crystalsubspaces)
quasistripmetalstate|>geteigenvectors
quasistripmetalstate|>getcrystal

function Zipper.wannierprojection(; crystalisometries::Dict{Momentum, <:Zipper.FockMap}, crystal::Crystal, crystalseeds::Dict{Momentum, <:Zipper.FockMap}, svdorthothreshold::Number = 1e-1)
    wannierunitcell::Subset{Offset} = Subset(mode |> getattr(:b) for mode in (crystalseeds |> first |> last).inspace |> orderedmodes)
    wanniercrystal::Crystal = Crystal(wannierunitcell, crystal.sizes)
    overlaps = ((k, isometry, isometry' * crystalseeds[k]) for (k, isometry) in crystalisometries)
    precarioussvdvalues::Vector = []
    function approximateisometry(k::Momentum, isometry::FockMap, overlap::FockMap)::FockMap
        U, Σ, Vt = overlap |> svd
        minsvdvalue::Number = minimum(v for (_, v) in Σ)
        if minsvdvalue < svdorthothreshold
            push!(precarioussvdvalues, minsvdvalue)
        end
        unitary::FockMap = U * Vt
        approximated::FockMap = isometry * unitary

        inspace=Zipper.FockSpace(approximated|>getinspace|>setattr(:k=>k)|>removeattr(:r)|>mapmodes(m->m))
        return FockMap(approximated, inspace=inspace, performpermute=false)
    end
    if (precarioussvdvalues |> length) > 0
        @warn "Precarious wannier projection with minimum svdvalue of $(precarioussvdvalues |> minimum)"
    end
    blocks::Dict = Dict((k, k)=>approximateisometry(k, isometry, overlap) for (k, isometry, overlap) in overlaps)
    return CrystalFockMap(crystal, wanniercrystal, blocks)
end

wanniermetalisometries = wannierprojection(crystalisometries=quasistripmetalstate|>geteigenvectors, crystal=quasistripmetalstate|>getcrystal, crystalseeds=crystalseeds)
wanniermetalisometries|>getinspace|>unitcellfock|>modeattrs

visualfock = (wannierregionfockR00 + (wannierregionfockR00 + normalvector * 0.5))|>RegionFock
visualfock|>getregion|>visualize
eigenfock = wanniermetalisometries|>getinspace|>unitcellfock|>RegionFock
leftfourier = fourier(wanniermetalisometries|>getoutspace, visualfock, homemappings) / sqrt(wanniermetalisometries|>getoutspace|>getcrystal|>vol)
rightfourier = fourier(wanniermetalisometries|>getinspace, eigenfock)
wanniermetalisometries
wannierstates = leftfourier' * FockMap(wanniermetalisometries) * rightfourier
visualize(wannierstates[:, 1:4]|>RegionState, markersizemultiplier=30, markersizescaling=0.6)
visualize(wannierstates[:, 5:8]|>RegionState, markersizemultiplier=30, markersizescaling=0.6)

wannierstates|>getoutspace|>getregion|>visualize
visualize(getsphericalregion(crystal=blockedcrystal, radius=2, metricspace=blockedcrystal|>getspace), wannierstates|>getoutspace|>getregion)

function Zipper.spatialremapper(regionfock::RegionFock, region::Region)
    basispoints = (p=>p|>basispoint for p in region)
    positions::Dict{Offset, Tuple} = Dict(p=>(p - b, b) for (p, b) in basispoints)
    remappingdata::Base.Generator = ((positions[mode|>getpos], mode) for mode in regionfock)
    remappedfock::RegionFock = RegionFock(mode|>setattr(:r=>r, :b=>b) for ((r, b), mode) in remappingdata)
    return idmap(remappedfock, regionfock)
end

remapper = spatialremapper(wannierstates|>getoutspace, getsphericalregion(crystal=blockedcrystal, radius=2, metricspace=blockedcrystal|>getspace))
wannierstates = remapper * wannierstates
wanniermetalprojector = crystalprojector(localisometry=wannierstates, crystalfock=blockedcrystalfock)
c4transform = c4 * blockedcrystalfock
c4transform|>typeof
fullwanniermetalprojector = FockMap(wanniermetalprojector) + c4transform * FockMap(wanniermetalprojector) * c4transform'
fullwanniermetalprojector|>crystalspectrum|>visualize

courierstates = groundstatespectrum(fullwanniermetalprojector|>crystalspectrum, perunitcellfillings=4)
courierstates|>visualize

courierprojector = courierstates|>crystalprojector
couriercorrelations = idmap(courierprojector|>getoutspace) - courierprojector


restrictregionfock = (couriercorrelations|>getoutspace|>unitcellfock|>RegionFock) + ([0.5, 0]∈(blockedcrystal|>getspace))
visualize(couriercorrelations|>getoutspace|>unitcellfock|>RegionFock|>getregion, restrictregionfock|>getregion)
restrict = fourier(couriercorrelations|>getoutspace, restrictregionfock) / (blockedcrystal|>vol|>sqrt)
localcorrelations = restrict' * couriercorrelations * restrict
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[4])|>collect
visualize(localstates[1], markersizemultiplier=20, markersizescaling=0.5)

visualize(c4 * localstates[1], markersizemultiplier=20, markersizescaling=0.5)

blockedcrystal|>sitepoints|>visualize

function regionalrestriction(crystalstate::FockMap, regionfock::RegionFock)::RegionState
    eigenmodes::Subset{Mode} = crystalstate |> getinspace |> unitcellfock |> orderedmodes

    function extractregionstate(mode::Mode)
        rightfourier::FockMap = fourier(crystalstate |> getinspace, mode |> RegionFock)
        leftfourier::FockMap = fourier(crystalstate |> getoutspace, regionfock)
        return leftfourier' * crystalstate * rightfourier
    end

    return Dict(mode => mode |> extractregionstate for mode in eigenmodes) |> RegionState{crystalstate |> getoutspace |> getcrystal |> dimension}
end

stripcorrelations|>crystalspectrum|>linespectrum|>visualize
wanniermetalisometries' * stripcorrelations * wanniermetalisometries |>crystalspectrum|>linespectrum|>visualize

function regionalrestriction(crystalstate::FockMap, regionfock::RegionFock)::RegionState
    eigenmodes::Subset{Mode} = crystalstate |> getinspace |> unitcellfock |> orderedmodes

    function extractregionstate(mode::Mode)
        rightfourier::FockMap = fourier(crystalstate |> getinspace, mode |> RegionFock)
        leftfourier::FockMap = fourier(crystalstate |> getoutspace, regionfock)
        return leftfourier' * crystalstate * rightfourier
    end

    return Dict(mode => mode |> extractregionstate for mode in eigenmodes) |> RegionState{crystalstate |> getoutspace |> getcrystal |> dimension}
end
