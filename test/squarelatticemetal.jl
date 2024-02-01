using Revise, LinearAlgebra, AbstractAlgebra, SmithNormalForm
using Zipper

square = euclidean(RealSpace, 2)
point = [1/2, 1/2] ∈ square
spatialsnappingcalibration([point])

kspace = convert(MomentumSpace, square)

c4 = pointgrouptransform([0 -1; 1 0])

unitcell = Subset(point)
crystal = Crystal(unitcell, [96, 96])
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
    stripunitcell = getcrosssection(crystal=blockedcrystal, normalvector=normalvector, radius=0.5)
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

mirror135 = pointgrouptransform([0 -1; -1 0], localspace=square)
wannierregion = getcrosssection(crystal=blockedcrystal, normalvector=normalvector * 1.5, radius=0.5, minbottomheight=0.15)
wannierregion|>visualize
scaledspace = stripcorrelations|>getoutspace|>getcrystal|>getspace
scaledunitcell = Subset(scaledspace * r for r in stripunitcell)
wannierregionfock = quantize(Subset(scaledspace * r for r in wannierregion), 1)
remapper = spatialremapper(wannierregionfock, offsets=Subset(-(scaledspace * normalvector), scaledspace|>getorigin, scaledspace * normalvector), unitcell=scaledunitcell)
wannierregionfock = remapper|>getoutspace

function getregionstates(; localcorrelations::FockMap{RegionFock, <:FockSpace}, grouping::Vector{<:Integer})
    localspectrum::EigenSpectrum = localcorrelations|>eigspech
    sortedmodes = [mode for (mode, _) in sort(localspectrum|>geteigenvalues|>collect, by=last)]
    currentindex = 1
    states::Vector{RegionState} = []
    for count in grouping
        currentfock = FockSpace(sortedmodes[currentindex:currentindex + count - 1])
        push!(states, geteigenvectors(localspectrum)[:, currentfock]|>RegionState)
        currentindex += count
    end
    return states
end

@info("Searching for Wannier seeds in strip metallic states at r=0...")
wannierregionfockR00 = wannierregionfock - normalvector*0.25
wannierregionfockR00|>getregion|>visualize
regionrestrictor = fourier(quasistripmetalcorrelations|>getoutspace, wannierregionfockR00, homemappings) / sqrt(quasistripmetalcorrelations|>getoutspace|>getcrystal|>vol)
regionrestrictor|>visualize
localcorrelations = regionrestrictor' * quasistripmetalcorrelations * regionrestrictor
localcorrelations|>eigspech|>visualize

r00state = getregionstates(localcorrelations=localcorrelations, grouping=[8])|>collect|>first
r00state|>getinspace|>modeattrs

visualize(c2 * r00state, markersizemultiplier=20, markersizescaling=0.1)
FockSpace(m for (m, _) in r00states[1])|>modeattrs
c2 * r00state|>FockMap|>getinspace|>modeattrs

localstateR00 = findlocalstate(localcorrelations=localcorrelations, symmetry=mirror135, bandgroupingthreshold=1e-3, bandpredicate=v -> v < 1e-2)
localstateR00 = localstateR00[:,1:2] + localstateR00[:,4]
localstateR00|>getinspace|>modeattrs
visualize(localstateR00|>RegionState, markersizemultiplier=20, markersizescaling=0.1)

function Zipper.findeigenfunction(transformation::AffineTransform; eigenvalue::Number = 1)::BasisFunction
    eigenvaluekey::Tuple = hashablecomplex(eigenvalue |> Complex, 2)
    if !haskey(transformation.eigenfunctions, eigenvaluekey)
        error("No eigenfunction is found for eigenvalue $(eigenvalue)!")
    end
    return transformation.eigenfunctions[eigenvaluekey]
end

@info("Searching for Wannier seeds in strip metallic states at r=0.5...")
wannierregionfockR05 = wannierregionfock - normalvector * 0.25
wannierregionfockR05|>getregion|>visualize
regionrestrictor = fourier(quasistripmetalcorrelations|>getoutspace, wannierregionfockR05, homemappings) / sqrt(quasistripmetalcorrelations|>getoutspace|>getcrystal|>vol)
regionrestrictor|>visualize
localcorrelations = regionrestrictor' * quasistripmetalcorrelations * regionrestrictor
localcorrelations|>eigspech|>visualize
localstateR05 = getregionstates(localcorrelations=localcorrelations, grouping=[8])|>collect|>first
localstateR05 = findlocalstate(localcorrelations=localcorrelations, symmetry=mirror135, bandgroupingthreshold=1e-3, bandpredicate=v -> v < 1e-2)
visualize(mirror135 * localstateR05, markersizemultiplier=20, markersizescaling=0.1)

seeds = r00state|>FockMap
seeds = FockMap(seeds, outspace=seeds|>getoutspace|>RegionFock, inspace=seeds|>getinspace|>RegionFock)

visualize(RegionState(seeds), markersizemultiplier=15, markersizescaling=0.08)

scaledspace = stripcorrelations|>getoutspace|>getcrystal|>getspace
scaledunitcell = Subset(scaledspace * r for r in stripunitcell)
offsets = Subset(normalvector*0, -normalvector, normalvector)
scaledoffsets = Subset(scaledspace * off for off in offsets)
remapper = spatialremapper(seeds|>getoutspace|>RegionFock, offsets=scaledoffsets, unitcell=scaledunitcell)
remappedseeds = remapper * seeds

remappedseeds|>getoutspace|>getregion|>visualize

seedtransform = fourier(quasistripmetalcorrelations|>getoutspace, seeds|>getoutspace, homemappings) / sqrt(quasistripmetalcorrelations|>getoutspace|>getcrystal|>vol)
seedtransform|>visualize
seedtransform' * seedtransform |>visualize
seedtransform|>visualize

crystalseeds = seedtransform * seeds
pseudoidentities = (crystalseeds[subspace, :]' * crystalseeds[subspace, :] for (_, subspace) in crystalseeds|>getoutspace|>crystalsubspaces)
(v for id in pseudoidentities for (_, v) in id|>eigvalsh)|>minimum
(crystalseeds' * crystalseeds)|>eigvalsh
minimum(v for (_, v) in pseudoidentity|>eigvalsh)
