using Revise, LinearAlgebra, AbstractAlgebra, SmithNormalForm
using Zipper

# a = [1, 1]
# b = [-1, -1]

# T = a * b'

# U, S, V = svd(T)
# M = U * V'
# det(M)
# M = map(v -> round(v, digits=10), U * V')
# U, _, V = svd(M)
# V * b
# V * [-24, 0]

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

function getcrosssection(; crystal::Crystal, normalvector::Offset, radius::Real, metricspace::RealSpace = crystal|>getspace|>orthospace, minbottomheight::Real = 0)
    height::Real = (*(metricspace, normalvector)|>norm)
    sphericalregion::Region = getsphericalregion(crystal=crystal, radius=sqrt(height^2 + radius^2)*2, metricspace=metricspace)

    normaldirection::Offset = *(metricspace, normalvector)|>normalize
    function crosssectionfilter(point::Point)::Bool
        metricpoint::Point = metricspace * point
        iswithinheight::Bool = minbottomheight <= dot(metricpoint, normaldirection) <= height
        orthoreminder::Point = metricpoint - dot(metricpoint, normaldirection) * normaldirection
        iswithinradius::Bool = norm(orthoreminder) < radius
        return iswithinheight && iswithinradius
    end

    rawregion::Region = sphericalregion|>filter(crosssectionfilter)
    proximityregion::Region = rawregion - normalvector

    return rawregion - intersect(rawregion, proximityregion)
end

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

truncationregioncorrelations = truncator' * stripfrozencorrelations * truncator

@info("Performing truncation on strip frozen correlations...")
@time begin
    truncationregionindices = Iterators.product(truncationregionfock, truncationregionfock)
    truncationregionbonds = Dict((getattr(tomode, :r) - getattr(frommode, :r), frommode|>removeattr(:r), tomode|>removeattr(:r)) => (frommode, tomode) for (frommode, tomode) in truncationregionindices)
    pruningindices = (index for (_, index) in truncationregionbonds)
    prunedcorrelations = extractindices(truncationregioncorrelations, pruningindices)
    stripfouriers = (truncator[subspace, :] for (_, subspace) in truncator|>getoutspace|>crystalsubspaces)
    stripcrystal = truncator|>getoutspace|>getcrystal
    truncatedstripfrozencorrelations = crystaldirectsum((transform * prunedcorrelations * transform' for transform in stripfouriers), outcrystal=stripcrystal, incrystal=stripcrystal)
    truncatedstripfrozencorrelationspectrum = truncatedstripfrozencorrelations|>crystalspectrum
end
truncatedstripfrozencorrelationspectrum|>linespectrum|>visualize

@info("Retrieving strip metallic states...")
@time begin
    stripmetalstates = @time groundstatespectrum(truncatedstripfrozencorrelationspectrum, perunitcellfillings=8)
    stripmetalprojector = stripmetalstates|>crystalprojector
    stripmetalcorrelations = idmap(stripmetalprojector|>getoutspace) - stripmetalprojector
end
stripmetalstates|>linespectrum|>visualize

function Base.:+(regionfock::RegionFock, offset::Offset)::RegionFock
    positiontomodes = ((getpos(mode) + offset)=>mode for mode in regionfock)
    return RegionFock(mode|>setattr(:r=>(pos - basispoint(pos)), :b=>basispoint(pos)) for (pos, mode) in positiontomodes)
end

Base.:-(regionfock::RegionFock, offset::Offset)::RegionFock = regionfock + (-offset)

mirror135 = pointgrouptransform([0 -1; -1 0], localspace=square)
wannierregion = getcrosssection(crystal=blockedcrystal, normalvector=normalvector * 0.875, radius=0.5, minbottomheight=0.15)
wannierregionfock = RegionFock(mode for mode in truncationregionfock if mode|>getpos in wannierregion)    

function findlocalstate(; localcorrelations::FockMap, symmetry::AffineTransform, bandgroupingthreshold::Real = 1e-3, bandpredicate::Function)
    localseedingspectrum::EigenSpectrum = eigspech(localcorrelations, groupingthreshold=bandgroupingthreshold)
    groupeigenvalues::Base.Generator = (
        subset => (localseedingspectrum|>geteigenvalues)[subset|>first]
        for subset in localseedingspectrum|>geteigenvectors|>getinspace|>sparsegrouping(:eigenindex)|>rep)
    selectedgroups = Iterators.filter(p -> p.second|>bandpredicate, groupeigenvalues)
    selectedgroup = selectedgroups|>first
    selectedfock = selectedgroup|>first|>FockSpace
    selectedlocalstate = (localseedingspectrum|>geteigenvectors)[:, selectedfock]

    localstateregion = selectedlocalstate|>getoutspace|>getregion
    localsymmetry = symmetry|>recenter(localstateregion|>getcenter)
    transform = selectedlocalstate * localsymmetry
    symmetriclocalstate = selectedlocalstate * transform

    return symmetriclocalstate * spatialmap(symmetriclocalstate)'
end

@info("Searching for Wannier seeds in strip metallic states at r=0...")

function Base.:*(fockmap::FockMap, symmetry::AffineTransform)::FockMap
    outspacerep::FockMap = *(symmetry, fockmap |> getoutspace)
    hassamespan(outspacerep |> getoutspace, fockmap |> getoutspace) || error("The symmetry action on the outspace in not closed!")
    inspacerep::FockMap = fockmap' * outspacerep * fockmap

    phasespectrum::EigenSpectrum = inspacerep |> eigspec
    inspace::FockSpace = FockSpace(
        m |> setattr(:orbital => findeigenfunction(symmetry, eigenvalue=(phasespectrum |> geteigenvalues)[m]))
        # TODO: What if there are more attributes exists in the home modes? Such as :flavor, :orbital, etc.
        #   |> removeattr(:eigenindex) # The :orbital can subsitute the :eigenindex.
        for m in phasespectrum |> geteigenvectors |> getinspace)
    return FockMap(phasespectrum |> geteigenvectors, inspace=inspace, performpermute=false)
end

    wannierregionfockR00 = wannierregionfock - normalvector*0.5
    wannierregionfockR00|>getregion|>visualize
    regionrestrictor = fourier(stripmetalcorrelations|>getoutspace, wannierregionfockR00, homemappings) / sqrt(stripmetalcorrelations|>getoutspace|>getcrystal|>vol)
    localseedingcorrelations = regionrestrictor' * stripmetalcorrelations * regionrestrictor
    localseedingspectrum::EigenSpectrum = eigspech(localseedingcorrelations, groupingthreshold=1e-3)
    localseedingspectrum|>visualize
    groupeigenvalues::Base.Generator = (
        subset => (localseedingspectrum|>geteigenvalues)[subset|>first]
        for subset in localseedingspectrum|>geteigenvectors|>getinspace|>sparsegrouping(:eigenindex)|>rep)
    selectedgroups = Iterators.filter(p -> p.second|>(v -> v < 1e-2), groupeigenvalues)
    selectedgroup = selectedgroups|>first
    selectedfock = selectedgroup|>first|>FockSpace
    selectedlocalstateR00 = (localseedingspectrum|>geteigenvectors)[:, selectedfock]

    localstateregion = selectedlocalstateR00|>getoutspace|>getregion
    mirror135 = mirror135|>recenter(localstateregion|>getcenter)
    transform = selectedlocalstateR00 * mirror135
    symmetriclocalstateR00 = selectedlocalstateR00 * transform

localseedingspectrum|>visualize

visualize(RegionState(symmetriclocalstateR00), markersizemultiplier=20, markersizescaling=0.1)
visualize(RegionState(symmetriclocalstateR05), markersizemultiplier=20, markersizescaling=0.1)

symmetriclocalstateR00 * spatialmap(symmetriclocalstateR00)' |>getinspace|>modeattrs

@info("Searching for Wannier seeds in strip metallic states at r=0.5...")
@time begin
    regionrestrictor = fourier(stripmetalcorrelations|>getoutspace, wannierregionfock, homemappings) / sqrt(stripmetalcorrelations|>getoutspace|>getcrystal|>vol)
    localseedingcorrelations = regionrestrictor' * stripmetalcorrelations * regionrestrictor
    localseedingspectrum::EigenSpectrum = eigspech(localseedingcorrelations, groupingthreshold=6e-6)
    groupeigenvalues::Base.Generator = (
        subset => (localseedingspectrum|>geteigenvalues)[subset|>first]
        for subset in localseedingspectrum|>geteigenvectors|>getinspace|>sparsegrouping(:eigenindex)|>rep)
    selectedgroups = Iterators.filter(p -> p.second|>(v -> v < 1e-2), groupeigenvalues)
    selectedgroup = selectedgroups|>first
    selectedfock = selectedgroup|>first|>FockSpace
    selectedlocalstateR05 = (localseedingspectrum|>geteigenvectors)[:, selectedfock]

    localstateregion = selectedlocalstateR05|>getoutspace|>getregion
    mirror135 = mirror135|>recenter(localstateregion|>getcenter)
    transform = selectedlocalstateR05 * mirror135
    symmetriclocalstateR05 = selectedlocalstateR05 * transform
end
localseedingspectrum|>visualize

visualize(selectedlocalstateR00|>getoutspace|>getregion, selectedlocalstateR05|>getoutspace|>getregion)

mirror135 = pointgrouptransform([-1 0; 0 -1], localspace=square)
localstateregion = selectedlocalstate|>getoutspace|>getregion
mirror135 = mirror135|>recenter(localstateregion|>getcenter)
selectedlocalstate * mirror135 |>getinspace|>modeattrs

visualize(getsphericalregion(crystal=blockedcrystal, radius=5, metricspace=blockedcrystal|>getspace), crosssection, crosssection + normalvector)

crosssectionfock = quantize(crosssection, 1)
crosssectionfourier = fourier(blockedcorrelations|>getoutspace, crosssectionfock)

extendedscale = extendedcrystalscale(crystal=blockedcrystal, generatingvector=normalvector)

scaledcrystal = extendedscale * blockedcrystal

scaledspace = scaledcrystal|>getspace

embeddedfock = FockSpace(m|>setattr(:R=>(scaledspace * getpos(m)))|>removeattr(:r, :b) for m in crosssectionfock)

scaledkspace = convert(MomentumSpace, scaledspace)

using OrderedCollections
momentummappings = (scaledkspace * k |> basispoint => k for k in blockedcrystal|>brillouinzone)
mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point, Vector{Point}}()) do d, (k, v)
    mergewith!(append!, d, LittleDict(k => [v]))
end
ksubsets = blockedcrystalfock|>crystalsubsets
scaledfock::FockSpace = ((ksubsets[k] for k in mappingpartitions[kscaled])|>subsetunion|>FockSpace
        for kscaled in mappingpartitions|>keys)|>fockspaceunion
scaledbz = Subset(sk for sk in mappingpartitions|>keys)

embeddedunitcell = Subset(mode|>getattr(:R) for mode in embeddedfock)

embeddedcrystal = Crystal(embeddedunitcell, scaledcrystal|>size)

volumeratio = vol(blockedcrystal) / vol(embeddedcrystal)

function repackfourierblocks(source::FockMap, kscaled::Momentum, partition::Subset{Mode})::FockMap
    partitionrows::FockMap = rows(source, partition |> FockSpace)
    inspace::FockSpace = FockSpace(setattr(mode, :r => kscaled, :b => scaledspace * (mode|>getpos)) for mode in partitionrows|>getinspace)
    return FockMap(partitionrows.outspace, inspace, partitionrows |> rep) / sqrt(volumeratio)
end

repackedblocks::Base.Generator = (
    repackfourierblocks(crosssectionfourier, kscaled, partition)
    for (kscaled, partition) in Iterators.zip(scaledbz, scaledfock|>rep))

Subset(m|>getattr(:k) for m in repackedblocks|>collect|>first|>getoutspace)|>collect
blockmap::FockMap = crystaldirectsum(repackedblocks, outcrystal=blockedcrystal, incrystal=embeddedcrystal)
blockmap = FockMap(blockmap, inspace=embeddedcrystalfock, outspace=blockedcrystalfock)

extendedrestrictedC = blockmap' * blockedcorrelations * blockmap
# extendedrestrictedC|>visualize

extendedCspec = extendedrestrictedC|>crystalspectrum

extendedCspec|>linespectrum|>visualize

extendedfrozenspec = distillation(extendedCspec, :frozen => v -> v < 1e-2 || v > 1 - 1e-2)[:frozen]
extendedfrozenprojector = extendedfrozenspec|>crystalprojector
extendedfrozencorrelations = idmap(extendedfrozenprojector|>getoutspace) - extendedfrozenprojector
extendedCspec = extendedfrozencorrelations|>crystalspectrum

extendedCspec|>linespectrum|>visualize

scaledspace = extendedrestrictedC|>getoutspace|>getcrystal|>getspace
truncregion = reduce(+, crosssection + normalvector * n for n in 0:2)
transformedcrosssection = Subset(scaledspace * point for point in truncregion)

visualize(blockedcrystal|>sitepoints, transformedcrosssection)
visualize(blockedcrystal|>sitepoints, blockedcrystal|>getunitcell)

using SparseArrays
function Zipper.fourier(crystal::Crystal, region::Region)
    bz::Subset{Momentum} = crystal|>brillouinzone
    barepoint::Offset = crystal|>getspace|>getorigin
    barecrystal::Crystal = Crystal(barepoint|>Subset, crystal|>size)
    outspace::CrystalFock = FockSpace((Mode(:r => k, :b => barepoint, :flavor => 1) for k in bz), reflected=barecrystal)
    
    latticesites::Region = Subset(p - (p|>basispoint) for p in region) # Removing all sub-lattice degrees of freedom.
    inspace::FockSpace{Region} = latticesites|>regionfock

    momentummatrix::Matrix = hcat((k|>euclidean|>vec for k in bz)...)
    offsetmatrix::Matrix = hcat((p|>euclidean|>vec for p in latticesites)...)

    fouriermatrix::SparseMatrixCSC = exp.(-1im * momentummatrix' * offsetmatrix)
    return FockMap(outspace, inspace, fouriermatrix)
end

function crystaltoregionsublatticemapping(crystalfock::CrystalFock, regionfock::RegionFock)
    crystalunitfock = crystalfock|>unitcellfock|>orderedmodes|>removeattr(:r)|>FockSpace
    regionunitfock = regionfock|>orderedmodes|>removeattr(:r)|>FockSpace
    if !hasamespan(crystalunitfock, regionunitfock)
        @assert(==(crystalunitfock|>dimension, regionunitfock|>dimension))
        @warn("Abstract mapping between sub-lattice degrees of freedom...")
    end
    return idmap(crystalunitfock, regionunitfock)
end

function tensorproduct(primary::FockSpace, secondary::FockSpace)
    return fockspaceunion(FockSpace(merge(pmode, smode) for smode in secondary) for pmode in primary)
end

function Zipper.fourier(crystal::Crystal, regionfock::RegionFock)
end

function Zipper.fourier(crystalfock::CrystalFock, regionfock::RegionFock, sublatticemapping::FockMap)
    baretransform::FockMap = fourier(crystalfock|>getcrystal, regionfock|>getregion)
    ftrep::SparseMatrixCSC = kron(baretransform|>rep, sublatticemapping|>rep)
    return FockSpace(crystalfock, regionfock, ftrep)
end

transformedcrosssection
Ft = fourier(extendedrestrictedC|>getoutspace|>getcrystal, transformedcrosssection)
Ft|>getinspace|>modeattrs
Ft|>visualize
# Compute unit-cell mapping FockMap
modefrombasis = Dict((mode|>getattr(:b)|>basispoint) => mode for mode in extendedrestrictedC|>getoutspace|>orderedmodes|>removeattr(:r))
truncunitfock = Subset(p|>basispoint for p in transformedcrosssection)|>regionfock
values = Dict((modefrombasis[rmode|>getattr(:b)], rmode) => 1 + 0im for rmode in truncunitfock)
modemap = FockMap(extendedrestrictedC|>getoutspace|>orderedmodes|>removeattr(:r)|>FockSpace, truncunitfock, values)
# After this perform tensor product

function tensorproduct(primary::FockSpace, secondary::FockSpace; keepprimaryattrs)
    secondarymodes = (removeattr(mode, keepprimaryattrs...) for mode in secondary)
    return fockspaceunion(FockSpace(merge(pmode, smode) for smode in secondarymodes) for pmode in primary)
end

Ftinspace = tensorproduct(Ft|>getinspace, modemap|>getinspace, keepprimaryattrs=[:r])|>FockSpace{Region}
Ftoutspace = tensorproduct(Ft|>getoutspace, modemap|>getoutspace, keepprimaryattrs=[:r])
Ftoutspace = FockSpace(Ftoutspace, reflected=extendedrestrictedC|>getoutspace|>getcrystal)
Ftrep = kron(Ft|>rep, modemap|>rep)
newFt = FockMap(Ftoutspace, Ftinspace, Ftrep)
rFt = columns(newFt, transformedcrosssection|>regionfock)

tiisometries = (rows(rFt, fockspace) for (_, fockspace) in rFt|>getoutspace|>crystalsubspaces)
tiprojector = directsum(isometry * isometry' for isometry in tiisometries)
tiprojector = FockMap(tiprojector, inspace=rFt|>getoutspace, outspace=rFt|>getoutspace)

tiprojector = rFt * rFt'

truncatedextendedfrozencorrelations = tiprojector * extendedfrozencorrelations * tiprojector'
truncextendedCspec = truncatedextendedfrozencorrelations|>crystalspectrum

truncextendedCspec|>linespectrum|>visualize

truncfilledCspec = groundstatespectrum(truncextendedCspec, perunitcellfillings=8)
truncfilledCspec|>linespectrum|>visualize

truncfilledprojector = truncfilledCspec|>crystalprojector
truncfilledcorrelations = idmap(truncfilledprojector|>getoutspace) - truncfilledprojector
truncfilledcorrelations|>crystalspectrum|>linespectrum|>visualize

crosssectionfock|>modeattrs
truncfilledcorrelations|>getoutspace|>modeattrs

wanniercrosssection = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.5, radius=0.5)
wannierregion = Subset(scaledspace * point for point in wanniercrosssection)
regioncorrelations(truncfilledcorrelations, wannierregion|>regionfock)|>eigspech|>visualize

wanniercrosssection2 = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.35, radius=0.5) + normalvector * 0.5
wannierregion2 = Subset(scaledspace * point for point in wanniercrosssection2)
regioncorrelations(truncfilledcorrelations, wannierregion2|>regionfock)|>eigspech|>visualize

visualize(crosssection, wanniercrosssection, wanniercrosssection2)

function _findlocalspstates(;
    statecorrelations::FockMap, regionfock::FockSpace,
    symmetry::AffineTransform = identitytransform(statecorrelations|>getoutspace|>getcrystal|>dimension),
    spectrumextractpredicate::Function = v -> v < 1e-2,
    linearindependencethreshold::Real = 5e-2,
    degeneracythreshold::Real = 1e-7)::Dict{Integer, FockMap}

    function lineardependencefilter(spstate::FockMap)::Bool
        crystalspstates::Dict{Momentum, FockMap} = crystalisometries(localisometry=spstate, crystalfock=statecorrelations.outspace)
        crystalspstate::FockMap = directsum(v for (_, v) in crystalspstates)
        pseudoidentity::FockMap = (crystalspstate' * crystalspstate)
        mineigenvalue = minimum(v for (_, v) in pseudoidentity |> eigvalsh)
        return mineigenvalue > linearindependencethreshold
    end

    localcorrelations::FockMap = regioncorrelations(statecorrelations, regionfock)
    localspectrum::EigenSpectrum = eigspech(localcorrelations, groupingthreshold=degeneracythreshold)
    groupeigenvalues::Base.Generator = (
        subset => (localspectrum |> geteigenvalues)[subset |> first]
        for subset in localspectrum |> geteigenvectors |> getinspace |> sparsegrouping(:eigenindex) |> rep)
    selectedgroups = Iterators.filter(p -> p.second |> spectrumextractpredicate, groupeigenvalues)

    selectedisometries = ((localspectrum |> geteigenvectors)[:, group.first |> FockSpace] for group in selectedgroups)
    orthogonalspstates = Iterators.filter(lineardependencefilter, selectedisometries)
    symmetricspstates = (state * *(state, symmetry) for state in orthogonalspstates)
    spstates = (state * spatialmap(state)' for state in symmetricspstates)

    return Dict(state |> getinspace |> dimension => state for state in spstates)
end

wregioncorr = regioncorrelations(truncfilledcorrelations, wannierregion|>regionfock)
localspectrum::EigenSpectrum = eigspech(wregioncorr, groupingthreshold=1e-4)
localspectrum|>visualize
groupeigenvalues::Base.Generator = (
    subset => (localspectrum |> geteigenvalues)[subset |> first]
    for subset in localspectrum |> geteigenvectors |> getinspace |> sparsegrouping(:eigenindex) |> rep)
[groupeigenvalues...]
selectedgroups = Iterators.filter(p -> p.second |> (v -> v < 1e-2), groupeigenvalues)
[selectedgroups...]
selectedisometries = ((localspectrum |> geteigenvectors)[:, (group.first|>FockSpace)] for group in selectedgroups)
selectedisometries|>collect
function lineardependencefilter(spstate::FockMap)::Bool
    
    
    crystalspstates::Dict{Momentum, FockMap} = crystalisometries(localisometry=spstate, crystalfock=truncfilledcorrelations|>getoutspace)
    crystalspstate::FockMap = directsum(v for (_, v) in crystalspstates)
    pseudoidentity::FockMap = (crystalspstate' * crystalspstate)
    mineigenvalue = minimum(v for (_, v) in pseudoidentity |> eigvalsh)
    println([v for (_, v) in pseudoidentity |> eigvalsh])
    return mineigenvalue > 5e-2
end
orthogonalspstates = Iterators.filter(lineardependencefilter, selectedisometries)
symmetricspstates = (state * *(state, symmetry) for state in orthogonalspstates)
spstates = (state * spatialmap(state)' for state in symmetricspstates)
Dict(state |> getinspace |> dimension => state for state in spstates)

localseed = findlocalspstates(statecorrelations=truncfilledcorrelations, regionfock=wannierregion|>regionfock, spectrumextractpredicate=v -> v < 5e-1, degeneracythreshold=0.2)



truncregion = crosssection
transformedcrosssection = Subset(scaledspace * point for point in truncregion)
Ft = fourier(extendedrestrictedC|>getoutspace|>getcrystal, transformedcrosssection)
# Compute unit-cell mapping FockMap
modefrombasis = Dict((mode|>getattr(:b)|>basispoint) => mode for mode in extendedrestrictedC|>getoutspace|>unitcellfock)
truncunitfock = Subset(p|>basispoint for p in transformedcrosssection)|>regionfock
values = Dict((modefrombasis[rmode|>getattr(:b)], rmode) => 1 + 0im for rmode in truncunitfock)
modemap = FockMap(extendedrestrictedC|>getoutspace|>unitcellfock, truncunitfock, values)
# After this perform tensor product

Ftinspace = tensorproduct(Ft|>getinspace, modemap|>getinspace, keepprimaryattrs=[:r])|>FockSpace{Region}
Ftoutspace = tensorproduct(Ft|>getoutspace, modemap|>getoutspace, keepprimaryattrs=[:r])
Ftoutspace = FockSpace(Ftoutspace, reflected=extendedrestrictedC|>getoutspace|>getcrystal)
Ftrep = kron(Ft|>rep, modemap|>rep)
newFt = FockMap(Ftoutspace, Ftinspace, Ftrep)
rFt = columns(newFt, transformedcrosssection|>regionfock)

rFt' * truncfilledcorrelations * rFt / 24 |>eigspech|>visualize
truncfilledprojector * extendedrestrictedC * truncfilledprojector |>crystalspectrum|>linespectrum|>visualize


comparefock = FockSpace(mode|>fixposition for mode in extendedrestrictedC|>getoutspace|>unitcellfock)



truncfock = FockSpace{Region}(mode|>setattr(:b => mode|>getpos)|>setattr(:r => mode|>getpos|>getspace|>getorigin) for mode in regionfock(transformedcrosssection))
truncfourier = fourier(extendedfrozencorrelations|>getoutspace, truncfock)
truncfourier|>visualize
truncator = truncfourier * truncfourier'

truncatedextendedfrozencorrelations = truncator * extendedfrozencorrelations * truncator'
truncextendedCspec = truncatedextendedfrozencorrelations|>crystalspectrum
cspec = hcat(([truncextendedCspec.eigenvalues[mode] for mode in modes]|>sort for (_, modes) in truncextendedCspec.eigenmodes)...)
using Plotly
plot([scatter(y=cspec[n, :]) for n in axes(cspec, 1)])

extendedfrozencorrelationprime = blockmap * extendedfrozencorrelations * blockmap'
extendedfrozencorrelationprime|>crystalspectrum|>visualize
truncrosssectionfourier = fourier(blockedcorrelations|>getoutspace, regionfock(crosssection)) / sqrt(vol(blockedcorrelations|>getoutspace|>getcrystal))

truncator = truncrosssectionfourier * truncrosssectionfourier'

trunextendedfrozencorrelationprime = truncator * extendedfrozencorrelationprime * truncator'
trunextendedfrozencorrelationprime|>crystalspectrum|>visualize

trunextendedfrozencorrelations = blockmap' * trunextendedfrozencorrelationprime * blockmap
trunextendedCspec = trunextendedfrozencorrelations|>crystalspectrum

cspec = hcat(([trunextendedCspec.eigenvalues[mode] for mode in modes]|>sort for (_, modes) in trunextendedCspec.eigenmodes)...)
using Plotly
plot([scatter(y=cspec[n, :]) for n in axes(cspec, 1)])

extendedscale * crystal |>getunitcell|>collect
extendedscale * crystal |>getspace|>rep

scale * crystal|>getunitcell|>collect
scale * crystal |>getspace|>rep

extendedscale * ([1, 1] ∈ square)
