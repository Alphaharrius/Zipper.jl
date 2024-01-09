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

m = quantize(:b, unitcell, 1) |> first

tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m, m |> setattr(:offset => [1, 0] ∈ square)) => tₙ,
    (m, m |> setattr(:offset => [0, 1] ∈ square)) => tₙ])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=0.5)
groundstates |> visualize
groundstateprojector = groundstates |> crystalprojector

C = idmap(groundstateprojector.outspace) - groundstateprojector

correlations = C

correlations |> crystalspectrum |> visualize

crystalfock = correlations |> getoutspace

scale = Scale([4 0; 0 4], square)
blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

H = energyspectrum |> FockMap
blockedH = blockresult[:transformer] * H * blockresult[:transformer]'
blockedH |> crystalspectrum |> visualize

blockedcrystal = blockresult[:crystal]

function extendedcrystalscale(; crystal::Crystal, generatingvector::Offset)::Scale
    snfinput::Matrix{Integer} = map(Integer, vcat(generatingvector |> vec |> transpose, crystal |> size |> diagm)) 
    _, S, Vd = snfinput |> dosnf
    return Scale((S |> diag |> diagm) * Vd |> transpose, crystal |> getspace)
end

blockedcorrelations = blockresult[:correlations]

blockedcrystal = blockresult[:crystal]

normalvector = (1, 1) ∈ (blockedcrystal|>getspace)

function getcrosssection(; crystal::Crystal, normalvector::Offset, radius::Real, metricspace::RealSpace = crystal|>getspace|>orthospace)
    height::Real = (*(metricspace, normalvector)|>norm)
    sphericalregion::Region = getsphericalregion(crystal=crystal, radius=sqrt(height^2 + radius^2)*2, metricspace=metricspace)

    normaldirection::Offset = *(metricspace, normalvector)|>normalize
    function crosssectionfilter(point::Point)::Bool
        metricpoint::Point = metricspace * point
        iswithinheight::Bool = 0 <= dot(metricpoint, normaldirection) <= height
        orthoreminder::Point = metricpoint - dot(metricpoint, normaldirection) * normaldirection
        iswithinradius::Bool = norm(orthoreminder) < radius
        return iswithinheight && iswithinradius
    end

    rawregion::Region = sphericalregion|>filter(crosssectionfilter)
    proximityregion::Region = rawregion - normalvector

    return rawregion - intersect(rawregion, proximityregion)
end

crosssection = getcrosssection(crystal=blockedcrystal, normalvector=normalvector, radius=0.5)

visualize(getsphericalregion(crystal=blockedcrystal, radius=5, metricspace=blockedcrystal|>getspace), crosssection, crosssection + normalvector)

function sitefock(site::Offset; flavorcount::Integer = 1)::FockSpace{Offset}
    basis::Offset = site|>basispoint
    offset::Offset = site - basis
    return FockSpace((Mode(:offset => offset, :b => basis, :flavor => f) for f in 1:flavorcount), reflected=site)
end

regionfock(region::Region; flavorcount::Integer = 1)::FockSpace{Region} = (
    fockspaceunion(sitefock(r, flavorcount=flavorcount) for r in region)|>FockSpace{Region})

crosssectionfock = regionfock(crosssection)
crosssectionfourier = fourier(blockedcorrelations|>getoutspace, crosssectionfock)

function extendedcrystalscale(; crystal::Crystal, generatingvector::Offset)
    snfinput::Matrix{Integer} = map(Integer, vcat(generatingvector|>vec|>transpose, crystal|>size|>diagm)) 
    _, S, Vd = snfinput|>dosnf
    return Scale((S|>diag|>diagm) * Vd |>transpose, crystal|>getspace)
end

extendedscale = extendedcrystalscale(crystal=blockedcrystal, generatingvector=normalvector)
extendedscale|>rep

scaledcrystal = extendedscale * blockedcrystal
scaledcrystal|>getspace|>rep
visualize(scaledcrystal|>sitepoints, scaledcrystal|>getunitcell)

scaledspace = scaledcrystal|>getspace

embeddedfock = FockSpace(mode|>setattr(:b=>(scaledspace * getpos(mode)))|>removeattr(:offset) for mode in crosssectionfock)
embeddedfock|>modeattrs

scaledkspace = convert(MomentumSpace, scaledspace)

blockedcrystalfock = blockedcorrelations|>getoutspace

using OrderedCollections
momentummappings = (scaledkspace * k |> basispoint => k for k in blockedcrystal|>brillouinzone)
mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point, Vector{Point}}()) do d, (k, v)
    mergewith!(append!, d, LittleDict(k => [v]))
end
ksubsets = blockedcrystalfock|>crystalsubsets
scaledfock::FockSpace = ((ksubsets[k] for k in mappingpartitions[kscaled])|>subsetunion|>FockSpace
        for kscaled in mappingpartitions|>keys)|>fockspaceunion
scaledbz = Subset(sk for sk in mappingpartitions|>keys)

embeddedunitcell = Subset(mode|>getattr(:b) for mode in embeddedfock)
embeddedunitcell|>collect

embeddedcrystal = Crystal(embeddedunitcell, scaledcrystal|>size)

volumeratio = vol(blockedcrystal) / vol(embeddedcrystal)

function repackfourierblocks(source::FockMap, kscaled::Momentum, partition::Subset{Mode})::FockMap
    partitionrows::FockMap = rows(source, partition |> FockSpace)
    inspace::FockSpace = FockSpace(setattr(mode, :offset => kscaled, :b => scaledspace * (mode|>getpos)) for mode in partitionrows|>getinspace)
    return FockMap(partitionrows.outspace, inspace, partitionrows |> rep) / sqrt(volumeratio)
end

repackedblocks::Base.Generator = (
    repackfourierblocks(crosssectionfourier, kscaled, partition)
    for (kscaled, partition) in Iterators.zip(scaledbz, scaledfock|>rep))
blockmap::FockMap = directsum(repackedblocks)
embeddedcrystalfock = FockSpace(blockmap|>getinspace, reflected=embeddedcrystal)
blockmap = FockMap(blockmap, inspace=embeddedcrystalfock, outspace=blockedcrystalfock)

extendedrestrictedC = blockmap' * blockedcorrelations * blockmap
# extendedrestrictedC|>visualize

extendedCspec = extendedrestrictedC|>crystalspectrum

extendedCspec|>linespectrum|>visualize

eigenmodes = sort([(k, modes) for (k, modes) in extendedCspec.eigenmodes], by=(v -> (v[1]|>vec)[2]))
cspec = hcat(([extendedCspec.eigenvalues[mode] for mode in modes]|>sort for (_, modes) in eigenmodes)...)
using Plotly
plot([scatter(y=cspec[n, :]) for n in axes(cspec, 1)])

extendedfrozenspec = distillation(extendedCspec, :frozen => v -> v < 1e-2 || v > 1 - 1e-2)[:frozen]
extendedfrozenprojector = extendedfrozenspec|>crystalprojector
extendedfrozencorrelations = idmap(extendedfrozenprojector|>getoutspace) - extendedfrozenprojector
extendedCspec = extendedfrozencorrelations|>crystalspectrum
eigenmodes = sort([(k, modes) for (k, modes) in extendedCspec.eigenmodes], by=(v -> (v[1]|>vec)[2]))
cspec = hcat(([extendedCspec.eigenvalues[mode] for mode in modes]|>sort for (_, modes) in eigenmodes)...)
using Plotly
plot([scatter(y=cspec[n, :]) for n in axes(cspec, 1)])

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
    outspace::CrystalFock = FockSpace((Mode(:offset => k, :b => barepoint, :flavor => 1) for k in bz), reflected=barecrystal)
    
    latticesites::Region = Subset(p - (p|>basispoint) for p in region) # Removing all sub-lattice degrees of freedom.
    inspace::FockSpace{Region} = latticesites|>regionfock

    momentummatrix::Matrix = hcat((k|>euclidean|>vec for k in bz)...)
    offsetmatrix::Matrix = hcat((p|>euclidean|>vec for p in latticesites)...)

    fouriermatrix::SparseMatrixCSC = exp.(-1im * momentummatrix' * offsetmatrix)
    return FockMap(outspace, inspace, fouriermatrix)
end

function crystaltoregionsublatticemapping(crystalfock::CrystalFock, regionfock::RegionFock)
    crystalunitfock = crystalfock|>unitcellfock|>orderedmodes|>removeattr(:offset)|>FockSpace
    regionunitfock = regionfock|>orderedmodes|>removeattr(:offset)|>FockSpace
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
modefrombasis = Dict((mode|>getattr(:b)|>basispoint) => mode for mode in extendedrestrictedC|>getoutspace|>unitcellfock)
truncunitfock = Subset(p|>basispoint for p in transformedcrosssection)|>regionfock
values = Dict((modefrombasis[rmode|>getattr(:b)], rmode) => 1 + 0im for rmode in truncunitfock)
modemap = FockMap(extendedrestrictedC|>getoutspace|>unitcellfock, truncunitfock, values)
# After this perform tensor product

function tensorproduct(primary::FockSpace, secondary::FockSpace; keepprimaryattrs)
    secondarymodes = (removeattr(mode, keepprimaryattrs...) for mode in secondary)
    return fockspaceunion(FockSpace(merge(pmode, smode) for smode in secondarymodes) for pmode in primary)
end

Ftinspace = tensorproduct(Ft|>getinspace, modemap|>getinspace, keepprimaryattrs=[:offset])|>FockSpace{Region}
Ftoutspace = tensorproduct(Ft|>getoutspace, modemap|>getoutspace, keepprimaryattrs=[:offset])
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

wannierregion = Subset(scaledspace * point for point in crosssection)
wannierregion|>regionfock

regioncorrelations(truncfilledcorrelations, wannierregion|>regionfock)|>eigspech|>visualize

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
selectedisometries = ((localspectrum |> geteigenvectors)[:, (group.first|>FockSpace)[1:8]] for group in selectedgroups)
selectedisometries|>collect
function lineardependencefilter(spstate::FockMap)::Bool
    crystalspstates::Dict{Momentum, FockMap} = crystalisometries(localisometry=spstate, crystalfock=truncfilledcorrelations|>getoutspace)
    crystalspstate::FockMap = directsum(v for (_, v) in crystalspstates)
    pseudoidentity::FockMap = (crystalspstate' * crystalspstate)
    mineigenvalue = minimum(v for (_, v) in pseudoidentity |> eigvalsh)
    println(mineigenvalue)
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

Ftinspace = tensorproduct(Ft|>getinspace, modemap|>getinspace, keepprimaryattrs=[:offset])|>FockSpace{Region}
Ftoutspace = tensorproduct(Ft|>getoutspace, modemap|>getoutspace, keepprimaryattrs=[:offset])
Ftoutspace = FockSpace(Ftoutspace, reflected=extendedrestrictedC|>getoutspace|>getcrystal)
Ftrep = kron(Ft|>rep, modemap|>rep)
newFt = FockMap(Ftoutspace, Ftinspace, Ftrep)
rFt = columns(newFt, transformedcrosssection|>regionfock)

rFt' * truncfilledcorrelations * rFt / 24 |>eigspech|>visualize
truncfilledprojector * extendedrestrictedC * truncfilledprojector |>crystalspectrum|>linespectrum|>visualize


comparefock = FockSpace(mode|>fixposition for mode in extendedrestrictedC|>getoutspace|>unitcellfock)



truncfock = FockSpace{Region}(mode|>setattr(:b => mode|>getpos)|>setattr(:offset => mode|>getpos|>getspace|>getorigin) for mode in regionfock(transformedcrosssection))
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
