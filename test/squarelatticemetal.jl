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

m = quantize(unitcell, 1) |> first

tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m, m |> setattr(:r => [1, 0] ∈ square)) => tₙ,
    (m, m |> setattr(:r => [0, 1] ∈ square)) => tₙ])


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

blocker = scale * crystalfock
blockedcorrelations = blocker * correlations * blocker'

H = energyspectrum |> FockMap
blockedH = blocker * H * blocker'
blockedH|>crystalspectrum|>visualize

blockedcrystal = blocker|>getoutspace|>getcrystal

function extendedcrystalscale(; crystal::Crystal, generatingvector::Offset)::Scale
    snfinput::Matrix{Integer} = map(Integer, vcat(generatingvector |> vec |> transpose, crystal |> size |> diagm)) 
    _, S, Vd = snfinput |> dosnf
    return Scale((S |> diag |> diagm) * Vd |> transpose, crystal |> getspace)
end

normalvector = (1, 1) ∈ (blockedcrystal|>getspace)

function getcrosssection(; crystal::Crystal, normalvector::Offset, radius::Real, metricspace::RealSpace = crystal|>getspace|>orthospace, minbottomheight::Real = 0)
    height::Real = (*(metricspace, normalvector)|>norm)
    sphericalregion::Region = getsphericalregion(crystal=crystal, radius=sqrt(height^2 + radius^2)*2, metricspace=metricspace)

    normaldirection::Offset = *(metricspace, normalvector)|>normalize
    function crosssectionfilter(point::Point)::Bool
        metricpoint::Point = metricspace * point
        iswithinheight::Bool = minbottomheight < dot(metricpoint, normaldirection) <= height
        orthoreminder::Point = metricpoint - dot(metricpoint, normaldirection) * normaldirection
        iswithinradius::Bool = norm(orthoreminder) <= radius
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

crosssectionfock = quantize(crosssection, 1)
crosssectionfourier = Zipper.fourier(blockedcorrelations|>getoutspace, crosssectionfock)

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
    inspace::FockSpace = FockSpace(setattr(mode, :k => kscaled, :b => scaledspace * (mode|>getpos))|>removeattr(:r) for mode in partitionrows|>getinspace)
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

extendedrestrictedC|>getoutspace|>getcrystal|>getspace|>rep

extendedfrozenspec = distillation(extendedCspec, :frozen => v -> v < 1e-2 || v > 1 - 1e-2)[:frozen]
extendedfrozenprojector = extendedfrozenspec|>crystalprojector
extendedfrozencorrelations = idmap(extendedfrozenprojector|>getoutspace) - extendedfrozenprojector
extendedCspec = extendedfrozencorrelations|>crystalspectrum

extendedCspec|>linespectrum|>visualize

scaledspace = extendedrestrictedC|>getoutspace|>getcrystal|>getspace
truncregion = reduce(+, crosssection + normalvector * n for n in 0:23)
transformedcrosssection = Subset(scaledspace * point for point in truncregion)

visualize(transformedcrosssection, crosssection)

# ======================================== Truncation of Strip correlations ========================================
Ft = Zipper.fourier(extendedrestrictedC|>getoutspace|>getcrystal, transformedcrosssection)
Ft|>visualize
# Compute unit-cell mapping FockMap
modefrombasis = Dict((mode|>getattr(:b)|>basispoint) => mode for mode in extendedrestrictedC|>getoutspace|>orderedmodes|>removeattr(:k))
truncregionfock = quantize(transformedcrosssection, 1)|>RegionFock
truncunitfock = truncregionfock|>unitcellfock
values = Dict((modefrombasis[rmode|>getattr(:b)], rmode) => 1 + 0im for rmode in truncunitfock)
modemap = FockMap(extendedrestrictedC|>getoutspace|>orderedmodes|>removeattr(:k)|>FockSpace, truncunitfock, values)
modemap|>visualize
# After this perform tensor product
rFt = Zipper.fourier(extendedrestrictedC|>getoutspace, truncregionfock, modemap) / 24
visualize(rFt)

realspaceextendedC = rFt' * extendedrestrictedC * rFt

keptmodepairs = ((a, b) for (a, b) in Iterators.product(realspaceextendedC|>getoutspace, realspaceextendedC|>getoutspace) if norm(euclidean(getpos(a) - getpos(b))) < 17)
keptmodepairs|>collect
trunccorrelationdata = spzero(realspaceextendedC|>getoutspace|>dimension, realspaceextendedC|>getoutspace|>dimension)


rFt * rFt' |>visualize

rFt[(rFt|>getoutspace|>rep)[10]|>FockSpace, :]' * rFt[(rFt|>getoutspace|>rep)[10]|>FockSpace, :] |>visualize

tiprojector = FockMap(directsum(rFt[subspace, :] * rFt[subspace, :]' for (k, subspace) in rFt|>getoutspace|>crystalsubspaces), inspace=rFt|>getoutspace, outspace=rFt|>getoutspace)
tiprojector|>visualize

truncatedextendedfrozencorrelations = tiprojector * extendedfrozencorrelations * tiprojector'
truncatedextendedfrozencorrelations|>visualize



truncextendedCspec = truncatedextendedfrozencorrelations|>crystalspectrum

truncextendedCspec|>linespectrum|>visualize

truncfilledCspec = groundstatespectrum(truncextendedCspec, perunitcellfillings=8)
truncfilledCspec|>linespectrum|>visualize

truncfilledprojector = truncfilledCspec|>crystalprojector
truncfilledcorrelations = idmap(truncfilledprojector|>getoutspace) - truncfilledprojector
truncfilledcorrelations|>crystalspectrum|>linespectrum|>visualize
# ======================================== Truncation of Strip correlations ========================================

# ========================================Strip Metal Wannierization========================================
normalvector*0.875|>norm
(normalvector|>norm)

p = [3.5,3.5]∈square
bp = getspace(blockedcrystal) * p
dot(bp, normalvector|>normalize)

function crosssectionlocalstate(wannierregion::Region, predicate::Function)
    wannierregionfock = quantize(wannierregion, 1)
    Ft = Zipper.fourier(truncfilledcorrelations|>getoutspace|>getcrystal, wannierregionfock|>getregion)
    rFt = kron(Ft, modemap)[:, wannierregionfock] / sqrt(Ft|>getoutspace|>dimension)
    wannierlocalcorrelations = rFt' * truncfilledcorrelations * rFt
    localspectrum = eigspech(wannierlocalcorrelations, groupingthreshold=1e-3)
    # localspectrum|>visualize

    groupeigenvalues::Base.Generator = (
        subset=>(localspectrum|>geteigenvalues)[subset|>first]
        for subset in localspectrum|>geteigenvectors|>getinspace|>sparsegrouping(:eigenindex)|>rep)
    selectedgroups = Iterators.filter(p -> p.second|>predicate, groupeigenvalues)
    selectedisometries = ((localspectrum|>geteigenvectors)[:, group.first|>FockSpace] for group in selectedgroups)
    selectedstate = selectedisometries|>first
    stateoutspace = selectedstate|>getoutspace
    stateoutregion = stateoutspace|>getregion
    # TODO: The recenter method probably need some look into.
    offmirror135 = mirror135|>recenter(stateoutregion|>getcenter)
    # visualize(stateoutregion, offmirror135 * stateoutregion)
    return selectedstate * offmirror135
end

function localcrosssectioncorrelations(wannierregion::Region)
    wannierregionfock = quantize(wannierregion, 1)
    Ft = Zipper.fourier(truncfilledcorrelations|>getoutspace|>getcrystal, wannierregionfock|>getregion)
    rFt = kron(Ft, modemap)[:, wannierregionfock] / sqrt(Ft|>getoutspace|>dimension)
    wannierlocalcorrelations = rFt' * truncfilledcorrelations * rFt
    return eigspech(wannierlocalcorrelations, groupingthreshold=1e-3)
end

# For r=1/2
basewanniercrosssection = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.875, radius=0.25, minbottomheight=0.125)
wanniercrosssection = basewanniercrosssection
trialwannierregion1 = Subset(scaledspace * point for point in wanniercrosssection)
visualize(crosssection, trialwannierregion1)

localcrosssectioncorrelations(trialwannierregion1)|>visualize
crosssectionlocalstate(trialwannierregion1, v -> v < 1e-3)|>getinspace|>modeattrs

basewanniercrosssection = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.87, radius=0.5, minbottomheight=0.25)
wanniercrosssection = basewanniercrosssection
trialwannierregion2 = Subset(scaledspace * point for point in wanniercrosssection)
visualize(crosssection, trialwannierregion2)

localcrosssectioncorrelations(trialwannierregion2)|>visualize
crosssectionlocalstate(trialwannierregion2, v -> v < 1e-3)|>getinspace|>modeattrs

basewanniercrosssection = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.75, radius=0.5, minbottomheight=0.5)
wanniercrosssection = basewanniercrosssection
trialwannierregion3 = Subset(scaledspace * point for point in wanniercrosssection)
visualize(crosssection, trialwannierregion3)

localcrosssectioncorrelations(trialwannierregion3)|>visualize
crosssectionlocalstate(trialwannierregion3, v -> v < 1e-3)|>getinspace|>modeattrs

trialwannierregionX = trialwannierregion1 + trialwannierregion2
visualize(crosssection, trialwannierregionX)

localcrosssectioncorrelations(trialwannierregionX)|>visualize

# For r=0
basewanniercrosssection = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.875, radius=0.6, minbottomheight=0.125)
wanniercrosssection = basewanniercrosssection - normalvector/2
trialwannierregion4 = Subset(scaledspace * point for point in wanniercrosssection)
visualize(crosssection, trialwannierregion4)

localcrosssectioncorrelations(trialwannierregion4)|>visualize
crosssectionlocalstate(trialwannierregion4, v -> v < 1e-3)|>getinspace|>modeattrs

basewanniercrosssection = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.87, radius=0.25, minbottomheight=0.25)
wanniercrosssection = basewanniercrosssection - normalvector/2
trialwannierregion5 = Subset(scaledspace * point for point in wanniercrosssection)
visualize(crosssection, trialwannierregion5)

localcrosssectioncorrelations(trialwannierregion5)|>visualize
crosssectionlocalstate(trialwannierregion5, true)|>getinspace|>modeattrs
crosssectionlocalstate(trialwannierregion5, false)|>getinspace|>modeattrs

basewanniercrosssection = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.75, radius=0.5, minbottomheight=0.5)
wanniercrosssection = basewanniercrosssection - normalvector/2
trialwannierregion6 = Subset(scaledspace * point for point in wanniercrosssection)
visualize(crosssection, trialwannierregion6)

crosssectionlocalstate(trialwannierregion6, true)|>getinspace|>modeattrs
crosssectionlocalstate(trialwannierregion6, false)|>getinspace|>modeattrs

offwanniercrosssection = wanniercrosssection + normalvector / 2
visualize(crosssection, offwanniercrosssection)
wannierregion = Subset(scaledspace * point for point in offwanniercrosssection)
wannierregionfock = quantize(wannierregion, 1)
Ft = Zipper.fourier(truncfilledcorrelations|>getoutspace|>getcrystal, wannierregionfock|>getregion)
rFt = kron(Ft, modemap)[:, wannierregionfock]
wannierlocalcorrelations2 = rFt' * truncfilledcorrelations * rFt
wannierlocalcorrelations2|>eigspech|>visualize

newwanniercrosssection = offwanniercrosssection - normalvector / 2
visualize(crosssection, newwanniercrosssection)

mirror135 = AffineTransform([-1 0; 0 -1], localspace=square)

visualize(crosssection, wanniercrosssection, offwanniercrosssection)

localspectrum2 = eigspech(wannierlocalcorrelations2, groupingthreshold=1e-3)
localspectrum2|>visualize

groupeigenvalues::Base.Generator = (
    subset=>(localspectrum|>geteigenvalues)[subset|>first]
    for subset in localspectrum|>geteigenvectors|>getinspace|>sparsegrouping(:eigenindex)|>rep)
selectedgroups = Iterators.filter(p -> p.second|>(v -> v > 0.98), groupeigenvalues)
selectedisometries = ((localspectrum|>geteigenvectors)[:, group.first|>FockSpace] for group in selectedgroups)
# offmirror135 = recenter(mirrow135, wanniercrosssection|>getcenter)
# visualize(wanniercrosssection, offmirror135 * wanniercrosssection)
# quantize(offmirror135 * wanniercrosssection, 1)|>getregion
# visualize(wanniercrosssection, quantize(offmirror135 * wanniercrosssection, 1)|>getregion)
# quantize(offmirror135 * wanniercrosssection, 1)|>modeattrs

function Base.:*(transformation::AffineTransform, regionfock::FockSpace{Region})::FockMap
    # This is used to correct the :b attribute, since the :b as a Point will be symmetrized,
    # which the basis point set might not include the symmetrized :b. Thus we would like to set
    # the :b to its corresponding basis point, and offload the difference to :r.
    function correctsymmetrizedmode(mode::Mode)::Mode
        actualposition::Offset = mode|>getattr(:R)
        basisposition::Offset = actualposition|>basispoint
        offset::Offset = actualposition - basisposition
        return mode|>setattr(:b=>basisposition)|>setattr(:r=>offset)|>removeattr(:R)
    end

    modemapping::Dict{Mode, Mode} = Dict()

    mergepositions(mode::Mode)::Mode = mode|>setattr(:R=>getpos(mode))|>removeattr(:r, :b)

    function modesymmetrize(mode::Mode)::Mode
        fixedmode = mode |> mergepositions
        newattrs::Dict{Symbol, Any} = Dict(fixedmode.attrs)
        # TODO: There are some attributes that are not meant to be transformed with the mode.
        filterpredicate = p -> hasmethod(*, Tuple{AffineTransform, p.second |> typeof})
        foreach(p -> newattrs[p.first] = transformation * p.second, Iterators.filter(filterpredicate, fixedmode.attrs))
        newmode::Mode = Mode(newattrs)|>correctsymmetrizedmode
        modemapping[newmode] = mode
        return newmode
    end

    connections::Dict{Tuple{Mode, Mode}, ComplexF64} = Dict()

    function rebaseorbital(mode::Mode)::Mode
        frommode::Mode = modemapping[mode]
        tomode::Mode = mode |> setattr(:orbital => (frommode |> getorbital()))
        connections[(tomode, frommode)] = (mode |> getorbital()) |> relativephase(frommode |> getorbital())
        return tomode
    end

    outmodes::Subset{Mode} = Subset(mode |> modesymmetrize |> rebaseorbital for mode in regionfock)
    
    return FockMap(outmodes |> FockSpace{Region}, regionfock, connections)
end

selectedstate = selectedisometries|>first
stateoutspace = state|>getoutspace
stateoutregion = stateoutspace|>getregion
visualize(stateoutregion)
# TODO: The recenter method probably need some look into.
offmirror135 = mirror135|>recenter(stateoutregion|>getcenter)
# visualize(stateoutregion, offmirror135 * stateoutregion)
selectedstate * offmirror135

symmetricspstates = (state * *(state, offmirror135) for state in selectedisometries)
spstates = (state * spatialmap(state)' for state in symmetricspstates)

localstate = Dict(state|>getinspace|>dimension => state for state in spstates)[4]
localstate|>getinspace|>modeattrs
localstate|>getoutspace|>modeattrs


groupeigenvalues::Base.Generator = (
    subset => (localspectrum2 |> geteigenvalues)[subset |> first]
    for subset in localspectrum2 |> geteigenvectors |> getinspace |> sparsegrouping(:eigenindex) |> rep)
selectedgroups = Iterators.filter(p -> p.second|>(v -> v < 1e-3), groupeigenvalues)

selectedisometries = ((localspectrum2|>geteigenvectors)[:, group.first|>FockSpace] for group in selectedgroups)
state = selectedisometries|>first
stateoutspace = state|>getoutspace
stateoutregion = stateoutspace|>getregion
# TODO: The recenter method probably need some look into.
offmirror135 = offmirror135|>recenter(stateoutregion|>getcenter)
visualize(stateoutregion, offmirror135 * stateoutregion)

*(state, offmirror135)

symmetricspstates = (state * *(state, offmirror135) for state in selectedisometries)
spstates = (state * spatialmap(state)' for state in symmetricspstates)

localstate = Dict(state|>getinspace|>dimension => state for state in spstates)[5]
localstate|>getinspace|>modeattrs
localstate|>getoutspace|>modeattrs


function findlocalstate(;
    regioncorrelations::FockMap{RegionFock, RegionFock},
    symmetry::AffineTransform = identitytransform(regioncorrelations|>getoutspace|>getspace|>dimension),
    spectrumfilterpredicate::Function,
    spectrumdegeneracythreshold::Real = 1e-3)

    localspectrum::EigenSpectrum = eigspech(regioncorrelations, groupingthreshold=spectrumdegeneracythreshold)
    groupeigenvalues::Base.Generator = (
        subset => (localspectrum|>geteigenvalues)[subset|>first]
        for subset in localspectrum|>geteigenvectors|>getinspace|>sparsegrouping(:eigenindex)|>rep)
    selectedgroups = Iterators.filter(p -> p.second|>spectrumfilterpredicate, groupeigenvalues)
    
    selectedisometries = ((localspectrum|>geteigenvectors)[:, group.first|>FockSpace] for group in selectedgroups)
    symmetricspstates = (state * *(state, symmetry) for state in selectedisometries)
    spstates = (state * spatialmap(state)' for state in symmetricspstates)

    return Dict(state|>getinspace|>dimension => state for state in spstates)
end
# ========================================Strip Metal Wannierization========================================


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
