module Renormalization

using LinearAlgebra, OrderedCollections
using ..Spaces, ..Geometries, ..Quantum, ..Transformations, ..QuantumTransformations

"""
    frozenselectionbythreshold(threshold::Float64)

A selection strategy for `localfrozenisometries` and `globaldistillerhamiltonian` that selects the modes with correlation eigenvalues γ
that is localized within the local region ∀ γ ≈ 0 or γ ≈ 1 up to a `threshold`.

### Input
- `threshold::Float64`: The threshold for the correlation eigenvalues.

### Output
A function representing the selection strategy with an input of the local correlations 𝐶ᵣ `FockMap` and outputs a `Dict{Symbol, FockMap}`
keyed by two grouping symbols with their associated local isometries selected from the unitary that diagonalizes 𝐶ᵣ. The grouping symbol
`:filled` represents the eigenmodes that corresponds to γ ≈ 0; `:empty` represents the eigenmodes that corresponds to γ ≈ 1.
"""
function frozenselectionbythreshold(threshold::Float64)::Function
    function frozenfockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
        spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
        filledmodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second < threshold, spectrum |> geteigenvalues)))
        emptymodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second > 1.0 - threshold, spectrum |> geteigenvalues)))
        return Dict(:filled => columns(spectrum |> geteigenvectors, FockSpace(filledmodes)), :empty => columns(spectrum |> geteigenvectors, FockSpace(emptymodes)))
    end
    return frozenfockmaps
end
export frozenselectionbythreshold

function frozenselectionbycount(count::Integer)
    function frozenfockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
        𝑈ᵣ::FockMap = eigvecsh(𝐶ᵣ)
        modes::Vector{Mode} = [orderedmodes(𝑈ᵣ.inspace)...]
        return Dict(:filled => columns(𝑈ᵣ, FockSpace(Subset(modes[1:count]))), :empty => columns(𝑈ᵣ, FockSpace(Subset(modes[(end - count + 1):end]))))
    end
    return frozenfockmaps
end
export frozenselectionbycount

"""
    regioncorrelations(correlations::FockMap, regionfock::FockSpace)::FockMap

Restrict the crystal correlations to a real space regional local correlations.

### Input
- `correlations::FockMap`: The crystal correlations.
- `regionfock::FockSpace`: The real space regional `FockSpace`.

### Output
The real space local correlation with `inspace` and `outspace` as the `regionfock`.
"""
function regioncorrelations(correlations::FockMap, regionfock::FockSpace)::FockMap
    fouriermap::FockMap = fourier(correlations.inspace, regionfock) / (correlations.inspace |> subspacecount |> sqrt)
    return fouriermap' * correlations * fouriermap
end
export regioncorrelations

localfrozenisometries(
    correlations::FockMap, regionfock::FockSpace;
    selectionstrategy::Function = frozenselectionbythreshold(1e-3))::Dict{Symbol, FockMap} = (
    regioncorrelations(correlations, regionfock) |> selectionstrategy)
export localfrozenisometries

blocking(parameters::Pair{Symbol}...)::Dict{Symbol} = blocking(Dict(parameters...))
export blocking

function blocking(parameters::Dict{Symbol})::Dict{Symbol}
    @assert(haskey(parameters, :action))
    @assert(haskey(parameters, :correlations))
    @assert(haskey(parameters, :crystal))
    result::Dict{Symbol, Any} = Dict()

    scale::Scale = parameters[:action]
    correlation::FockMap = parameters[:correlations]
    crystal::Crystal = parameters[:crystal]

    scaling::FockMap = scale * correlation.inspace
    result[:action] = scaling
    result[:transformer] = scaling
    result[:correlations] = scaling * correlation * scaling'
    result[:crystal] = scale * crystal

    return result
end

function crystalisometries(; localisometry::FockMap, crystalfock::FockSpace{Crystal},
    addinspacemomentuminfo::Bool = false)::Dict{Momentum, FockMap}

    crystal::Crystal = getcrystal(crystalfock)
    fouriermap::FockMap = fourier(crystalfock, localisometry.outspace) / (crystal |> vol |> sqrt)
    momentumfouriers::Base.Generator = rowsubmaps(fouriermap)
    bz::Subset{Momentum} = brillouinzone(crystal)

    function preprocesslocalisometry(k::Momentum)::FockMap
        if !addinspacemomentuminfo
            return localisometry
        end
        inspace::FockSpace = localisometry.inspace |> orderedmodes |> setattr(:offset => k) |> FockSpace
        return FockMap(localisometry, inspace=inspace, performpermute=false)
    end

    return Dict(k => kfourier * preprocesslocalisometry(k) for (k, kfourier) in zip(bz, momentumfouriers))
end
export crystalisometries

function crystalisometry(; localisometry::FockMap, crystalfock::FockSpace{Crystal})::FockMap
    isometries::Dict{Momentum, FockMap} = crystalisometries(
        localisometry=localisometry, crystalfock=crystalfock, addinspacemomentuminfo=true)
    isometryunitcell::Subset{Offset} = Subset(mode |> getattr(:pos) for mode in localisometry.inspace |> orderedmodes)
    isometrycrystal::Crystal = Crystal(isometryunitcell, crystal.sizes)
    isometry::FockMap = (isometry for (_, isometry) in isometries) |> directsum
    isometrycrystalfock::CrystalFock = FockSpace(isometry.inspace, reflected=isometrycrystal)
    return FockMap(isometry, inspace=isometrycrystalfock, outspace=crystalfock)
end
export crystalisometry

function crystalprojector(; localisometry::FockMap, crystalfock::FockSpace{Crystal})
    momentumisometries::Dict{Point, FockMap} = crystalisometries(localisometry=localisometry, crystalfock=crystalfock)
    crystal::Crystal = getcrystal(crystalfock)
    bz::Subset{Momentum} = brillouinzone(crystal)
    globalprojector::FockMap = map(k -> momentumisometries[k] * momentumisometries[k]', bz) |> directsum
    return FockMap(
        globalprojector,
        outspace=FockSpace(globalprojector.outspace, reflected=crystal),
        inspace=FockSpace(globalprojector.inspace, reflected=crystal))
end
export crystalprojector

function crystalprojector(spectrum::CrystalSpectrum)::FockMap
    fockmap::FockMap = directsum(spectrum.eigenvectors[k] * spectrum.eigenvectors[k]' for (k, _) in spectrum.eigenmodes)
    basismodes::Subset{Mode} = spectrum |> unitcellfock |> orderedmodes
    fockspace::CrystalFock = crystalfock(basismodes, spectrum |> getcrystal)
    if hassamespan(fockmap |> getoutspace, fockspace)
        return FockMap(fockmap, inspace=fockspace, outspace=fockspace, performpermute=true)
    end
    # If some of the modes of the crystalfock is not included within the spectrum, we have to include them back into the nullspace of the projector.
    return zerosmap(fockspace, fockspace) + fockmap
end

function globaldistillerhamiltonian(;
    correlations::FockMap, restrictspace::FockSpace, localisometryselectionstrategy::Function, manualeigenenergies::Dict{Symbol, <:Number} = Dict(:filled => -1, :empty => 1),
    symmetry::AffineTransform)

    localisometries::Dict{Symbol} = localfrozenisometries(correlations, restrictspace, selectionstrategy=localisometryselectionstrategy)
    crystalprojectors::Dict{Symbol, FockMap} = Dict(
        name => crystalprojector(localisometry=localisometries[name], crystalfock=correlations.inspace)
        for (name, isometry) in localisometries)
    globaldistillhamiltonian::FockMap = reduce(+, manualeigenenergies[name] * crystalprojector for (name, crystalprojector) in crystalprojectors)

    order::Integer = symmetry |> pointgrouporder
    transformers::Vector{FockMap} = repeat([symmetry * globaldistillhamiltonian.outspace], order) |> cumprod

    return reduce(+, (transformer * globaldistillhamiltonian * transformer' for transformer in transformers))
end
export globaldistillerhamiltonian

function generategroupingfunction(grouppredicates)
    labels::Vector{Symbol} = [label for (label, _) in grouppredicates]
    predicates::Vector{Function} = [predicate for (_, predicate) in grouppredicates]
    function f(value::Number)::Symbol
        return labels[findfirst(p -> value |> p, predicates)]
    end
    return f
end

function distillation(spectrum::CrystalSpectrum, bandpredicates...)::Dict{Symbol, CrystalSpectrum}
    groupingfunction::Function = bandpredicates |> Renormalization.generategroupingfunction
    labeled::Base.Generator = ((mode |> getattr(:offset), v |> groupingfunction) => mode for (mode, v) in spectrum |> geteigenvalues)
    momentumgroups::Dict{Tuple, Vector{Mode}} = foldl(labeled; init=Dict{Tuple, Vector{Mode}}()) do d,(k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end

    bands::Dict{Symbol, Dict{Momentum, FockMap}} = Dict()
    for ((k, band), modes) in momentumgroups
        !haskey(bands, band) && (bands[band] = Dict())
        bands[band][k] = columns((spectrum |> geteigenvectors)[k], modes |> FockSpace)
    end

    function repacktospectrum(isometries::Dict{Momentum, FockMap})::CrystalSpectrum
        eigenmodes::Dict{Momentum, Subset{Mode}} = Dict(k => fockmap |> getinspace |> orderedmodes for (k, fockmap) in isometries)
        eigenvalues::Dict{Mode, Number} = Dict(mode => (spectrum |> geteigenvalues)[mode] for (_, modes) in eigenmodes for mode in modes)
        return CrystalSpectrum(spectrum.crystal, eigenmodes, eigenvalues, isometries)
    end

    return Dict(band => isometries |> repacktospectrum for (band, isometries) in bands)
end
export distillation

function findlocalspstates(;
    statecorrelations::FockMap, regionfock::FockSpace,
    symmetry::AffineTransform = identitytransform(statecorrelations |> getcrystal |> dimension),
    spectrumextractpredicate::Function = v -> v < 1e-2,
    linearindependencethreshold::Real = 5e-2,
    degeneracythreshold::Real = 1e-7)

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
    symmetricspstates = (state * getsymmetrizer(symmetry, state) for state in orthogonalspstates)
    spstates = (state * spatialmap(state)' for state in symmetricspstates)

    return Dict(state |> getinspace |> dimension => state for state in spstates)
end
export findlocalspstates

function wannierprojection(; crystalisometries::Dict{Momentum, FockMap}, crystal::Crystal, crystalseeds::Dict{Momentum, FockMap}, svdorthothreshold::Number = 1e-1)
    wannierunitcell::Subset{Offset} = Subset(mode |> getattr(:pos) for mode in (crystalseeds |> first |> last).inspace |> orderedmodes)
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

        return FockMap(approximated, inspace=FockSpace(approximated.inspace |> orderedmodes |> setattr(:offset => k)), performpermute=false)
    end
    if (precarioussvdvalues |> length) > 0
        @warn "Precarious wannier projection with minimum svdvalue of $(precarioussvdvalues |> minimum)"
    end
    wannierisometry::FockMap = directsum(approximateisometry(k, isometry, overlap) for (k, isometry, overlap) in overlaps)
    inspace::FockSpace = FockSpace(wannierisometry.inspace, reflected=wanniercrystal)
    outspace::FockSpace = FockSpace(wannierisometry.outspace, reflected=crystal)
    return FockMap(wannierisometry, inspace=inspace, outspace=outspace, performpermute=false)
end
export wannierprojection

struct RegionState{Dim} <: Element{FockMap}
    rep::FockMap
end
export RegionState

Base.:convert(::Type{FockMap}, source::RegionState)::FockMap = source.rep

function regionalrestriction(crystalstate::FockMap, regionfocks::FockSpace...)::RegionState
    positionmodes::Dict{Offset, Mode} = Dict(
        getattr(mode, :pos) => mode for mode in crystalstate |> getinspace |> unitcellfock |> orderedmodes)

    function regionstate(regionfock::FockSpace)::FockMap
        statecenter::Offset = Subset(mode |> pos for mode in regionfock |> orderedmodes) |> center
        haskey(positionmodes, statecenter) || error("Region center $statecenter is not a center in the crystal state!")
        mode::Mode = positionmodes[statecenter]
        rightfourier::FockMap = fourier(crystalstate |> getinspace, mode |> FockSpace)
        leftfourier::FockMap = fourier(crystalstate |> getoutspace, regionfock)
        return leftfourier' * crystalstate * rightfourier
    end

    localstate::FockMap = reduce(+, regionfock |> regionstate for regionfock in regionfocks)
    return localstate |> RegionState{crystalstate |> getoutspace |> getcrystal |> dimension}
end
export regionalrestriction

Quantum.:getinspace(state::RegionState) = state |> rep |> getinspace
Quantum.:getoutspace(state::RegionState) = state |> rep |> getoutspace

Base.:iterate(state::RegionState, i...) = iterate([columns(state |> rep, mode |> FockSpace) for mode in state |> getinspace |> orderedmodes], i...)
Base.:length(state::RegionState) = state |> getinspace |> dimension

end
