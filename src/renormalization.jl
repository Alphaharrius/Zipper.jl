module Renormalization

using LinearAlgebra, OrderedCollections
using ..Spaces, ..Geometries, ..Quantum, ..Transformations

function frozenselectionbythreshold(threshold::Float64)::Function
    function frozenfockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
        spectrum, 𝑈ᵣ::FockMap = eigh(𝐶ᵣ)
        filledmodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second < threshold, spectrum)))
        emptymodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second > 1.0 - threshold, spectrum)))
        return Dict(:filled => columns(𝑈ᵣ, FockSpace(filledmodes)), :empty => columns(𝑈ᵣ, FockSpace(emptymodes)))
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
    isometryunitcell::Subset{Position} = Subset(mode |> getattr(:pos) for mode in localisometry.inspace |> orderedmodes)
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
    crystalfock::FockSpace = FockSpace(fockmap.inspace, reflected=spectrum.crystal)
    return FockMap(fockmap, inspace=crystalfock, outspace=crystalfock, performpermute=false)
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

function regionalwannierseeding(statecorrelations::FockMap, regionspace::FockSpace;
    symmetry::AffineTransform,
    seedingthreshold::Number = 1e-2, seedsgroupingprecision::Number = 1e-5, linearindependencethreshold::Number = 5e-2)

    localcorrelations::FockMap = regioncorrelations(statecorrelations, regionspace)
    localspectrum::EigenSpectrum = localcorrelations |> eigspech
    validgroups = Iterators.filter(p -> p.first <= seedingthreshold, groupbyeigenvalues(localspectrum, groupingthreshold=seedsgroupingprecision))

    function symmetrizeseed(seedisometry::FockMap)::Tuple{Base.Generator, FockMap}
        transform::FockMap = symmetry * seedisometry.outspace
        eigensymmetryrep::FockMap = seedisometry' * transform * seedisometry
        eigenvalues, unitary = eigensymmetryrep |> eigen
        return (eigenvalues, seedisometry * unitary)
    end

    regioncenter::Point = Subset(mode |> pos for mode in regionspace |> getmodes) |> center

    function extractglobalseed(group::Pair{<:Number, Subset{Mode}})
        phases, seed = columns(localspectrum.eigenvectors, group.second |> FockSpace) |> symmetrizeseed
        crystalseeds::Dict{Momentum, FockMap} = crystalisometries(localisometry=seed, crystalfock=statecorrelations.inspace)
        crystalseed::FockMap = directsum(v for (_, v) in crystalseeds)

        # Check if linear independent.
        pseudoidentity::FockMap = (crystalseed' * crystalseed)
        mineigenvalue = min(map(p -> p.second, pseudoidentity |> eigvalsh)...)
        if mineigenvalue < linearindependencethreshold
            return nothing
        end

        dim::Integer = statecorrelations.inspace |> getcrystal |> dimension
        seedfock::FockSpace = Subset(
            mode |> setattr(:orbital => findeigenfunction(symmetry; dimensionrange=0:dim, eigenvalue=phase))
                 |> setattr(:pos => regioncenter)
                 |> setattr(:offset => regioncenter |> getspace |> origin)
            for (mode, phase) in phases) |> FockSpace

        return FockMap(seed; inspace=seedfock, performpermute=false)
    end

    return Iterators.filter(v -> !(v isa Nothing), extractglobalseed(group) for group in validgroups)
end
export regionalwannierseeding

function wannierprojection(; crystalisometries::Dict{Momentum, FockMap}, crystal::Crystal, crystalseeds::Dict{Momentum, FockMap}, svdorthothreshold::Number = 1e-1)
    wannierunitcell::Subset{Position} = Subset(mode |> getattr(:pos) for mode in (crystalseeds |> first |> last).inspace |> orderedmodes)
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
    positionmodes::Dict{Position, Mode} = Dict(
        getattr(mode, :pos) => mode for mode in crystalstate |> getinspace |> unitcellfock |> orderedmodes)

    function regionstate(regionfock::FockSpace)::FockMap
        statecenter::Position = Subset(mode |> pos for mode in regionfock |> orderedmodes) |> center
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
