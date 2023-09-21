module Renormalization

using LinearAlgebra, OrderedCollections
using ..Spaces, ..Geometries, ..Quantum, ..Transformations

function frozenselectionbythreshold(threshold::Float64)
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

localfrozenisometries(
    correlations::FockMap, regionfock::FockSpace;
    selectionstrategy = frozenselectionbythreshold(1e-3))::Dict{Symbol, FockMap} = (
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

    crystal::Crystal = crystalof(crystalfock)
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

function crystalprojector(; localisometry::FockMap, crystalfock::FockSpace{Crystal})
    momentumisometries::Dict{Point, FockMap} = crystalisometries(localisometry=localisometry, crystalfock=crystalfock)
    crystal::Crystal = crystalof(crystalfock)
    bz::Subset{Momentum} = brillouinzone(crystal)
    globalprojector::FockMap = map(k -> momentumisometries[k] * momentumisometries[k]', bz) |> directsum
    return FockMap(
        globalprojector,
        outspace=FockSpace(globalprojector.outspace, reflected=crystal),
        inspace=FockSpace(globalprojector.inspace, reflected=crystal))
end
export crystalprojector

function globaldistillerhamiltonian(;
    correlations::FockMap, restrictspace::FockSpace, localisometryselectionstrategy, manualeigenenergies::Dict{Symbol, <:Number} = Dict(:filled => -1, :empty => 1),
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

function distillation(globaldistiller::FockMap; courierenergythreshold::Number = 1e-5)
    spectrum::CrystalSpectrum = globaldistiller |> crystalspectrum
    partitionrule = v -> v > courierenergythreshold ? :empty : v < -courierenergythreshold ? :filled : :courier
    labelled::Base.Generator = ((mode |> getattr(:offset), v |> partitionrule) => mode for (mode, v) in spectrum.eigenvalues)
    isometrygroups::Dict{Tuple, Vector{Mode}} = foldl(labelled; init=Dict{Tuple, Vector{Mode}}()) do d,(k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end
    filledcrystalisometries::Dict{Momentum, FockMap} = Dict()
    emptycrystalisometries::Dict{Momentum, FockMap} = Dict()
    couriercrystalisometries::Dict{Momentum, FockMap} = Dict()
    for (k, unitary) in spectrum.eigenvectors
        filledfock::FockSpace = isometrygroups[(k, :filled)] |> Subset |> FockSpace
        filledcrystalisometries[k] = columns(unitary, filledfock)
        emptyfock::FockSpace = isometrygroups[(k, :empty)] |> Subset |> FockSpace
        emptycrystalisometries[k] = columns(unitary, emptyfock)
        courierfock::FockSpace = isometrygroups[(k, :courier)] |> Subset |> FockSpace
        couriercrystalisometries[k] = columns(unitary, courierfock)
    end

    function repacktospectrum(isometries::Dict{Momentum, FockMap})::CrystalSpectrum
        eigenmodes::Dict{Momentum, Subset{Mode}} = Dict(k => fockmap.inspace |> orderedmodes for (k, fockmap) in isometries)
        eigenvalues::Dict{Mode, Number} = Dict(mode => spectrum.eigenvalues[mode] for (_, modes) in eigenmodes for mode in modes)
        return CrystalSpectrum(spectrum.crystal, eigenmodes, eigenvalues, isometries)
    end

    return Dict(
        :empty => emptycrystalisometries |> repacktospectrum,
        :filled => filledcrystalisometries |> repacktospectrum,
        :courier => couriercrystalisometries |> repacktospectrum)
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

        dimension::Integer = correlations.inspace |> crystalof |> dimension
        seedfock::FockSpace = Subset(
            mode |> setattr(:orbital => findeigenfunction(symmetry; dimensionrange=0:dimension, eigenvalue=phase))
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

end
