module Zer

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
    symmetries::Vector{PointGroupTransformation})

    localisometries::Dict{Symbol} = localfrozenisometries(correlations, restrictspace, selectionstrategy=localisometryselectionstrategy)
    crystalprojectors::Dict{Symbol, FockMap} = Dict(
        name => crystalprojector(localisometry=localisometries[name], crystalfock=correlations.inspace)
        for (name, isometry) in localisometries)
    globaldistillhamiltonian::FockMap = reduce(+, manualeigenenergies[name] * crystalprojector for (name, crystalprojector) in crystalprojectors)

    transformations = (transformation for symmetry in symmetries for transformation in symmetry |> pointgroupelements)

    function transformdistillhamiltonian(transformation::PointGroupTransformation)::FockMap
        transformer::FockMap = transformation * globaldistillhamiltonian.outspace
        return transformer * globaldistillhamiltonian * transformer'
    end

    return reduce(+, (transformation |> transformdistillhamiltonian for transformation in transformations))
end
export globaldistillerhamiltonian

end
