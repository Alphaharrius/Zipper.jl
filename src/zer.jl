module Zer

using ..Spaces, ..Geometries, ..Quantum, ..Transformations
export DistillRegion
export localfrozenisometries, frozenselectionbythreshold, frozenselectionbycount
export blocking

function frozenselectionbythreshold(threshold::Float64)
    function frozenfockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
        spectrum, 𝑈ᵣ::FockMap = eigh(𝐶ᵣ)
        filledmodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second < threshold, spectrum)))
        emptymodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second > 1.0 - threshold, spectrum)))
        return Dict(:filled => columns(𝑈ᵣ, FockSpace(filledmodes)), :empty => columns(𝑈ᵣ, FockSpace(emptymodes)))
    end
    return frozenfockmaps
end

function frozenselectionbycount(count::Integer)
    function frozenfockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
        𝑈ᵣ::FockMap = eigvecsh(𝐶ᵣ)
        modes::Vector{Mode} = [orderedmodes(𝑈ᵣ.inspace)...]
        return Dict(:filled => columns(𝑈ᵣ, FockSpace(Subset(modes[1:count]))), :empty => columns(𝑈ᵣ, FockSpace(Subset(modes[(end - count + 1):end]))))
    end
    return frozenfockmaps
end

function localfrozenisometries(correlations::FockMap, restrictedfock::FockSpace;
                               selectionstrategy = frozenselectionbythreshold(1e-3))::Dict{Symbol, FockMap}
    fouriermap::FockMap = fourier(correlations.inspace, restrictedfock) / (correlations.inspace |> subspacecount |> sqrt)
    restrictedcorrelations::FockMap = fouriermap' * correlations * fouriermap
    return selectionstrategy(restrictedcorrelations)
end

blocking(parameters::Pair{Symbol}...)::Dict{Symbol} = blocking(Dict(parameters...))

function blocking(parameters::Dict{Symbol})::Dict{Symbol}
    @assert(haskey(parameters, :scale))
    @assert(haskey(parameters, :correlations))
    @assert(haskey(parameters, :crystal))
    result::Dict{Symbol, Any} = Dict()

    scale::Scale = parameters[:scale]
    correlation::FockMap = parameters[:correlations]
    crystal::Crystal = parameters[:crystal]

    scaling::FockMap = scale * correlation.inspace
    result[:action] = scaling
    result[:correlations] = scaling * correlation * scaling'
    result[:crystal] = scale * crystal

    return result
end

function fourierisometries(; localisometry::FockMap, crystalfock::FockSpace{Crystal})::Dict{Point, FockMap}
    crystal::Crystal = crystalof(crystalfock)
    fouriermap::FockMap = fourier(crystalfock, localisometry.outspace) / (crystal |> vol |> sqrt)
    momentumfouriers::Vector{FockMap} = colsubmaps(fouriermap)
    bz::Subset{Point} = brillouinzone(crystal)
    return Dict(k => fourier_k * localisometry for (k, fourier_k) in zip(bz, momentumfouriers))
end

function isometryglobalprojector(; localisometry::FockMap, crystalfock::FockSpace{Crystal})
    momentuisometries::Dict{Point, FockMap} = fourierisometries(localisometry=localisometry, crystalfock=crystalfock)
    crystal::Crystal = crystalof(crystalfock)
    bz::Subset{Point} = brillouinzone(crystal)
    globalprojector::FockMap = focksum(map(k -> momentuisometries[k] * momentuisometries[k]', bz))
    return FockMap(
        globalprojector,
        outspace=FockSpace(globalprojector.outspace, reflected=crystal),
        inspace=FockSpace(globalprojector.inspace, reflected=crystal))
end

end
