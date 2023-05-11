module Zer

using ..Spaces, ..Geometries, ..Quantum, ..Transformations
export DistillRegion
export frozenisometries, frozenselectionbythreshold, frozenselectionbycount
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

function frozenisometries(correlations::FockMap, restrictedfock::FockSpace;
                          selectionstrategy = frozenselectionbythreshold(1e-3))::Dict{Symbol, FockMap}
    𝐹ₖ::FockMap = fourier(correlations.inspace, restrictedfock) / sqrt(subspacecount(correlations.inspace))
    𝐶ᵣ::FockMap = 𝐹ₖ' * correlations * 𝐹ₖ
    return selectionstrategy(𝐶ᵣ)
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

    scaling::FockMap = scale * Recipient(correlation.inspace, :crystal => crystal)
    result[:action] = scaling
    result[:correlations] = scaling * correlation * scaling'
    result[:crystal] = scale * crystal

    return result
end

end
