if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end

module Zer

using ..Spaces, ..Geometries, ..Quantum

export DistillRegion
export frozenisometries, frozenselectionbythreshold, frozenselectionbycount

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

end
