if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end

module Zer

using ..Spaces, ..Geometries, ..Quantum

export DistillRegion
export frozenisometries, frozenselectionbythreshold, frozenselectionbycount

struct DistillRegion
    center::Point
    modes::Subset{Mode}
end

function frozenselectionbythreshold(threshold::Float64)
    function frozenfocks(𝐶ᵣ::FockMap)::Dict{Symbol, FockSpace}
        corrspec::Vector{Pair{Mode, Float64}} = eigvalsh(𝐶ᵣ)
        filledmodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second < threshold, corrspec)))
        emptymodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second > 1.0 - threshold, corrspec)))
        return Dict(:filled => FockSpace(filledmodes), :empty => FockSpace(emptymodes))
    end
    return frozenfocks
end

function frozenselectionbycount(count::Integer)
    function frozenfocks(𝐶ᵣ::FockMap)::Dict{Symbol, FockSpace}
        𝑈ᵣ::FockMap = eigvecsh(𝐶ᵣ)
        modes::Subset{Mode} = orderedmodes(𝑈ᵣ.inspace)
        return Dict(:filled => FockSpace(modes[1:count]), :empty => FockSpace(modes[(end - count):end]))
    end
    return frozenfocks
end

function frozenisometries(crystal::Crystal, correlations::FockMap, region::DistillRegion;
                          selectionstrategy = frozenselectionbythreshold(1e-3))::Dict{Symbol, FockMap}
    𝐹ₖ::FockMap = fourier(brillouin_zone(crystal), FockSpace(region.modes)) / sqrt(vol(crystal))
    𝐶ᵣ::FockMap = 𝐹ₖ' * correlations * 𝐹ₖ
    frozenfocks::Dict{Symbol, FockSpace} = selectionstrategy(𝐶ᵣ)
    return Dict(:filled => columns(𝐶ᵣ, frozenfocks[:filled]), :empty => columns(𝐶ᵣ, frozenfocks[:empty]))
end

end
