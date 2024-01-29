"""
    modeselectionbythreshold(threshold::Float64)

A selection strategy for `localisometries` that selects the modes with correlation eigenvalues γ
that is localized within the local region ∀ γ ≈ 0 or γ ≈ 1 up to a `threshold`. Extension of "frozenselectionbythreshold" in renormalization by including courier modes

### Input
- `threshold::Float64`: The threshold for the correlation eigenvalues.

### Output
A function representing the selection strategy with an input of the local correlations 𝐶ᵣ `FockMap` and outputs a `Dict{Symbol, FockMap}`
keyed by three grouping symbols with their associated local isometries selected from the unitary that diagonalizes 𝐶ᵣ. The grouping symbol
`:filled` represents the eigenmodes that corresponds to γ ≈ 0; `:empty` represents the eigenmodes that corresponds to γ ≈ 1; 
`:courier` represents the eigenmodes that corresponds to 0 ≈< γ <≈ 1
"""
function modeselectionbythreshold(threshold::Float64)::Function
    function modesfockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
        spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
        filledmodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second < threshold, spectrum |> geteigenvalues)))
        emptymodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second > 1.0 - threshold, spectrum |> geteigenvalues)))
        couriermodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> threshold <= p.second <= 1.0 - threshold, spectrum |> geteigenvalues)))
        return Dict(:filled => columns(spectrum |> geteigenvectors, FockSpace(filledmodes)), :empty => columns(spectrum |> geteigenvectors, FockSpace(emptymodes)),
        :courier => columns(spectrum |> geteigenvectors, FockSpace(couriermodes)))
    end
    return modesfockmaps
end
export modeselectionbythreshold

function modeselectionbycount(count::Integer)
    function modefockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
        𝑈ᵣ::FockMap = eigvecsh(𝐶ᵣ)
        modes::Vector{Mode} = [orderedmodes(𝑈ᵣ.inspace)...]
        return Dict(:filled => columns(𝑈ᵣ, FockSpace(Subset(modes[1:count]))), :empty => columns(𝑈ᵣ, FockSpace(Subset(modes[(end - count + 1):end]))),
        :courier => columns(𝑈ᵣ, FockSpace(Subset(modes[count+1:(end - count)]))))
    end
    return modefockmaps
end
export modeselectionbycount

localisometries(
    correlations::FockMap, regionfock::FockSpace;
    selectionstrategy::Function = modeselectionbythreshold(1e-3))::Dict{Symbol, FockMap} = (
    regioncorrelations(correlations, regionfock) |> selectionstrategy)
export localisometries
