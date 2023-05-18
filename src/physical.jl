if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end

module Physical

using OrderedCollections
using ..Spaces, ..Quantum, ..Geometries

export bloch, bondmap, espec, hamiltonian, groundstates, groundstatecorrelation

function bondmap(bonds::Vector{Pair{Tuple{Mode, Mode}, ComplexF64}})::FockMap
    fockspace::FockSpace = FockSpace(Subset(mode for bond in bonds for mode in bond.first))
    half::FockMap = FockMap(fockspace, fockspace, Dict(bonds))
    return half + half'
end

function espec(bonds::FockMap, momentums::Vector{Point})::Vector{Pair{Mode, Float64}}
    ret::Vector{Pair{Mode, Float64}} = []
    bondmodes::Subset{Mode} = flatten(rep(bonds.inspace))

    for momentum in momentums
        Fₖ::FockMap = fourier(Subset(momentum), bondmodes)
        bloch::FockMap = Fₖ * bonds * Fₖ'
        foreach(p -> push!(ret, p), eigvalsh(bloch, :offset => momentum))
    end

    return ret
end

function hamiltonian(crystal::Crystal, bondmap::FockMap)::FockMap
    𝐵𝑍::Subset{Point} = brillouinzone(crystal)
    bondmodes::Subset{Mode} = flatten(rep(bondmap.outspace))
    ∑𝐹ₖ = Iterators.map(𝑘 -> fourier(Subset(𝑘), FockSpace(bondmodes)), 𝐵𝑍)
    ∑𝐻ₖ = Iterators.map(𝐹ₖ -> 𝐹ₖ * bondmap * 𝐹ₖ', ∑𝐹ₖ)
    fockmap::FockMap = focksum([∑𝐻ₖ...])
    return FockMap(FockSpace(fockmap.outspace, T=CrystalFock), FockSpace(fockmap.inspace, T=CrystalFock), rep(fockmap))
end

function groundstates(hamiltonian::FockMap)::FockMap
    bloch_partitions::Subset{Subset{Mode}} = rep(hamiltonian.inspace)
    ∑𝑈₀::Vector{FockMap} = Vector{FockMap}(undef, length(bloch_partitions))
    filledgroup::ModeGroup = ModeGroup(quantized, "filled")
    for (n, partition) in enumerate(bloch_partitions)
        𝑘::Point = getattr(first(partition), :offset)
        fockspace::FockSpace = FockSpace(partition)
        𝐻ₖ::FockMap = restrict(hamiltonian, fockspace, fockspace)
        𝔘::FockMap = eigvecsh(𝐻ₖ, :offset => 𝑘, :groups => [filledgroup])
        filledmodes::Vector{Mode} = map(p -> p.first, filter(p -> p.second < 0, eigvalsh(𝐻ₖ, :offset => 𝑘, :groups => [filledgroup])))
        𝑈₀::FockMap = columns(𝔘, FockSpace(Subset(filledmodes)))
        ∑𝑈₀[n] = 𝑈₀
    end
    𝔘₀::FockMap = focksum(∑𝑈₀)
    return FockMap(FockSpace(𝔘₀.outspace, T=CrystalFock), 𝔘₀.inspace, rep(𝔘₀))
end

function groundstatecorrelation(hamiltonian::FockMap)::FockMap
    𝔘₀::FockMap = groundstates(hamiltonian)
    𝐼::FockMap = idmap(𝔘₀.outspace, 𝔘₀.outspace)
    𝐶::FockMap = 𝐼 - 𝔘₀ * 𝔘₀'
    return FockMap(FockSpace(𝐶.outspace, T=CrystalFock), FockSpace(𝐶.inspace, T=CrystalFock), rep(𝐶))
end

end
