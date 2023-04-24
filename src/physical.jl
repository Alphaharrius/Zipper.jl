if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end

module Physical

using OrderedCollections
using ..Spaces, ..Quantum, ..Geometries

function bondmap(bonds::Vector{Pair{Tuple{Mode, Mode}, ComplexF64}})::FockMap
    modes::OrderedSet{Mode} = OrderedSet([mode for bond in bonds for mode in bond.first])
    fockspace::FockSpace = FockSpace(Subset(modes))
    half::FockMap = FockMap(fockspace, fockspace, Dict(bonds))
    return half + half'
end

function espec(bonds::FockMap, momentums::Vector{Point})::Vector{Pair{Mode, Float64}}
    ret::Vector{Pair{Mode, Float64}} = []
    bondmodes::Subset{Mode} = flatten(rep(bonds.inspace))

    for momentum in momentums
        Fₖ::FockMap = fourier(Subset([momentum]), bondmodes)
        bloch::FockMap = Fₖ * bonds * Fₖ'
        foreach(p -> push!(ret, p), eigvalsh(bloch, :offset => momentum))
    end

    return ret
end

function hamiltonian(crystal::Crystal, bondmap::FockMap)::FockMap
    𝐵𝑍::Subset{Point} = brillouin_zone(crystal)
    bondmodes::Subset{Mode} = flatten(rep(bondmap.outspace))
    sqrtvol::Float64 = sqrt(vol(crystal))
    ∑𝐹ₖ = Iterators.map(𝑘 -> fourier(Subset([𝑘]), bondmodes) / sqrtvol, 𝐵𝑍)
    ∑𝐻ₖ = Iterators.map(𝐹ₖ -> 𝐹ₖ * bondmap * 𝐹ₖ', ∑𝐹ₖ)
    return directsum([∑𝐻ₖ...])
end

function filledstates(hamiltonian::FockMap)::FockMap
    bloch_partitions::OrderedSet{Subset{Mode}} = partitions(hamiltonian.inspace)
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
    return directsum(∑𝑈₀)
end

function ground_state_correlation(hamiltonian::FockMap)::FockMap
    spectrum::Vector{Pair{Mode, Float64}} = @time eigvalsh(hamiltonian) # Since the eigenmodes here is intermediate within this function, we don't have to identify them.
    contributing_modes::Vector{Mode} = @time map(p -> p.first, filter(p -> p.second <= 0, spectrum)) # Extract all eigenmodes that has negative eigenenergy.
    𝖀::FockMap = @time eigvecsh(hamiltonian) # Unitary.
    𝑈₀::FockMap = @time columns(𝖀, FockSpace(Subset(contributing_modes))) # Extract the eigenmode representations with negative eigenenergies.
    return @time 𝑈₀ * 𝑈₀' # Using convention of 𝐶 ≔ ⟨𝑐𝑐†⟩
end

export Bond
export bloch, bondmap, espec, hamiltonian, filledstates, ground_state_correlation

end
