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
    return half + dagger(half)
end

function espec(bonds::FockMap, momentums::Vector{Point})::Vector{Pair{Mode, Float64}}
    ret::Vector{Pair{Mode, Float64}} = []
    bondmodes::Subset{Mode} = flatten(rep(bonds.inspace))

    for momentum in momentums
        Fₖ::FockMap = fourier(Subset([momentum]), bondmodes)
        bloch::FockMap = Fₖ * bonds * dagger(Fₖ)
        foreach(p -> push!(ret, p), eigvalsh(bloch, :offset => momentum))
    end

    return ret
end

function hamiltonian(crystal::Crystal, bondmap::FockMap)::FockMap
    bzone::Subset{Point} = brillouin_zone(crystal)
    bondmodes::Subset{Mode} = flatten(rep(bondmap.outspace))
    fmaps = Iterators.map(momentum -> fourier(Subset([momentum]), bondmodes), bzone)
    bloch_hamiltonians = Iterators.map(fmap -> fmap * bondmap * dagger(fmap), fmaps)
    return directsum([bloch_hamiltonians...])
end

function ground_state_correlation(hamiltonian::FockMap)::FockMap
    spectrum::Vector{Pair{Mode, Float64}} = eigvalsh(hamiltonian) # Since the eigenmodes here is intermediate within this function, we don't have to identify them.
    contributing_modes::Vector{Mode} = map(p -> p.first, filter(p -> p.second <= 0, spectrum)) # Extract all eigenmodes that has negative eigenenergy.
    unitary::FockMap = eigvecsh(hamiltonian)
    groundstate_unitary::FockMap = columns(unitary, Subset(contributing_modes)) # Extract the eigenmode representations with negative eigenenergies.
    return groundstate_unitary * dagger(groundstate_unitary) # Using convention of C = <cc†>
end

export Bond
export bloch, bondmap, espec, hamiltonian, ground_state_correlation

end
