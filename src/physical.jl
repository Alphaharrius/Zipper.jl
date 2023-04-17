if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end

module Physical

using ..Spaces, ..Quantum, ..Geometries

function bondmap(bonds::Vector{Pair{Tuple{Mode, Mode}, ComplexF64}})::FockMap
    modes::Set{Mode} = Set([mode for bond in bonds for mode in bond.first])
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
        foreach(p -> push!(ret, setattr(p.first, :offset => momentum) => p.second), eigvalsh(bloch))
    end

    return ret
end

export Bond
export bloch, bondmap, espec

end
