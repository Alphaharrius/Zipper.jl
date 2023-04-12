if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end

module Physical

using ..Spaces, ..Quantum, ..Geometries

struct Bond
    modes::Tuple{Mode, Mode}
    offset::Point
    strength::Number
end

function bloch(bonds::Set{Bond}, k::Point, crystal::Crystal)::FockMap
    all_modes::Set{Mode} = Set([mode for bond in bonds for mode in bond.modes])
    fock_space::FockSpace = FockSpace(Subset(all_modes))
    fock_dict::Dict{Tuple{Mode, Mode}, ComplexF64} = Dict()
    foreach(bond -> fock_dict[bond.modes] = get(fock_dict, bond.modes, 0im) + fourier_coef(k, bond.offset, vol(crystal)) * bond.strength, bonds)
    upper_triangular::FockMap = FockMap(fock_space, fock_space, fock_dict)
    return upper_triangular + dagger(upper_triangular)
end

export Bond
export bloch, test

end
