module Physical

using OrderedCollections
using ..Spaces, ..Quantum, ..Geometries


"""
    bondmap(bonds::Vector{Pair{Tuple{Mode, Mode}, ComplexF64}})::FockMap

Generate a bonding information `FockMap` from the given `bonds` between real space `Mode` objects, the `inspace` & `outspace`
of the `FockMap` will be spanned from all provided `Mode` objects given in `bonds`.

### Input

- `bonds`   A vector of pairs of `Mode` objects `(frommode, tomode)` and their corresponding complex hopping amplitudes.
- `crystal` The `Crystal` object to compute the energy spectrum on.

### Output
The `FockMap` object that contains the bonding information.

### Example

```julia-repl
# Create the bonding information for a honeycomb lattice based Dirac semi-metal.
tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (modeA, modeB) => tₙ,
    (modeA, modeB |> setattr(:offset => triangularspace & [-1, 0])) => tₙ,
    (modeA, modeB |> setattr(:offset => triangularspace & [0, 1])) => tₙ])
```
"""
function bondmap(bonds::Vector{Pair{Tuple{Mode, Mode}, ComplexF64}})::FockMap
    fockspace::FockSpace = FockSpace(Subset(mode for bond in bonds for mode in bond.first))
    half::FockMap = FockMap(fockspace, fockspace, Dict(bonds))
    diag::FockMap = FockMap(fockspace, fockspace, Dict(filter(p -> p.first |> first == p.first |> last, bonds)))
    return half + half' - diag
end
export bondmap

"""
    computeenergyspectrum(bonds::FockMap; crystal::Crystal)::CrystalSpectrum

Compute the momentum space energy spectrum within a `Crystal` brillouin zone for a given `bonds` generated from `bondmap`.
"""
function computeenergyspectrum(bonds::FockMap; crystal::Crystal)::CrystalSpectrum
    bz::Subset{Momentum} = crystal |> brillouinzone
    bondmodes::Subset{Mode} = bonds.outspace |> orderedmodes
    fouriers::Base.Generator = (k => fourier(k |> Subset, bondmodes |> FockSpace) for k in bz)
    momentumhamiltonians::Base.Generator = (k => fourier * bonds * fourier' for (k, fourier) in fouriers)
    return crystalspectrum(momentumhamiltonians, crystal=crystal)
end
export computeenergyspectrum

"""
    groundstatespectrum(energyspectrum::CrystalSpectrum; perunitcellfillings::Integer, energyresolution::Real = 1e-7)::CrystalSpectrum

Filter out the high energy degrees of freedom in the crystal energy spectrum, determined by the filling information, and return a new
spectrum by filling all fermionic degrees of freedom below the Fermi energy.

### Input
- `energyspectrum`      The `CrystalSpectrum` to corresponds to the crystal Hamiltonian.
- `perunitcellfillings` The number of fermionic degrees of freedom to fill per unit cell.
- `energyresolution`    The energy resolution to use when grouping the energy spectrum degeneracies, if this value is set too high it
                        might group bands that are not supposed to be degenerated.

### Output
The `CrystalSpectrum` object that contains the ground state energy spectrum.
"""
function groundstatespectrum(energyspectrum::CrystalSpectrum; perunitcellfillings::Integer, energyresolution::Real = 1e-7)::CrystalSpectrum
    perunitcellfillings > 0 || error("Filling must be positive!")
    
    unitcellvacancies::Integer = energyspectrum |> geteigenvectors |> values |> first |> getinspace |> dimension
    perunitcellfillings <= unitcellvacancies || error("Unable to fill more than $unitcellvacancies per unit cell!")
    totalfillings::Integer = perunitcellfillings * (energyspectrum |> getcrystal |> vol)

    groupedeigenvalues::Base.Generator = groupbyeigenvalues(energyspectrum, groupingthreshold=energyresolution)
    cumfillings::Vector = cumsum(modeset |> length for (_, modeset) in groupedeigenvalues)
    groupcollectcount::Integer = findfirst(v -> v >= totalfillings, cumfillings)
    contributions = Iterators.take(groupedeigenvalues, groupcollectcount)
    groundstatemodesets::Base.Generator = (modeset for (_, modeset) in contributions)
    groundstatemodes::Subset{Mode} = reduce(+, groundstatemodesets)

    decoratedmodes::Base.Generator = ((m |> getattr(:offset)) => m for m in groundstatemodes)
    momentummodes::Dict = foldl(decoratedmodes; init=Dict()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end

    crystaleigenmodes::Dict{Momentum, Subset{Mode}} = Dict(k => modes |> Subset for (k, modes) in momentummodes)
    crystaleigenvalues::Dict{Mode, Number} = Dict(m => (energyspectrum |> geteigenvalues)[m] for m in groundstatemodes)
    crystaleigenvectors::Dict{Momentum, FockMap} = Dict(
        k => columns((energyspectrum |> geteigenvectors)[k], crystaleigenmodes[k] |> FockSpace) for k in crystaleigenmodes |> keys)
    return CrystalSpectrum(energyspectrum |> getcrystal, crystaleigenmodes, crystaleigenvalues, crystaleigenvectors)
end
export groundstatespectrum

"""
    roundingpurification(correlationspectrum::CrystalSpectrum)::CrystalSpectrum

Round the correlation eigenvalues of the `correlationspectrum` to `0` or `1`, this is
useful when the eigenvalues in the pure state are close to `0` or `1` but not exactly.
"""
function roundingpurification(correlationspectrum::CrystalSpectrum)::CrystalSpectrum
    function fixeigenvalues(eigenvalue::Number)
        if isapprox(eigenvalue, 1, atol=1e-1)
            return 1
        end
        if isapprox(eigenvalue, 0, atol=1e-1)
            return 0
        end
        @warn("Eigenvalue $(eigenvalue) is not close to 0 or 1!")
    end
    eigenvalues::Dict{Mode, Number} = Dict(m => v |> fixeigenvalues for (m, v) in correlationspectrum.eigenvalues)
    return CrystalSpectrum(correlationspectrum.crystal, correlationspectrum.eigenmodes, eigenvalues, correlationspectrum.eigenvectors)
end
export roundingpurification

end
