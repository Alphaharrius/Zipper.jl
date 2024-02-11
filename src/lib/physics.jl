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
    (modeA, modeB |> setattr(:r => [-1, 0] ∈ triangularspace)) => tₙ,
    (modeA, modeB |> setattr(:r => [0, 1] ∈ triangularspace)) => tₙ])
```
"""
function bondmap(bonds::Vector{Pair{Tuple{Mode, Mode}, T}})::FockMap where {T <: Complex}
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
    bondmodes::Subset{Mode} = bonds|>getoutspace|>orderedmodes
    basismodes::Subset{Mode} = bonds|>getoutspace|>RegionFock|>unitcellfock|>orderedmodes
    transform::FockMap = fourier(getcrystalfock(basismodes, crystal), bondmodes|>RegionFock)
    function compute(k)
        subspace = getsubspace(transform|>getoutspace, k)
        ktransform = transform[subspace, :]
        return k=>(ktransform*bonds*ktransform')
    end
    momentumhamiltonians = paralleltasks(
        name="computeenergyspectrum",
        tasks=(()->compute(k) for k in crystal|>brillouinzone),
        count=crystal|>vol)|>parallel
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
# TODO: This function requires revisit as the performance is bad, but the results are correct.
function groundstatespectrum(energyspectrum::CrystalSpectrum; perunitcellfillings::Number, energyresolution::Real = 1e-7)::CrystalSpectrum
    perunitcellfillings > 0 || error("Filling must be positive!")
    
    unitcellvacancies::Integer = energyspectrum |> geteigenvectors |> values |> first |> getinspace |> dimension
    perunitcellfillings <= unitcellvacancies || error("Unable to fill more than $unitcellvacancies per unit cell!")
    totalfillings::Integer = (perunitcellfillings * (energyspectrum |> getcrystal |> vol)) |> round |> Int

    groupedeigenvalues::Base.Generator = groupbyeigenvalues(energyspectrum, groupingthreshold=energyresolution)
    cumfillings::Vector = cumsum(modeset |> length for (_, modeset) in groupedeigenvalues)
    groupcollectcount::Integer = findfirst(v -> v >= totalfillings, cumfillings)
    contributions = Iterators.take(groupedeigenvalues, groupcollectcount)
    groundstatemodesets::Base.Generator = (modeset for (_, modeset) in contributions)
    groundstatemodes::Subset{Mode} = groundstatemodesets |> subsetunion

    decoratedmodes::Base.Generator = ((m |> getattr(:k)) => m for m in groundstatemodes)
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

function entanglemententropy(correlationspectrum::CrystalSpectrum)
    sumresult = 0
    for (_, v) in correlationspectrum |> geteigenvalues
        if isapprox(v, 0, atol=1e-7) || isapprox(v, 1, atol=1e-7)
            continue
        end
        
        sumresult += v * log(v) + (1 - v) * log(1 - v)
    end
    return -sumresult
end
export entanglemententropy

function momentumoccupations(correlations::FockMap)::FockMap
    kcorrelations::Base.Generator = correlations |> crystalsubmaps
    crystal::Crystal = correlations |> getoutspace |> getcrystal
    center::Offset = crystal |> getspace |> getorigin
    tracecrystal::Crystal = Crystal(center |> Subset, crystal |> size)
    mode::Mode = Mode(:b => center)

    function tracing(k::Momentum, corr::FockMap)::FockMap
        space::FockSpace = mode |> setattr(:k => k) |> FockSpace
        return FockMap(space, space, [corr |> tr][:, :] |> SparseMatrixCSC) / dimension(corr |> getoutspace)
    end

    occupations::FockMap = directsum(tracing(k, corr) for (k, corr) in kcorrelations)
    fockspace::FockSpace = FockSpace(occupations |> getoutspace, reflected=tracecrystal)
    return FockMap(occupations, inspace=fockspace, outspace=fockspace)
end
export momentumoccupations
