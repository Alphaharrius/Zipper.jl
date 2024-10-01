"""
    frozenselectionbythreshold(threshold::Float64)

A selection strategy for `localfrozenisometries` and `globaldistillerhamiltonian` that selects the modes with correlation eigenvalues γ
that is localized within the local region ∀ γ ≈ 0 or γ ≈ 1 up to a `threshold`.

### Input
- `threshold::Float64`: The threshold for the correlation eigenvalues.

### Output
A function representing the selection strategy with an input of the local correlations 𝐶ᵣ `FockMap` and outputs a `Dict{Symbol, FockMap}`
keyed by two grouping symbols with their associated local isometries selected from the unitary that diagonalizes 𝐶ᵣ. The grouping symbol
`:filled` represents the eigenmodes that corresponds to γ ≈ 0; `:empty` represents the eigenmodes that corresponds to γ ≈ 1.
"""
function frozenselectionbythreshold(threshold::Float64)::Function
    function frozenfockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
        spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
        filledmodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second < threshold, spectrum |> geteigenvalues)))
        emptymodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second > 1.0 - threshold, spectrum |> geteigenvalues)))
        return Dict(:filled => columns(spectrum |> geteigenvectors, FockSpace(filledmodes)), :empty => columns(spectrum |> geteigenvectors, FockSpace(emptymodes)))
    end
    return frozenfockmaps
end
export frozenselectionbythreshold

function frozenselectionbycount(count::Integer)
    function frozenfockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
        𝑈ᵣ::FockMap = eigvecsh(𝐶ᵣ)
        modes::Vector{Mode} = [orderedmodes(𝑈ᵣ.inspace)...]
        return Dict(:filled => columns(𝑈ᵣ, FockSpace(Subset(modes[1:count]))), :empty => columns(𝑈ᵣ, FockSpace(Subset(modes[(end - count + 1):end]))))
    end
    return frozenfockmaps
end
export frozenselectionbycount

"""
    regioncorrelations(correlations::FockMap, regionfock::FockSpace)::FockMap

Restrict the crystal correlations to a real space regional local correlations.

### Input
- `correlations::FockMap`: The crystal correlations.
- `regionfock::RegionFock`: The real space regional `RegionFock`.

### Output
The real space local correlation with `inspace` and `outspace` as the `regionfock`.
"""
function regioncorrelations(correlations::FockMap, regionfock::RegionFock)::FockMap
    crystalfock::CrystalFock = correlations|>getinspace
    transform = fourier(crystalfock, regionfock) / (crystalfock|>subspacecount|>sqrt)
    return transform' * correlations * transform
end
export regioncorrelations

localfrozenisometries(
    correlations::FockMap, regionfock::FockSpace;
    selectionstrategy::Function = frozenselectionbythreshold(1e-3))::Dict{Symbol, FockMap} = (
    regioncorrelations(correlations, regionfock) |> selectionstrategy)
export localfrozenisometries

blocking(parameters::Pair{Symbol}...)::Dict{Symbol} = blocking(Dict(parameters...))
export blocking

function blocking(parameters::Dict{Symbol})::Dict{Symbol}
    @assert(haskey(parameters, :action))
    @assert(haskey(parameters, :correlations))
    @assert(haskey(parameters, :crystal))
    result::Dict{Symbol, Any} = Dict()

    scale::Scale = parameters[:action]
    correlation::FockMap = parameters[:correlations]
    crystal::Crystal = parameters[:crystal]

    scaling::FockMap = scale * correlation.inspace
    result[:action] = scaling
    result[:transformer] = scaling
    result[:correlations] = scaling * correlation * scaling'
    result[:crystal] = scale * crystal

    return result
end

"""
    crystalisometries(; localisometry::FockMap, crystalfock::CrystalFock, addinspacemomentuminfo::Bool = false)::Dict{Momentum, FockMap}

Given a `localisometry` that is defined within a `RegionFock`, perform fourier transform to obtain the k-space representation of the isometry
and store as a dictionary of bloch isometries keyed by the momentum of the `crystalfock`.

### Input
- `localisometry::FockMap`: The local isometry defined within a `RegionFock`.
- `crystalfock::CrystalFock`: The crystal `FockSpace`.
- `addinspacemomentuminfo::Bool`: Whether to add the momentum information into the `inspace` as attribute `:k` of the returned `FockMap`.
"""
function crystalisometries(; localisometry::FockMap, crystalfock::CrystalFock,
    addinspacemomentuminfo::Bool = false)

    crystal::Crystal = getcrystal(crystalfock)
    transform::FockMap = fourier(crystalfock, localisometry|>getoutspace|>RegionFock) / (crystal |> vol |> sqrt)

    function preprocesslocalisometry(k::Momentum)::FockMap
        if !addinspacemomentuminfo
            return localisometry
        end
        inspace::FockSpace = localisometry.inspace |> orderedmodes |> setattr(:k => k) |> removeattr(:r) |> FockSpace
        return FockMap(localisometry, inspace=inspace, permute=false)
    end

    isometries = paralleltasks(
        name="crystalisometries",
        tasks=(()->(k=>transform[getsubspace(crystalfock, k), :]*preprocesslocalisometry(k)) for k in crystal|>brillouinzone),
        count=crystal|>vol)|>parallel

    return isometries
end
export crystalisometries

function crystalprojector(; localisometry::FockMap, crystalfock::CrystalFock)::FockMap
    momentumisometries = crystalisometries(localisometry=localisometry, crystalfock=crystalfock)
    crystal::Crystal = getcrystal(crystalfock)

    projectors::Dict = paralleltasks(
        name="crystalprojector",
        tasks=(()->((k, k)=>isometry * isometry') for (k, isometry) in momentumisometries),
        count=crystal|>vol)|>parallel|>Dict

    return crystalfockmap(crystal, crystal, projectors)
end
export crystalprojector

function crystalprojector(spectrum::CrystalSpectrum)::FockMap
    blocks = paralleltasks(
        name="crystalprojector",
        tasks=(()->((k, k)=>u*u') for (k, u) in spectrum|>geteigenvectors),
        count=spectrum|>getcrystal|>vol)|>parallel|>Dict
    return crystalfockmap(spectrum|>getcrystal, spectrum|>getcrystal, blocks)
end

function globaldistillerhamiltonian(;
    correlations::FockMap, restrictspace::FockSpace,
    localisometryselectionstrategy::Function, manualeigenenergies::Dict{Symbol, <:Number} = Dict(:filled => -1, :empty => 1))::FockMap
    @info("globaldistillerhamiltonian: Computing localisometries...")
    localisometries::Dict{Symbol} = localfrozenisometries(correlations, restrictspace, selectionstrategy=localisometryselectionstrategy)
    crystalfock::CrystalFock = correlations|>getinspace
    @info("globaldistillerhamiltonian: Computing crystalprojectors...")
    crystalprojectors = Dict(
        name=>crystalprojector(localisometry=isometry, crystalfock=crystalfock) for (name, isometry) in localisometries)
        @info("globaldistillerhamiltonian: finalizing...")
    return reduce(+, manualeigenenergies[name] * crystalprojector for (name, crystalprojector) in crystalprojectors)
end
export globaldistillerhamiltonian

function generategroupingfunction(grouppredicates)
    labels::Vector{Symbol} = [label for (label, _) in grouppredicates]
    predicates::Vector{Function} = [predicate for (_, predicate) in grouppredicates]
    function f(value::Number)::Symbol
        # This modification is added to support grouping the eigenvalues by rules that
        # are not fully include the entire spectrum, all remaining elements will be grouped
        # into the `:others` group.
        index = findfirst(p -> value |> p, predicates)
        return index|>isnothing ? :others : labels[index]
    end
    return f
end

function distillation(spectrum::CrystalSpectrum, bandpredicates...)::Dict{Symbol, CrystalSpectrum}
    groupingfunction::Function = bandpredicates |> generategroupingfunction
    labeled::Base.Generator = ((mode |> getattr(:k), v |> groupingfunction) => mode for (mode, v) in spectrum |> geteigenvalues)
    momentumgroups::Dict{Tuple, Vector{Mode}} = Dict()
    
    watchprogress(desc="distillation: labelling")
    for (key, mode) in labeled
        !haskey(momentumgroups, key) && (momentumgroups[key] = [])
        push!(momentumgroups[key], mode)
        updateprogress()
    end
    unwatchprogress()

    bands::Dict{Symbol, Dict{Momentum, FockMap}} = Dict()
    watchprogress(desc="distillation: band grouping")
    for ((k, band), modes) in momentumgroups
        !haskey(bands, band) && (bands[band] = Dict())
        bands[band][k] = columns((spectrum |> geteigenvectors)[k], modes |> FockSpace)
        updateprogress()
    end
    unwatchprogress()

    function repacktospectrum(isometries::Dict{Momentum, FockMap})::CrystalSpectrum
        eigenmodes::Dict{Momentum, Subset{Mode}} = Dict(k => fockmap |> getinspace |> orderedmodes for (k, fockmap) in isometries)
        eigenvalues::Dict{Mode, Number} = Dict(mode => (spectrum |> geteigenvalues)[mode] for (_, modes) in eigenmodes for mode in modes)
        return CrystalSpectrum(spectrum.crystal, eigenmodes, eigenvalues, isometries)
    end

    watchprogress(desc="distillation: finalizing")
    result = Dict()
    for (band, isometries) in bands
        result[band] = isometries|>repacktospectrum
        updateprogress()
    end
    unwatchprogress()
    return result
end
export distillation

"""
    findlocalspstates(; kwargs...)::Dict{Integer, FockMap}

Search for local degenerate groups of single particle states within the restricted local correlations.

### Input
- `statecorrelations::FockMap`: The source crystal correlations.
- `regionfock::FockSpace`: The restriction fockspace.
- `symmetry::AffineTransform`: The symmetry of the region defined in `regionfock`.
- `spectrumextractpredicate::Function`: The predicate to extract the spectrum values.
- `linearindependencethreshold::Real`: The threshold for the linear independence of the local single particle states.
- `degeneracythreshold::Real`: The threshold for the correlation value degeneracy of the local single particle states.

### Output
A dictionary of the local single particle states keyed by their group sizes.
"""
function findlocalspstates(;
    statecorrelations::FockMap, regionfock::FockSpace,
    symmetry::AffineTransform = identitytransform(statecorrelations|>getoutspace|>getcrystal|>dimension),
    spectrumextractpredicate::Function = v -> v < 1e-2,
    linearindependencethreshold::Real = 5e-2,
    degeneracythreshold::Real = 1e-7,
    statecrystalfock::CrystalFock = statecorrelations|>getoutspace)::Dict{Integer, FockMap}

    function lineardependencefilter(spstate::FockMap)::Bool
        crystalspstates = crystalisometries(localisometry=spstate, crystalfock=statecrystalfock)
        crystalspstate::FockMap = directsum(v for (_, v) in crystalspstates)
        pseudoidentity::FockMap = (crystalspstate' * crystalspstate)
        mineigenvalue = minimum(v for (_, v) in pseudoidentity |> eigvalsh)
        return mineigenvalue > linearindependencethreshold
    end

    localcorrelations::FockMap = regioncorrelations(statecorrelations, regionfock)
    localspectrum::EigenSpectrum = eigspech(localcorrelations, groupingthreshold=degeneracythreshold)
    groupeigenvalues::Base.Generator = (
        subset => (localspectrum |> geteigenvalues)[subset |> first]
        for subset in localspectrum |> geteigenvectors |> getinspace |> sparsegrouping(:eigenindex))
    selectedgroups = Iterators.filter(p -> p.second |> spectrumextractpredicate, groupeigenvalues)

    selectedisometries = ((localspectrum |> geteigenvectors)[:, group.first |> FockSpace] for group in selectedgroups)
    orthogonalspstates = Iterators.filter(lineardependencefilter, selectedisometries)
    symmetricspstates = (state * *(state, symmetry) for state in orthogonalspstates)
    spstates = (state * spatialmap(state) for state in symmetricspstates)

    return Dict(state |> getinspace |> dimension => state for state in spstates)
end
export findlocalspstates

function wannierprojection(;
    crystalisometries::Dict{Momentum, <:FockMap}, crystal::Crystal, crystalseeds::Dict{Momentum, <:FockMap}, svdorthothreshold::Number = 1e-1)

    wannierunitcell::Subset{Offset} = Subset(mode |> getattr(:b) for mode in (crystalseeds |> first |> last).inspace |> orderedmodes)
    wanniercrystal::Crystal = Crystal(wannierunitcell, crystal|>getbc)
    precarioussvdvalues::Vector = []
    function approximateisometry(k::Momentum, isometry::FockMap, overlap::FockMap)::FockMap
        U, Σ, Vt = overlap |> svd
        minsvdvalue::Number = minimum(v for (_, v) in Σ)
        if minsvdvalue < svdorthothreshold
            push!(precarioussvdvalues, minsvdvalue)
        end
        unitary::FockMap = U * Vt
        approximated::FockMap = isometry * unitary

        inspace=FockSpace(approximated|>getinspace|>mapmodes(m->m|>setattr(:k=>k)|>removeattr(:r)))
        return FockMap(approximated, inspace=inspace, permute=false)
    end
    if (precarioussvdvalues |> length) > 0
        @warn "Precarious wannier projection with minimum svdvalue of $(precarioussvdvalues |> minimum)"
    end

    blocks = paralleltasks(
        name="wannierprojection",
        tasks=(()->(
            (k, k)=>approximateisometry(k, isometry, isometry' * crystalseeds[k])) 
            for (k, isometry) in crystalisometries),
        count=crystalisometries|>length)|>parallel|>Dict
    
    return crystalfockmap(crystal, wanniercrystal, blocks)
end
export wannierprojection

function wannierprojection(targetbands::CrystalSpectrum, seedstates::RegionState)
    localseeds = seedstates|>FockMap

    outspace = targetbands|>crystalprojector|>getoutspace
    transform = fourier(outspace, localseeds|>getoutspace)
    crystalseeds = transform * localseeds
    crystalseeds = Dict(k=>crystalseeds[k, :] for k in targetbands|>getcrystal|>brillouinzone)
    pseudoidens = (u' * u for (_, u) in crystalseeds)
    lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
    lineardepmetric < 0.5 && @warn "Wannierization local seeds linear-dependence metric: $lineardepmetric"

    return wannierprojection(
        crystalisometries=targetbands|>geteigenvectors,
        crystal=targetbands|>getcrystal, crystalseeds=crystalseeds)
end

function getlocalstates(crystalisometry, region::Region)::RegionState
    regionfock = getregionfock(crystalisometry|>getoutspace, region)
    leftrestrict = fourier(crystalisometry|>getoutspace, regionfock)
    rightfock = crystalisometry|>getinspace|>unitcellfock|>RegionFock
    rightrestrict = fourier(crystalisometry|>getinspace, rightfock)
    return leftrestrict'*crystalisometry*rightrestrict|>RegionState|>normalize
end
export getlocalstates

function getlocalstates(correlations; region::Region, count::Real)
    localfock = getregionfock(correlations|>getoutspace, region)
    restrict = fourier(correlations|>getoutspace, localfock)
    crystalvol = correlations|>getoutspace|>getcrystal|>vol
    localcorrelations = restrict'*correlations*restrict/crystalvol
    localspectrum = localcorrelations|>eigspech
    states = []
    for i in 1:count
        localvector = geteigenvectors(localspectrum)[:, i]
        m = localvector|>getinspace|>first
        v = geteigenvalues(localspectrum)[m]
        push!(states, (localvector|>RegionState, v))
    end

    @info "Showing local states with associated eigenvalues..."
    header = ["Index", "Eigenvalue"]
    rowlength = length(header)
    rows = (reshape([Int(n), val], (1, rowlength)) for (n, (_, val)) in states|>enumerate)
    pretty_table(vcat(rows...), header=header)

    return sum(state for (state, _) in states)
end

"""
    svdrebase(inputbands::CrystalSpectrum, targetbands::CrystalSpectrum)

Map the `inputbands` to the basis spanning the `targetbands` using the SVD projection.

### Output
The purified `CrystalSpectrum` object with the same span of the input bands.
"""
function svdrebase(inputbands::CrystalSpectrum, targetbands::CrystalSpectrum)
    function rebase(k, ψ, Φ)
        ϕ = Φ[k]
        p = ϕ'ψ
        u, s, vd = p|>svd
        @assert all([v for (_, v) in s].>0.9)
        m = u*vd
        return k=>ϕ*m
    end

    purified = paralleltasks(
        name="svdrebase",
        tasks=(
            ()->rebase(k, ψ, targetbands|>geteigenvectors) 
            for (k, ψ) in inputbands|>geteigenvectors),
        count=inputbands|>geteigenvectors|>length)|>parallel|>Dict

    return CrystalSpectrum{inputbands|>getcrystal|>dimension}(
        inputbands|>getcrystal, inputbands|>geteigenmodes, inputbands|>geteigenvalues, 
        purified, inputbands.bandcount)
end
export svdrebase
