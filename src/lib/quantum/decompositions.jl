# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ EigenSpectrum definition ◆
""" Packaging the result computed from eigenvalue decomposition. """
struct EigenSpectrum
    eigenvalues::Dict{Mode, Number}
    eigenvectors::FockMap
end
export EigenSpectrum
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ EigenSpectrum APIs ◆
function geteigenmodes(spectrum::EigenSpectrum)::Subset{Mode}
    return spectrum.eigenvectors|>getinspace |> orderedmodes
end
export geteigenmodes

geteigenvalues(spectrum::EigenSpectrum)::Dict{Mode, Number} = spectrum.eigenvalues
export geteigenvalues

geteigenvectors(spectrum::EigenSpectrum)::FockMap = spectrum.eigenvectors
export geteigenvectors

"""
    eigspech(hermitian::FockMap, attrs::Pair{Symbol}...; groupingthreshold::Real = 1e-7)::EigenSpectrum

Perform Hermitian eigenvalue decomposition to find the eigenvalues and eigenvectors simultaneously, the 
corresponding eigenmodes will have the attributes of `:eigenindex` which corresponds to the degenerate 
group associated to a eigenvalue, and in ascending order of the eigenvalues; `:flavor` which indicates 
the individual degrees of freedom within the degenerate group; along with the attributes supplied by `attrs`.

### Input
- `hermitian`           The source of the decomposition, must be a Hermitian.
- `attrs`               Attributes to be inserted to the generated eigenmodes.
- `groupingthreshold`   The threshold for grouping degenerated eigenvalues.

### Output
An `EigenSpectrum` object containing all the computed information.
"""
function eigspech(hermitian::FockMap, attrs::Pair{Symbol}...; groupingthreshold::Real = 1e-7)::EigenSpectrum
    vals, U = hermitian |> rep |> Matrix |> Hermitian |> eigen
    eigenvalues::Base.Iterators.Flatten = digesteigenvalues(Rational, Real, vals, groupingthreshold, attrs...)
    eigenvectors::FockMap = FockMap(hermitian |> getoutspace, FockSpace(m for (m, _) in eigenvalues), U)
    return EigenSpectrum(eigenvalues |> Dict, eigenvectors)
end
export eigspech

"""
    eigspec(fockmap::FockMap, attrs::Pair{Symbol}...; groupingthreshold::Real = 1e-7)::EigenSpectrum

Perform eigenvalue decomposition to find the eigenvalues and eigenvectors simultaneously, the corresponding 
eigenmodes will have the attributes of `:eigenindex` which corresponds to the degenerate group associated to 
a eigenvalue; `:flavor` which indicates the individual degrees of freedom within the degenerate group; along 
with the attributes supplied by `attrs`.

### Input
- `fockmap`             The source of the decomposition.
- `attrs`               Attributes to be inserted to the generated eigenmodes.
- `groupingthreshold`   The threshold for grouping degenerated eigenvalues, since the eigenvalues are complex 
                        numbers, this threshold will be applied to the real and imaginary parts separately.
"""
function eigspec(fockmap::FockMap, attrs::Pair{Symbol}...; groupingthreshold::Real = 1e-7)::EigenSpectrum
    vals, U = fockmap |> rep |> Matrix |> eigen
    eigenvalues::Base.Iterators.Flatten = digesteigenvalues(Tuple, Complex, vals, groupingthreshold, attrs...)
    eigenvectors::FockMap = FockMap(fockmap |> getoutspace, FockSpace(m for (m, _) in eigenvalues), U)
    return EigenSpectrum(eigenvalues |> Dict, eigenvectors)
end
export eigspec

"""
    groupbyeigenvalues(spectrum; groupingthreshold::Number = 1e-7)

Given a spectrum, attempt to group the eigenmodes based on their corresponding eigenvalues with a 
eigenvalue grouping threshold.

### Output
A generator yielding `Pair{Number, Subset{Mode}}` objects, with the eigenvalues as keys and the 
corresponding eigenmodes as values.
"""
function groupbyeigenvalues(spectrum; groupingthreshold::Number = 1e-7)::Base.Generator
    denominator::Integer = (1 / groupingthreshold)|>round|>Integer
    actualvalues::Dict{Rational, Number} = Dict(
        hashablereal(v, denominator)=>v for (_, v) in spectrum|>geteigenvalues)
    items::Base.Generator = (hashablereal(v, denominator)=>m for (m, v) in spectrum|>geteigenvalues)
    groups::Dict{Rational, Vector{Mode}} = foldl(items; init=Dict{Rational, Vector{Mode}}()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k=>[v]))
    end
    sortedrationals::Vector{Rational} = sort([(groups|>keys)...])
    return (actualvalues[r]=>groups[r]|>Subset for r in sortedrationals)
end
export groupbyeigenvalues
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Internal methods ◆
""" Internal method that generates the `eigenmode => eigenvalue` pairs. """
function digesteigenvalues(H::Type, V::Type, vals, groupingthreshold::Real, attrs::Pair{Symbol}...)
    denominator = (1 / groupingthreshold) |> round |> Integer
    valtable::Dict{H, V} = Dict(hashablenumber(v |> V, denominator) => v for v in vals)
    items::Base.Generator = (hashablenumber(v |> V, denominator) => n for (n, v) in vals |> enumerate)
    groups::Dict{H, Vector} = foldl(items; init=Dict{H, Vector}()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end
    sortedgroups::Vector = sort([groups...], by=(g -> g.first))

    return (
        Mode(:eigenindex => n, :flavor => f, attrs...) => valtable[group.first]
        for (n, group) in sortedgroups |> enumerate
        for f in group.second |> eachindex)
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ EigenSpectrum display ◆
Base.:show(io::IO, spectrum::EigenSpectrum) = print(
    io, string("$(spectrum |> typeof)(entries=$(spectrum.eigenvalues |> length))"))
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalSpectrum definition ◆
"""
Decomposition of `FockMap` with `inspace` and `outspace` of type `CrystalFock` of same span,
and store the underlying information into eigen value decomposed form indexed by the `Momentum` of the
brillouin zone of the crystal. The number of eigenmodes per `Momentum` is allowed to be different from
the amount of modes in the unit cell `FockSpace` of the `CrystalFock`, which corresponds to a projector
if being converted back to a `FockMap`, and the corresponding eigenvalues are `1`.

Packing information into a `CrystalSpectrum` allows visualization of eigen spectrum in a band diagram.
"""
struct CrystalSpectrum{Dim}
    crystal::Crystal
    eigenmodes::Dict{Momentum, Subset{Mode}}
    eigenvalues::Dict{Mode, Number}
    eigenvectors::Dict{Momentum, FockMap}
    bandcount::Tuple
end
export CrystalSpectrum

function CrystalSpectrum(crystal::Crystal, crystaleigenmodes, args...)
    maxbands = 0
    minbands = typemax(Int)
    watchprogress(desc="CrystalSpectrum")
    for k in crystal|>brillouinzone
        if !haskey(crystaleigenmodes, k)
            minbands = 0
            continue
        end
        modes = crystaleigenmodes[k]
        count = modes|>length
        maxbands = max(maxbands, count)
        minbands = min(minbands, count)
        updateprogress()
    end
    unwatchprogress()
    return CrystalSpectrum{crystal|>dimension}(crystal, crystaleigenmodes, args..., (minbands, maxbands))
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalSpectrum APIs ◆
"""
    crystalspectrum(momentumfockmaps; crystal::Crystal)::CrystalSpectrum

Given a collection of Hermitian `FockMap` objects each associated with a specific `Momentum` from the brillouin zone,
pack into a `CrystalSpectrum` object.
"""
function crystalspectrum(momentumfockmaps; crystal::Crystal)::CrystalSpectrum
    function compute(k, fockmap)
        eigenspectrum::EigenSpectrum = eigspech(fockmap, :k => k)
        eigenmodes = Subset(m for (m, _) in eigenspectrum |> geteigenvalues)
        eigenvectors = eigenspectrum|>geteigenvectors
        eigenvalues = (mode=>v for (mode, v) in eigenspectrum|>geteigenvalues)
        return k=>(eigenmodes, eigenvectors, eigenvalues)
    end
    tasks = (()->compute(k, fockmap) for (k, fockmap) in momentumfockmaps)
    eigendatas = paralleltasks(name="crystalspectrum", tasks=tasks, count=crystal|>vol)|>parallel

    crystaleigenmodes::Dict{Momentum, Subset{Mode}} = Dict(k=>kmodes for (k, (kmodes, _, _)) in eigendatas)
    crystaleigenvalues::Dict{Mode, Number} = Dict(m=>v for (k, (_, _, pairs)) in eigendatas for (m, v) in pairs)
    crystaleigenvectors::Dict{Momentum, FockMap} = Dict(k=>kvecs for (k, (_, kvecs, _)) in eigendatas)
    return CrystalSpectrum(crystal, crystaleigenmodes, crystaleigenvalues, crystaleigenvectors)
end
export crystalspectrum

"""
    getbandcount(spectrum::CrystalSpectrum)

Get the band count of the spectrum, if the band count is not uniform 
over the brillouin zone, an error will be thrown.
"""
function getbandcount(spectrum::CrystalSpectrum)
    if first(spectrum.bandcount) != last(spectrum.bandcount)
        error("Band count is not uniform over the brillouin zone!")
    end
    return first(spectrum.bandcount)
end
export getbandcount

"""
    crystalspectrum(fockmap::FockMap)::CrystalSpectrum

Given a Hermitian `FockMap` with `inspace` and `outspace` of type 
`CrystalFock` of same span, pack into a `CrystalSpectrum` object.
"""
crystalspectrum(fockmap::FockMap)::CrystalSpectrum = crystalspectrum(
    fockmap|>crystalsubmaps, crystal=fockmap|>getinspace|>getcrystal)

crystalspectrum(fockmap::CrystalFockMap)::CrystalSpectrum = crystalspectrum(
    fockmap|>crystalsubmaps, crystal=fockmap|>getoutspace|>getcrystal)

"""
    linespectrum(spectrum::CrystalSpectrum)::CrystalSpectrum{1}

Since some spectrum might be defined in a crystal that is implicitly on 1-D space, such as a line in 
the 2D plane, inorder to visualize the spectrum properly instead of using the plot in the original 
dimension, we can embed the spectrum into a 1-D space for plotting.
"""
function linespectrum(spectrum::CrystalSpectrum)
    crystal::Crystal = spectrum|>getcrystal
    nontrivialbasis = Iterators.filter(
        v -> (v[1]) == 1, zip(crystal|>getboundsize, crystal|>getspace|>getbasisvectors))|>collect
    @assert(nontrivialbasis|>collect|>length == 1, "More than 1 non-trivial basis in this crystal.")
    embeddedbasis::Real = nontrivialbasis[1][2]|>norm
    embeddedspace::RealSpace = RealSpace([embeddedbasis][:, :])
    unitcelllength = crystal|>getunitcell|>length
    dummyunitcell::Region = Subset([r / unitcelllength] ∈ embeddedspace for r in (0:unitcelllength - 1))
    embeddedcrystal::Crystal = Crystal(dummyunitcell, crystal|>getbc)
    return CrystalSpectrum(
        embeddedcrystal, spectrum|>geteigenmodes, spectrum|>geteigenvalues, spectrum|>geteigenvectors)
end
export linespectrum

"""
    groupbands(spectrum::CrystalSpectrum)

Group the bands of the `CrystalSpectrum` by a specified band grouping strategy.

### Keyword Arguments
- `bandgrouping` The grouping strategy specified by `[count1, count2, ...]`.
- `frombottom`   Whether to group the bands from the bottom of the spectrum, if set to `false` this will group 
                 from the top of the spectrum.
"""
function groupbands(spectrum::CrystalSpectrum; bandgrouping=[spectrum.bandcount|>last], frombottom::Bool = true)
    kvalues = (
        k=>(m=>geteigenvalues(spectrum)[m] for m in (spectrum|>geteigenmodes)[k]) 
        for k in spectrum|>getcrystal|>brillouinzone)
    
    function grouping(k, iter)
        sorted = sort(iter|>collect, by=last)
        frombottom || reverse!(sorted)
        groupedvalues = []
        n = 1
        for g in bandgrouping
            group = sorted[n:n+g-1]
            push!(groupedvalues, group)
            n += g
        end
        groupedmodes = [Subset(m for (m, _) in group) for group in groupedvalues]
        groupedvectors = [geteigenvectors(spectrum)[k][:, modes|>FockSpace] for modes in groupedmodes]
        return k=>(groupedvalues, groupedmodes, groupedvectors)
    end

    kgrouped = paralleltasks(name="groupbands",
        tasks=(()->grouping(k, iter) for (k, iter) in kvalues),
        count=spectrum|>getcrystal|>vol)|>parallel|>Dict
    
    function getgroupspectrum(n)
        eigenvalues = Dict()
        eigenmodes = Dict()
        eigenvectors = Dict()
        for k in spectrum|>getcrystal|>brillouinzone
            values, modes, vectors = kgrouped[k]
            eigenmodes[k] = modes[n]
            eigenvectors[k] = vectors[n]
            for (m, v) in values[n]
                eigenvalues[m] = v
            end
        end
        return CrystalSpectrum(spectrum|>getcrystal, eigenmodes, eigenvalues, eigenvectors)
    end

    return paralleltasks(name="groupbands",
        tasks=(()->getgroupspectrum(n) for n in 1:length(bandgrouping)),
        count=length(bandgrouping))|>parallel|>collect
end
export groupbands

"""
    groupbands(spectrum::CrystalSpectrum, bandpredicates...)

Group the bands of the `CrystalSpectrum` by a specified band grouping strategy.

### Input
- `spectrum`        The source of the grouping.
- `bandpredicates`  Predicates to group the bands, each predicate should be a function that 
                    takes a single eigenvalue and returns a boolean value, and the structure of 
                    this argument is a `Pair` `bandname=>predicate`.

### Output
A dictionary of `bandname=>CrystalSpectrum` objects, each containing the eigenmodes and eigenvalues, 
the bands that is excluded from all the required bands will be grouped as `:others`.
"""
function groupbands(spectrum::CrystalSpectrum, bandpredicates...)::Dict{Symbol, CrystalSpectrum}
    labels::Vector{Symbol} = [label for (label, _) in bandpredicates]
    predicates::Vector{Function} = [predicate for (_, predicate) in bandpredicates]

    function classifier(value::Number)::Symbol
        # This modification is added to support grouping the eigenvalues by rules that
        # are not fully include the entire spectrum, all remaining elements will be grouped
        # into the `:others` group.
        index = findfirst(p -> value |> p, predicates)
        return index|>isnothing ? :others : labels[index]
    end

    labeled::Base.Generator = ((mode|>getattr(:k), v|>classifier) => mode for (mode, v) in spectrum|>geteigenvalues)
    momentumgroups::Dict{Tuple, Vector{Mode}} = Dict()
    
    watchprogress(desc="groupbands: labelling")
    for (key, mode) in labeled
        !haskey(momentumgroups, key) && (momentumgroups[key] = [])
        push!(momentumgroups[key], mode)
        updateprogress()
    end
    unwatchprogress()

    bands::Dict{Symbol, Dict{Momentum, FockMap}} = Dict()
    watchprogress(desc="groupbands: band grouping")
    for ((k, band), modes) in momentumgroups
        !haskey(bands, band) && (bands[band] = Dict())
        bands[band][k] = columns((spectrum |> geteigenvectors)[k], modes |> FockSpace)
        updateprogress()
    end
    unwatchprogress()

    function repacktospectrum(isometries::Dict{Momentum, FockMap})::CrystalSpectrum
        eigenmodes::Dict{Momentum, Subset{Mode}} = Dict(
            k => fockmap |> getinspace |> orderedmodes for (k, fockmap) in isometries)
        eigenvalues::Dict{Mode, Number} = Dict(
            mode => (spectrum |> geteigenvalues)[mode] for (_, modes) in eigenmodes for mode in modes)
        return CrystalSpectrum(spectrum.crystal, eigenmodes, eigenvalues, isometries)
    end

    watchprogress(desc="groupbands: finalizing")
    result = Dict()
    for (band, isometries) in bands
        result[band] = isometries|>repacktospectrum
        updateprogress()
    end
    unwatchprogress()
    return result
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalSpectrum implementation of spatial interface ◆
""" Get the associated crystal of a `CrystalSpectrum` object. """
Zipper.:getcrystal(spectrum::CrystalSpectrum)::Crystal = spectrum.crystal
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalSpectrum extension to EigenSpectrum APIs ◆
"""
    geteigenmodes(spectrum::CrystalSpectrum)::Dict{Momentum, Subset{Mode}}

Get the eigenmodes of a `CrystalSpectrum` object, indexed by the `Momentum` of 
the brillouin zone of the crystal.

### Output
A dictionary with keys of `Momentum` and values of `Subset{Mode}` which contains 
the momentum indexed unitcell modes within the crystal.
"""
geteigenmodes(spectrum::CrystalSpectrum)::Dict{Momentum, Subset{Mode}} = spectrum.eigenmodes

"""
    geteigenvalues(spectrum::CrystalSpectrum)::Dict{Mode, Number}

Get the eigenvalues of a `CrystalSpectrum` object, indexed by all the momentum 
indexed `Mode` objects of the `CrystalFock`.
"""
geteigenvalues(spectrum::CrystalSpectrum)::Dict{Mode, Number} = spectrum.eigenvalues

"""
    geteigenvectors(spectrum::CrystalSpectrum)::Dict{Momentum, FockMap}

Get the eigenvectors associated with each indexing `Momentum` within the brillouin zone, 
each with `inspace` corresponds to the returned `Subset{Mode}` of the same indexing `Momentum` 
from `geteigenmodes(spectrum)`.
"""
geteigenvectors(spectrum::CrystalSpectrum)::Dict{Momentum, FockMap} = spectrum.eigenvectors
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalSpectrum extension to other interfaces ◆
""" Shorthand to retrieve the unitcell fockspace from a `CrystalSpectrum`. """
function unitcellfock(spectrum::CrystalSpectrum)
    sourcefock::FockSpace = spectrum |> geteigenvectors |> first |> last |> getoutspace
    return sourcefock|>orderedmodes|>removeattr(:k)|>FockSpace
end

"""
    Fockmap(spectrum::CrystalSpectrum)::FockMap

Pack the `CrystalSpectrum` into a `FockMap` with `outspace` and `inspace` of type `CrystalFock` 
of same span, the unitcell fockspace of the packed `FockMap` should corresponds directly to the 
`outspace` of individual engenvectors in the `CrystalSpectrum`.
"""
function FockMap(crystalspectrum::CrystalSpectrum)::FockMap
    function momentumfockmap(k::Momentum)
        modes::Subset{Mode} = crystalspectrum.eigenmodes[k]
        eigenfock::FockSpace = modes |> FockSpace
        diagonal::FockMap = FockMap(
            eigenfock, eigenfock, Dict((m, m) => crystalspectrum.eigenvalues[m] |> ComplexF64 for m in modes))
        return crystalspectrum.eigenvectors[k] * diagonal * crystalspectrum.eigenvectors[k]'
    end
    fockmap::FockMap = directsum(k |> momentumfockmap for (k, _) in crystalspectrum.eigenmodes)
    crystalfock::FockSpace = FockSpace(fockmap|>getinspace, reflected=crystalspectrum.crystal)
    return FockMap(fockmap, inspace=crystalfock, outspace=crystalfock, permute=false)
end

function CrystalFockMap(crystalspectrum::CrystalSpectrum)::CrystalFockMap
    eigenvectors = crystalspectrum|>geteigenvectors

    function compute(k::Momentum)
        modes::Subset{Mode} = crystalspectrum.eigenmodes[k]
        eigenfock::FockSpace = modes |> FockSpace
        diagonal::FockMap = FockMap(
            eigenfock, eigenfock,
            Dict((m, m) => crystalspectrum.eigenvalues[m]|>ComplexF64 for m in modes))
        return (k, k)=>(eigenvectors[k] * diagonal * eigenvectors[k]')
    end

    blocks::Dict = paralleltasks(
        name="CrystalFockMap from CrystalSpectrum",
        tasks=(()->compute(k) for (k, _) in crystalspectrum.eigenmodes),
        count=crystalspectrum|>getcrystal|>vol)|>parallel|>Dict

    return CrystalFockMap(crystalspectrum|>getcrystal, crystalspectrum|>getcrystal, blocks)
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalSpectrum display ◆
function Base.:show(io::IO, spectrum::CrystalSpectrum)
    minbands = spectrum.bandcount[1]
    maxbands = spectrum.bandcount[2]
    bandcountcontext = minbands == maxbands ? "$minbands" : "($minbands->$maxbands)"
    print(io, string(
        "$(spectrum |> typeof)(bandcount=$bandcountcontext, entries=$(spectrum.eigenvalues|>length))"))
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Decomposition APIs ◆
"""
    svd(fockmap::FockMap)::Tuple{FockMap, Base.Generator, FockMap}

Perform singular value decomposition on a FockMap and return the left singular vectors,
singular values, and right singular vectors.
"""
function LinearAlgebra.:svd(fockmap::FockMap)::Tuple{FockMap, Base.Generator, FockMap}
    leftmodes::Subset{Mode} = Subset(Mode([:svdindex => n]) for n in 1:dimension(fockmap|>getinspace))
    rightmodes::Subset{Mode} = Subset(Mode([:svdindex => n]) for n in 1:dimension(fockmap|>getoutspace))
    svdmodes = length(leftmodes) < length(rightmodes) ? leftmodes : rightmodes
    U, Σ, Vt = fockmap |> rep |> Matrix |> svd
    svdvalues::Base.Generator = (
        (Mode([:svdindex => n]), Mode([:svdindex => n])) => Σ[n] for n in 1:length(svdmodes))
    return (
        FockMap(fockmap|>getoutspace, svdmodes|>FockSpace, U),
        svdvalues,
        FockMap(svdmodes|>FockSpace, fockmap|>getinspace, Vt'))
end

""" Shorthand for retrieving the eigenvectors from the `eigspech` function. """
eigvecsh(hermitian::FockMap, attrs::Pair{Symbol}...)::FockMap = eigspech(hermitian, attrs...)|>geteigenvectors
export eigvecsh

""" Shorthand for retrieving the eigenvalues from the `eigspech` function. """
eigvalsh(hermitian::FockMap, attrs::Pair{Symbol}...)::Dict{Mode, Real} = eigspech(hermitian, attrs...)|>geteigenvalues
export eigvalsh
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
