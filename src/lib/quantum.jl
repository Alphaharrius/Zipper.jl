"""
    colmap(inmode::Mode, rowdata::Vector{Pair{Mode, ComplexF64}})::FockMap

Create a column `FockMap` with the specified complex value entries associated with the modes which forms the outspace, with a dimension `1` inspace.

### Input
- `inmode`  The single `Mode` that will form the inspace of the fock map.
- `rowdata` The `Mode` to `ComplexF64` pairs which the modes forms the outspace of the fock map and the complex values as the entries.

### Output
A column `FockMap` with the outspace ordered by the order of `rowdata`.
"""
# TODO: Remove trivial methods.
function colmap(inmode::Mode, rowdata::Vector{Pair{Mode, ComplexF64}})::FockMap
    outfock::FockSpace = FockSpace(Subset(p.first for p in rowdata))
    spmat::SparseMatrixCSC{ComplexF64, Int64} = spzeros(dimension(outfock), 1)
    foreach(n -> spmat[n, 1] += rowdata[n].second, 1:dimension(outfock))
    return FockMap(outfock, FockSpace(Subset(inmode)), spmat)
end
export colmap

"""
    rowsubmaps(fockmap::FockMap)::Base.Generator

Partition the `fockmap` into smaller `fockmaps` by the subspaces of its `outspace`.
"""
# TODO: Trivial methods
rowsubmaps(fockmap::FockMap)::Base.Generator = (rows(fockmap, subspace) for subspace in fockmap|>getoutspace |> subspaces)
export rowsubmaps

"""
    colsubmaps(fockmap::FockMap)::Base.Generator

Partition the `fockmap` into smaller `fockmaps` by the subspaces of its `inspace`.
"""
# TODO: Trivial methods
colsubmaps(fockmap::FockMap)::Base.Generator = (columns(fockmap, subspace) for subspace in fockmap|>getinspace |> subspaces)
export colsubmaps

function mapunitcellfock(to::FockSpace, from::FockSpace)
    tohomefock::FockSpace = to|>unitcellfock
    fromhomefock::FockSpace = from|>unitcellfock
    
    if hassamespan(tohomefock, fromhomefock) || issubspace(tohomefock, fromhomefock)
        return Dict(mode=>mode for mode in fromhomefock)
    end
    if issubspace(fromhomefock, tohomefock)
        return Dict(mode=>mode for mode in tohomefock)
    end

    error("This method requires the two fockspaces to have intersections!")
end
export mapunitcellfock

function rulemapunitcellfock(to::FockSpace, from::FockSpace, rule::Function = (a, b) -> (a == b))
    tohomefock::FockSpace = to|>unitcellfock
    fromhomefock::FockSpace = from|>unitcellfock

    pairs = product(tohomefock, fromhomefock)
    return Dict(tomode=>frommode for (tomode, frommode) in pairs if rule(tomode, frommode))
end
export rulemapunitcellfock

function fourier(crystalfock::CrystalFock, regionfock::RegionFock, unitcellfockmapping::Dict{Mode, Mode} = mapunitcellfock(crystalfock, regionfock))
    momentummatrix::SparseMatrixCSC = crystalfock|>getcrystal|>computemomentummatrix

    momentumhomefock = crystalfock|>unitcellfock
    values::Array = zeros(Complex, momentumhomefock|>length, size(momentummatrix, 2), regionfock|>dimension)
    
    function fillvalues(n, homemode, m, inmode)
        if !haskey(unitcellfockmapping, homemode) || unitcellfockmapping[homemode] != inmode|>removeattr(:r) return end
        offsetvector::Vector = inmode|>getattr(:r)|>euclidean|>vec
        values[n, :, m] = exp.(-1im * momentummatrix' * offsetvector)
    end

    # Since each (n, m) only corresponds to one entry, thus this is thread-safe.
    paralleltasks(
        name="fourier $(crystalfock|>dimension)×$(regionfock|>dimension)",
        # TODO: Remove test code for isolating the issue.
        tasks=(()->fillvalues(n, homemode, m, inmode) for ((n, homemode), (m, inmode)) in Iterators.product(momentumhomefock|>enumerate, regionfock|>enumerate)),
        count=dimension(momentumhomefock)*dimension(regionfock))|>parallel

    data::SparseMatrixCSC = reshape(values, (length(momentumhomefock) * size(momentummatrix, 2), regionfock|>dimension))|>SparseMatrixCSC
    return FockMap(crystalfock, regionfock, data)
end
export fourier

"""
    crystalsubmaps(fockmap::FockMap)

Given a `FockMap` with `inspace` and `outspace` of type `CrystalFock` of same span, partition the `FockMap`
into smaller `FockMap` objects by the subspaces indexed by the `Momentum` attribute.

### Output
A generator yielding `Pair{Momentum, FockMap}` objects, with the momentums corresponds to thw brillouin zone.
"""
function crystalsubmaps(fockmap::FockMap)::Base.Generator
    (fockmap|>getinspace isa CrystalFock && fockmap|>getoutspace isa CrystalFock && hassamespan(fockmap|>getinspace, fockmap|>getoutspace) ||
        error("The in/out spaces of the fock map must be the same crystal fock-spaces!"))
    hassamespan(fockmap|>getinspace, fockmap|>getoutspace) || error("Required a FockMap with the same in/out CrystalFock!")
    return (k => restrict(fockmap, fockspace, fockspace) for (k, fockspace) in fockmap|>getinspace|>crystalsubspaces|>Dict)
end
export crystalsubmaps



"""
    columnspec(fockmap::FockMap)::Vector{Pair{Mode, ComplexF64}}

Extract the column values as a `Mode` to `ComplexF64` pair from a `N×1` `FockMap`, this is used when you need to visualize the column spectrum.
"""
function columnspec(fockmap::FockMap)::Vector{Pair{Mode, ComplexF64}}
    @assert(dimension(fockmap|>getinspace) == 1)
    mat::SparseMatrixCSC{ComplexF64, Int64} = rep(fockmap)
    return [outmode => mat[(fockmap|>getoutspace)[outmode], 1] for outmode in orderedmodes(fockmap|>getoutspace)]
end

""" Implemented to support Fourier transform of `CrystalFockMap` objects. """
function Base.:*(fouriertransform::FockMap{RegionFock, CrystalFock}, fockmap::CrystalFockMap)
    crystal::Crystal = fouriertransform|>getinspace|>getcrystal
    realspace::RealSpace = crystal|>getspace
    kspace::MomentumSpace = convert(MomentumSpace, realspace)
    Γ::Momentum = kspace|>getorigin
    singularcrystal::Crystal = Crystal(realspace|>getorigin|>Subset, realspace|>dimension|>ones)

    blocks::Dict = paralleltasks(
        name="FockMap{RegionFock, CrystalFock} * CrystalFockMap",
        tasks=(()->((Γ, k)=>fouriertransform[:, getsubspace(fouriertransform|>getinspace, k)]) for k in crystal|>brillouinzone),
        count=crystal|>vol)|>parallel|>Dict
    crystaltransform::CrystalFockMap = CrystalFockMap(singularcrystal, crystal, blocks)

    transformed::CrystalFockMap = crystaltransform * fockmap
    outspace::RegionFock = fouriertransform|>getoutspace
    inspace::CrystalFock = fouriertransform|>getinspace

    function compute(blocks)
        data::SparseMatrixCSC = spzeros(Complex, outspace|>dimension, inspace|>dimension)
        for block in blocks
            blockinspace::FockSpace = block|>getinspace
            inorder::UnitRange = inspace[blockinspace|>first]:inspace[blockinspace|>last]
            data[:, inorder] += block|>rep
        end
        return data
    end

    # Without doing the followings this step is unbearablelly slow for large system size.
    batchsize::Integer = (length(transformed.blocks) / getmaxthreads())|>ceil
    summingbatches = Iterators.partition((block for (_, block) in transformed.blocks), batchsize)
    reps = paralleltasks(
        name="FockMap{RegionFock, CrystalFock} * CrystalFockMap",
        tasks=(()->compute(partition) for partition in summingbatches),
        count=getmaxthreads())|>parallel
    
    watchprogress(desc="FockMap{RegionFock, CrystalFock} * CrystalFockMap")
    spdata::SparseMatrixCSC = spzeros(Complex, outspace|>dimension, inspace|>dimension)
    for rep in reps
        spdata[:, :] += rep
        updateprogress()
    end
    unwatchprogress()

    return FockMap(outspace, inspace, spdata)
end
