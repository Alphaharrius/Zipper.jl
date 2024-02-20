



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
    makezero(fockmap::FockMap, eps::Number = 1e-7)::FockMap

Round all numerical zeros within the `FockMap` to actual zero with a given tolerance `eps`.
"""
makezero(fockmap::FockMap, eps::Number = 1e-7)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, map(v -> abs(v |> real) < eps && abs(v |> imag) < eps ? 0im : v, fockmap |> rep))
""" Shorthand returning a function that performs `makezero` with a given tolerance `eps` on parameter `fockmap`. """
makezero(eps::Number = 1e-7)::Function = fockmap -> makezero(fockmap, eps)
export makezero

"""
    isapproxzero(fockmap::FockMap, eps::Number = 1e-7)::Bool

Check if all data within the `FockMap` is numerical zero within a given tolerance `eps`.
"""
isapproxzero(fockmap::FockMap, eps::Number = 1e-7)::Bool = fockmap |> makezero(eps) |> iszero
""" Shorthand returning a function that performs `isapproxzero` with a given tolerance `eps` on parameter `fockmap`. """
isapproxzero(eps::Number = 1e-7)::Function = fockmap::FockMap -> approxzero(fockmap, eps)
export isapproxzero

""" Get the maximum data from the `FockMap` by absolute value. """
Base.:maximum(fockmap::FockMap)::Complex = (fockmap |> rep)[map(v -> v |> abs, fockmap |> rep) |> findmax |> last]
""" Get the minimum data from the `FockMap` by absolute value. """
Base.:minimum(fockmap::FockMap)::Complex = (fockmap |> rep)[map(v -> v |> abs, fockmap |> rep) |> findmin |> last]

"""
    commutation(a::FockMap, b::FockMap)::FockMap

Get the commutator of two `FockMap` objects, which is defined as `a * b - b * a`. 
"""
commutation(a::FockMap, b::FockMap)::FockMap = (a * b - b * a)
export commutation

"""
    commutation(a::FockMap, b::FockMap; eps::Number = 1e-7)::Bool

Check if two `FockMap` objects commutes with each other, up to a numerical zero tolerance `eps`.
"""
commute(a::FockMap, b::FockMap; eps::Number = 1e-7)::Bool = commutation(a, b) |> makezero(eps) |> iszero
export commute

"""
    columnspec(fockmap::FockMap)::Vector{Pair{Mode, ComplexF64}}

Extract the column values as a `Mode` to `ComplexF64` pair from a `N×1` `FockMap`, this is used when you need to visualize the column spectrum.
"""
function columnspec(fockmap::FockMap)::Vector{Pair{Mode, ComplexF64}}
    @assert(dimension(fockmap|>getinspace) == 1)
    mat::SparseMatrixCSC{ComplexF64, Int64} = rep(fockmap)
    return [outmode => mat[(fockmap|>getoutspace)[outmode], 1] for outmode in orderedmodes(fockmap|>getoutspace)]
end

"""
    spatialremapper(regionfock::RegionFock; offsets::Region, unitcell::Region)

This function remaps the positions of modes in a `RegionFock` object based on provided offsets and a unit cell. This function is used when 
the unit-cell defined in the original `regionfock` does not match with what you desire, most of the case the `unitcell` is not defined in 
the positive parallelogram spanned by the lattice basis vectors.

### Input
- `regionfock`  The `RegionFock` object to be remapped.

### Output
It creates a new `RegionFock` object with the remapped modes and returns an identity map between the new and original `RegionFock` objects.
"""
function spatialremapper(regionfock::RegionFock; offsets::Region, unitcell::Region)
    positions::Dict{Offset, Tuple} = Dict((r + b)=>(r, b) for (r, b) in Iterators.product(offsets, unitcell))
    remappingdata::Base.Generator = ((positions[mode|>getpos], mode) for mode in regionfock)
    remappedfock::RegionFock = RegionFock(mode|>setattr(:r=>r, :b=>b) for ((r, b), mode) in remappingdata)
    return idmap(remappedfock, regionfock)
end
export spatialremapper

"""
    mapmodes(mapper::Function, modes)::Subset{Mode}

A standard way to apply a mapping function onto an iterator of `Mode` that might ends up with mode duplication.
All duplications will be mapped with incremental `:flavor` indices starting from `1`. It is adviced that all code 
that looks like `Subset(mode|>mapper for mode in modes)` to be replaced with this function.

### Input
- `mapper`  The mapping function to be applied to each `Mode` object, it should only be taking one `Mode` as argument.
- `modes`   The iterator of `Mode` objects to be mapped.
"""
function mapmodes(mapper::Function, modes)::Subset{Mode}
    # :flavor is removed to ensure there is not overlooked degeneracy lifted.
    mapped::Base.Generator = (mode|>mapper|>removeattr(:flavor) for mode in modes)
    degenerates::Dict = Dict()
    flavoredmodes::Vector{Mode} = []
    for mode in mapped
        if !haskey(degenerates, mode)
            degenerates[mode] = 1
        else
            degenerates[mode] += 1
        end
        push!(flavoredmodes, mode|>setattr(:flavor=>degenerates[mode]))
    end
    return Subset(flavoredmodes)
end
export mapmodes

""" A shorthand for `mapmodes`. """
mapmodes(mapper::Function)::Function = input -> mapmodes(mapper, input)

"""
    weightedspatialmap(localisometry::FockMap)::FockMap{RegionFock, RegionFock}

Given `localisometry` with an `outspace` of `FockMap{Region}`, determine the center position of the column 
function by the weights of the eigenfunction and generate a identity map that transforms the `inspace` of 
`localisometry` to include the actual physical attribute of `:r` and `:b`.

### Output
The transformer `FockMap` with `inspace` of `localisometry` and the spatially decorated version as the `outspace`.
"""
function weightedspatialmap(localisometry::FockMap)::FockMap{RegionFock, RegionFock}
    function mapper(mode::Mode)::Mode
        modeisometry::FockMap = localisometry[:, mode]
        absmap::FockMap = modeisometry|>abs
        weights::FockMap = absmap / (absmap|>rep|>collect|>real|>sum)
        pos::Offset = reduce(+, (outmode|>getpos) * (weights[outmode, mode]|>real) for outmode in weights|>getoutspace)
        b::Offset = pos|>basispoint
        r::Offset = pos - b
        # removeattr(:eigenindex) is added since localisometry is coming from a EigenSpectrum.
        return mode|>setattr(:r=>r, :b=>b)|>removeattr(:eigenindex)
    end
    return idmap(localisometry|>getinspace, localisometry|>getinspace|>mapmodes(mapper)|>RegionFock)
end
export weightedspatialmap

function spatialmap(localisometry::FockMap)::FockMap
    function mapper(mode::Mode)::Mode
        pos::Offset = localisometry|>getoutspace|>getregion|>getcenter
        b::Offset = pos|>basispoint
        r::Offset = pos - b
        # removeattr(:eigenindex) is added since localisometry is coming from a EigenSpectrum.
        return mode|>setattr(:r=>r, :b=>b)|>removeattr(:eigenindex)
    end
    return idmap(localisometry|>getinspace, localisometry|>getinspace|>mapmodes(mapper)|>RegionFock)
end
export spatialmap

struct RegionState{Dim} <: Element{FockMap}
    spstates::Dict{Mode, FockMap}
end
export RegionState

Base.:show(io::IO, state::RegionState) = print(io, string("$(typeof(state))(count=$(state |> getinspace |> dimension))"))

Base.:convert(::Type{FockMap}, source::RegionState)::FockMap = reduce(+, spstate for (_, spstate) in spstates)

function regionalrestriction(crystalstate::FockMap, regionfock::RegionFock)::RegionState
    eigenmodes::Subset{Mode} = crystalstate |> getinspace |> unitcellfock |> orderedmodes

    function extractregionstate(mode::Mode)
        rightfourier::FockMap = fourier(crystalstate |> getinspace, mode |> RegionFock)
        leftfourier::FockMap = fourier(crystalstate |> getoutspace, regionfock)
        return leftfourier' * crystalstate * rightfourier
    end

    return Dict(mode => mode |> extractregionstate for mode in eigenmodes) |> RegionState{crystalstate |> getoutspace |> getcrystal |> dimension}
end
export regionalrestriction

"""" Convert a state like `FockMap` with outspace type of `RegionFock` to a `RegionState`. """
function RegionState(localstates::SparseFockMap{RegionFock, <:FockSpace})
    decorated::FockMap = localstates * spatialmap(localstates)
    dim::Integer = decorated|>getoutspace|>getregion|>getspace|>dimension
    return Dict(mode=>decorated[:, mode] for mode in decorated|>getinspace)|>RegionState{dim}
end

function FockMap(regionstate::RegionState)
    statemap::FockMap = reduce(+, state for (_, state) in regionstate.spstates)
    return FockMap(statemap, outspace=statemap|>getoutspace|>RegionFock, inspace=statemap|>getinspace|>RegionFock)
end

getinspace(state::RegionState) = FockSpace(m for (m, _) in state.spstates)
getoutspace(state::RegionState) = state.spstates |> first |> last |> getoutspace

Base.:iterate(state::RegionState, i...) = iterate(state.spstates, i...)
Base.:length(state::RegionState) = state.spstates |> length

function Base.:+(a::RegionState, b::RegionState)
    dim::Integer = a|>getoutspace|>getregion|>getspace|>dimension
    combinedstates::Vector = [a.spstates..., b.spstates...]
    combinedinmodes::Base.Generator = (mode for (mode, _) in combinedstates)
    # This step ensures that any duplicated modes from a and b will be mapped to different :flavor.
    mergedmodes::Subset{Mode} = combinedinmodes|>mapmodes(m -> m)
    mappedstates = (mode=>FockMap(state, inspace=mode|>FockSpace, performpermute=false) for (mode, (_, state)) in zip(mergedmodes, combinedstates))
    return RegionState{dim}(mappedstates|>Dict)
end

"""
    CrystalFockMap(outcrystal::Crystal, incrystal::Crystal, blocks::Dict{Tuple{Momentum, Momentum}, FockMap})

An implementation of `FockMap` designed to pack `FockMap{CrystalFock, CrystalFock}` more efficiently by storing blocks of sub-fockmaps into a `Dict`
indexed by the outspace and inspace `Momentum` of the `FockMap` blocks. Noted that when using this type, the performance for generating `CrystalSpectrum`
will be significantly faster than the `SparseFockMap`, yet for algebric operations it will be slower. This type is used when the generated data for
are blocks of `k -> k'` `FockMap` that requires a `directsum` to finalize the packaging, which is very inefficient when using `SparseFockMap`.

### Input
- `outcrystal`  The crystal of the `outspace` of the `FockMap`.
- `incrystal`   The crystal of the `inspace` of the `FockMap`.
- `blocks`      A dictionary indexed by the `Momentum` of the `outspace` and `inspace` of the `FockMap` blocks, with values of the `FockMap` blocks.
"""
struct CrystalFockMap <: FockMap{CrystalFock, CrystalFock}
    outcrystal::Crystal
    incrystal::Crystal
    blocks::Dict{Tuple{Momentum, Momentum}, FockMap}
end
export CrystalFockMap

function Base.:show(io::IO, fockmap::CrystalFockMap)
    example::FockMap = fockmap.blocks|>first|>last
    incrystalvol::Integer = vol(fockmap.incrystal)
    inspacedim::Integer = incrystalvol * dimension(example|>getinspace)
    inspacestring::String = "CrystalFock(sub=$(incrystalvol), dim=$(inspacedim))"
    outcrystalvol::Integer = vol(fockmap.outcrystal)
    outspacedim::Integer = outcrystalvol * dimension(example|>getoutspace)
    outspacestring::String = "CrystalFock(sub=$(outcrystalvol), dim=$(outspacedim))"
    print(io, string("$(inspacestring) => $(outspacestring)"))
end

""" Generates the input `CrystalFock` of the `CrystalFockMap`, if this requires frequent access, it is recommended to store the result. """
function Zipper.:getinspace(fockmap::CrystalFockMap)::CrystalFock
    basismodes::Subset{Mode} = fockmap.blocks|>first|>last|>getinspace|>removeattr(:k)
    return getcrystalfock(basismodes, fockmap.incrystal)
end

""" Generates the output `CrystalFock` of the `CrystalFockMap`, if this requires frequent access, it is recommended to store the result. """
function Zipper.:getoutspace(fockmap::CrystalFockMap)::CrystalFock
    basismodes::Subset{Mode} = fockmap.blocks|>first|>last|>getoutspace|>removeattr(:k)
    return getcrystalfock(basismodes, fockmap.outcrystal)
end

""" Convert a `CrystalFockMap` into a `SparseFockMap`. """
function Zipper.FockMap(fockmap::CrystalFockMap)::SparseFockMap{CrystalFock, CrystalFock}
    outspace::CrystalFock = fockmap|>getoutspace
    inspace::CrystalFock = fockmap|>getinspace

    function compute(batch)
        data::SparseMatrixCSC{ComplexF64, Int64} = spzeros(outspace|>dimension, inspace|>dimension)
        for (_, block) in batch
            blockoutspace::FockSpace = block|>getoutspace
            outorder::UnitRange = outspace[blockoutspace|>first]:outspace[blockoutspace|>last]
            blockinspace::FockSpace = block|>getinspace
            inorder::UnitRange = inspace[blockinspace|>first]:inspace[blockinspace|>last]
            data[outorder, inorder] += block|>rep
        end
        return data
    end

    batchsize::Integer = (length(fockmap.blocks) / getmaxthreads())|>ceil
    batches = Iterators.partition(fockmap.blocks, batchsize)
    spdata = paralleltasks(
        name="FockMap(::CrystalFockMap)",
        tasks=(()->compute(batch) for batch in batches),
        count=getmaxthreads())|>parallel|>sum

    return FockMap(outspace, inspace, spdata)
end

function Base.:+(a::CrystalFockMap, b::CrystalFockMap)::CrystalFockMap
    function compute(data)
        pair, block = data
        return pair=>(haskey(b.blocks, pair) ? block + b.blocks[pair] : block)
    end

    blocks::Dict{Any, Any} = paralleltasks(
        name="CrystalFockMap + CrystalFockMap",
        tasks=(()->compute(data) for data in a.blocks),
        count=length(a.blocks))|>parallel|>Dict

    watchprogress(desc="CrystalFockMap + CrystalFockMap")
    for (pair, block) in b.blocks
        updateprogress()
        if haskey(blocks, pair) continue end
        blocks[pair] = block
    end
    unwatchprogress()

    return CrystalFockMap(a.outcrystal, a.incrystal, blocks)
end

function Base.:-(target::CrystalFockMap)::CrystalFockMap
    blocks::Dict = Dict(mapping=>(-block) for (mapping, block) in target.blocks)
    return CrystalFockMap(target.outcrystal, target.incrystal, blocks)
end

Base.:-(a::CrystalFockMap, b::CrystalFockMap) = a + (-b)

function Base.:*(left::CrystalFockMap, right::CrystalFockMap)
    # [First]
    # Restructure the right blocks so that they can be keyed by the outspace momentum as multiplication 
    # involves the inspace momentum of the left blocks and the outspace momentum of the right blocks.

    function grouping(blocks)
        grouped = Dict()
        for (outk, ink, block) in blocks
            !haskey(grouped, outk) && (grouped[outk] = [])
            # The value is (ink, block) since ink is needed during the multiplication procedure.
            push!(grouped[outk], (ink, block))
        end
        return grouped
    end

    # Divide the blocks into batches and conquer the grouping process using multithreading.
    batchsize::Integer = (length(right.blocks) / getmaxthreads())|>ceil
    batches = Iterators.partition(((outk, ink, block) for ((outk, ink), block) in right.blocks), batchsize)
    groupedbatches = paralleltasks(
        name="CrystalFockMap * CrystalFockMap",
        tasks=(()->grouping(batch) for batch in batches),
        count=getmaxthreads())|>parallel

    # After batched grouping we have to merge the batches back into a single Dict.
    # We will not append the Vector from all batches since it is computationally expensive, 
    # thus we will just push the entire Vector into the Dict value entry.
    rightblocks = Dict()
    watchprogress(desc="CrystalFockMap * CrystalFockMap")
    for groupedbatch in groupedbatches, (outk, blocks) in groupedbatch
        !haskey(rightblocks, outk) && (rightblocks[outk] = [])
        push!(rightblocks[outk], blocks)
        updateprogress()
    end
    # Then in this step we will convert the Vector{Vector} for each key into an iterator of blocks.
    for (k, blockbatches) in rightblocks
        rightblocks[k] = (block for blockbatch in blockbatches for block in blockbatch)
        updateprogress()
    end
    unwatchprogress()

    # [Second]
    # We can then perform the multiplication procedure using multithreading.
    # Each left block with input subspace keyed by `k` will multiply with all right blocks 
    # with the output subspace keyed by the same `k`.

    function compute(outk, ink, block)
        batch = []
        if !haskey(rightblocks, ink) return batch end
        for (rink, rightblock) in rightblocks[ink]
            push!(batch, (outk, rink)=>block*rightblock)
        end
        return batch
    end

    multiplied = paralleltasks(
        name="CrystalFockMap * CrystalFockMap",
        tasks=(()->compute(outk, ink, block) for ((outk, ink), block) in left.blocks),
        count=length(left.blocks))|>parallel

    # Since the results are batched into a Vector for each batch, we have to flatten it 
    # before merging the result into a single Dict.
    multiplied = (el for batch in multiplied for el in batch)

    # [Third]
    # We will merge the multiplied results using batched multithreading, then we will perform 
    # one final merge to get the final result.

    function mergemultiplied(batch)
        merged = Dict()
        for (index, block) in batch
            haskey(merged, index) ? (merged[index] += block) : (merged[index] = block)
        end
        return merged
    end

    count = 0
    watchprogress(desc="CrystalFockMap * CrystalFockMap")
    for _ in multiplied
        count += 1
        if count % 512 == 0
            updateprogress()
        end
    end
    unwatchprogress()

    batchsize = (count / getmaxthreads())|>ceil
    mergedbatch = paralleltasks(
        name="CrystalFockMap * CrystalFockMap",
        tasks=(()->mergemultiplied(batch) for batch in Iterators.partition(multiplied, batchsize)),
        count=getmaxthreads())|>parallel

    blocks = Dict()
    watchprogress(desc="CrystalFockMap * CrystalFockMap")
    for merged in mergedbatch, (index, block) in merged
        haskey(blocks, index) ? (blocks[index] += block) : (blocks[index] = block)
        updateprogress()
    end
    unwatchprogress()

    return CrystalFockMap(left.outcrystal, right.incrystal, blocks)
end

function Base.:*(fockmap::CrystalFockMap, num::Number)
    blocks::Dict = paralleltasks(
        name="CrystalFockMap * Number",
        tasks=(()->(pair=>(num * block)) for (pair, block) in fockmap.blocks),
        count=length(fockmap.blocks))|>parallel|>Dict
    return CrystalFockMap(fockmap.outcrystal, fockmap.incrystal, blocks)
end

Base.:*(num::Number, fockmap::CrystalFockMap) = fockmap * num

function Base.:abs(fockmap::CrystalFockMap)::CrystalFockMap
    blocks::Dict = paralleltasks(
        name="abs CrystalFockMap",
        tasks=(()->(pair=>(abs(block))) for (pair, block) in fockmap.blocks),
        count=length(fockmap.blocks))|>parallel|>Dict
    return CrystalFockMap(fockmap.outcrystal, fockmap.incrystal, blocks)
end

function Base.:adjoint(fockmap::CrystalFockMap)
    blocks::Dict = paralleltasks(
        name="adjoint",
        tasks=(()->((ink, outk)=>block') for ((outk, ink), block) in fockmap.blocks),
        count=length(fockmap.blocks))|>parallel|>Dict
    return CrystalFockMap(fockmap.incrystal, fockmap.outcrystal, blocks)
end

function Base.:transpose(fockmap::CrystalFockMap)::CrystalFockMap
    blocks::Dict = paralleltasks(
        name="transpose",
        tasks=(()->((ink, outk)=>transpose(block)) for ((outk, ink), block) in fockmap.blocks),
        count=length(fockmap.blocks))|>parallel|>Dict
    return CrystalFockMap(fockmap.incrystal, fockmap.outcrystal, blocks)
end

Zipper.:crystalsubmaps(fockmap::CrystalFockMap)::Base.Generator = (outk=>block for ((outk, ink), block) in fockmap.blocks if outk == ink)

Zipper.:crystalspectrum(fockmap::CrystalFockMap)::CrystalSpectrum = crystalspectrum(fockmap|>crystalsubmaps, crystal=fockmap|>getoutspace|>getcrystal)

function CrystalFockMap(crystalspectrum::CrystalSpectrum)::CrystalFockMap
    function compute(k::Momentum)
        modes::Subset{Mode} = crystalspectrum.eigenmodes[k]
        eigenfock::FockSpace = modes |> FockSpace
        diagonal::FockMap = FockMap(eigenfock, eigenfock, Dict((m, m) => crystalspectrum.eigenvalues[m] |> ComplexF64 for m in modes))
        return (k, k)=>(crystalspectrum.eigenvectors[k] * diagonal * crystalspectrum.eigenvectors[k]')
    end

    blocks::Dict = paralleltasks(
        name="CrystalFockMap from CrystalSpectrum",
        tasks=(()->compute(k) for (k, _) in crystalspectrum.eigenmodes),
        count=crystalspectrum|>getcrystal|>vol)|>parallel|>Dict

    return CrystalFockMap(crystalspectrum|>getcrystal, crystalspectrum|>getcrystal, blocks)
end

"""
    crystaldirectsum(kfockmaps; outcrystal::Crystal, incrystal::Crystal)::CrystalFockMap

If the `FockMap` blocks are direct summed to form a `FockMap{CrystalFock, CrystalFock}`, this will be a more efficient way than using `directsum`.
"""
function crystaldirectsum(kfockmaps; outcrystal::Crystal, incrystal::Crystal)::CrystalFockMap
    blocks::Dict = Dict((commonattr(block|>getoutspace, :k), commonattr(block|>getinspace, :k))=>block for block in kfockmaps)
    return CrystalFockMap(outcrystal, incrystal, blocks)
end
export crystaldirectsum

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

function idmap(fockspace::CrystalFock)
    blocks::Dict = paralleltasks(
        name="idmap",
        tasks=(()->((k, k)=>idmap(getsubspace(fockspace, k))) for k in fockspace|>getcrystal|>brillouinzone),
        count=fockspace|>getcrystal|>vol)|>parallel|>Dict
    return CrystalFockMap(fockspace|>getcrystal, fockspace|>getcrystal, blocks)
end
