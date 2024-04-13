# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalFockMap definitions ◆
"""
    CrystalFockMap(outcrystal::Crystal, incrystal::Crystal, blocks::Dict{Tuple{Momentum, Momentum}, FockMap})

An implementation of `FockMap` designed to pack `FockMap{CrystalFock, CrystalFock}` more efficiently by storing 
blocks of sub-fockmaps into a `Dict` indexed by the outspace and inspace `Momentum` of the `FockMap` blocks. 
Noted that when using this type, the performance for generating `CrystalSpectrum` will be significantly faster 
than the `SparseFockMap`, yet for algebric operations it will be slower. This type is used when the generated 
data for are blocks of `k -> k'` `FockMap` that requires a `directsum` to finalize the packaging, which is very 
inefficient when using `SparseFockMap`.

### Input
- `outcrystal`  The crystal of the `outspace` of the `FockMap`.
- `incrystal`   The crystal of the `inspace` of the `FockMap`.
- `blocks`      A dictionary indexed by the `Momentum` of the `outspace` and `inspace` of the `FockMap` blocks, 
                with values of the `FockMap` blocks.
"""
struct CrystalFockMap <: FockMap{CrystalFock, CrystalFock}
    outcrystal::Crystal
    incrystal::Crystal
    blocks::Dict{Tuple{Momentum, Momentum}, FockMap}
end
export CrystalFockMap
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockMap interfaces implementation ◆
Base.:convert(::Type{SparseMatrixCSC{ComplexF64, Int64}}, v::CrystalFockMap) = v|>FockMap|>rep

"""
Generates the input `CrystalFock` of the `CrystalFockMap`, if this requires frequent access, 
it is recommended to store the result.
"""
function Zipper.:getinspace(fockmap::CrystalFockMap)::CrystalFock
    basismodes::Subset{Mode} = fockmap.blocks|>first|>last|>getinspace|>removeattr(:k)
    return getcrystalfock(basismodes, fockmap.incrystal)
end

"""
Generates the output `CrystalFock` of the `CrystalFockMap`, if this requires frequent access, 
it is recommended to store the result.
"""
function Zipper.:getoutspace(fockmap::CrystalFockMap)::CrystalFock
    basismodes::Subset{Mode} = fockmap.blocks|>first|>last|>getoutspace|>removeattr(:k)
    return getcrystalfock(basismodes, fockmap.outcrystal)
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalFockMap APIs ◆
"""
    outspacesubmaps(right::CrystalFockMap)::Dict{Momentum, Vector}

Given a `CrystalFockMap`, from the perspective of its `outspace`, get all the blocks indexed by the 
`outspace` momentum in the format of `Dict{Momentum, Iterators.Flatten}` for the `Vector` contains 
element in format of `(inspacemomentum::Momentum, block::FockMap)`. This method can assist in implementing 
multiplication algebra for `CrystalFockMap`.
"""
function outspacesubmaps(right::CrystalFockMap)::Dict
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
        name="outspacesubmaps",
        tasks=(()->grouping(batch) for batch in batches),
        count=getmaxthreads())|>parallel

    function merger(batch)
        merged = Dict()
        for v in batch, (outk, blocks) in v
            !haskey(merged, outk) && (merged[outk] = [])
            append!(merged[outk], blocks)
            updatedivideconquer()
        end
        return merged
    end

    return paralleldivideconquer(merger, groupedbatches, count=getmaxthreads(), desc="outspacesubmaps")
end
export outspacesubmaps
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalFockMap arithmetics ◆
function Base.:+(a::CrystalFockMap, b::CrystalFockMap)::CrystalFockMap
    @assert hassamespan(a|>getoutspace, b|>getoutspace)
    @assert hassamespan(a|>getinspace, b|>getinspace)

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
    @assert hassamespan(left|>getinspace, right|>getoutspace)

    rightblocks::Dict = outspacesubmaps(right)

    # We can then perform the multiplication procedure using multithreading.
    # Each left block with input subspace keyed by `k` will multiply with all right blocks 
    # with the output subspace keyed by the same `k`.

    function compute(outk, ink, block)
        batch = Dict()
        if !haskey(rightblocks, ink) return batch end
        for (rink, rightblock) in rightblocks[ink]
            batch[(outk, rink)] = block*rightblock
        end
        return batch
    end

    multiplied = paralleltasks(
        name="CrystalFockMap * CrystalFockMap",
        tasks=(()->compute(outk, ink, block) for ((outk, ink), block) in left.blocks),
        count=length(left.blocks))|>parallel

    # We will merge the multiplied results using batched multithreading, then we will perform 
    # one final merge to get the final result.

    function mergemultiplied(batch)
        merged = Dict()
        for v in batch, (index, block) in v
            if !iszero(block)
                haskey(merged, index) ? (merged[index] += block) : (merged[index] = block)
            end
            updatedivideconquer()
        end
        return merged
    end

    blocks = paralleldivideconquer(
        mergemultiplied, multiplied, count=left.blocks|>length, desc="CrystalFockMap * CrystalFockMap")

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
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalFockMap APIs ◆
crystalfockmapimpl = CrystalFockMap
crystalfockmap(args...) = crystalfockmapimpl(args...)
export crystalfockmap

crystalsubmaps(fockmap::CrystalFockMap)::Base.Generator = (
    outk=>block for ((outk, ink), block) in fockmap.blocks if outk == ink)
export crystalsubmaps

"""
    crystaldirectsum(kfockmaps; outcrystal::Crystal, incrystal::Crystal)::CrystalFockMap

If the `FockMap` blocks are direct summed to form a `FockMap{CrystalFock, CrystalFock}`, this will be a more 
efficient way than using `directsum`.
"""
function crystaldirectsum(kfockmaps; outcrystal::Crystal, incrystal::Crystal)::CrystalFockMap
    blocks::Dict = Dict(
        (commonattr(block|>getoutspace, :k), commonattr(block|>getinspace, :k))=>block for block in kfockmaps)
    return CrystalFockMap(outcrystal, incrystal, blocks)
end
export crystaldirectsum

function truncatetoregion(hermitian::CrystalFockMap, region::Region)
    regionfock = getregionfock(hermitian|>getoutspace, region)
    restrict = fourier(hermitian|>getoutspace, regionfock)
    localcorrelations = restrict'*hermitian*restrict/(hermitian|>getoutspace|>getcrystal|>vol)
    localindices = Iterators.product(localcorrelations|>getoutspace, localcorrelations|>getoutspace)
    pruningindices = Dict(
        (getattr(to, :r) - getattr(from, :r), from|>removeattr(:r), to|>removeattr(:r))=>(from, to) 
        for (from, to) in localindices)
    pruningindices = (index for (_, index) in pruningindices)
    pruned = extractindices(localcorrelations, pruningindices)
    return restrict * pruned .* restrict'
end
export truncatetoregion
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalFockMap conversions ◆
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
    datas = paralleltasks(
        name="FockMap(::CrystalFockMap)",
        tasks=(()->compute(batch) for batch in batches),
        count=getmaxthreads())|>parallel

    spdata = paralleldivideconquer(
        getconquerer(+), datas, count=getmaxthreads(), desc="FockMap(::CrystalFockMap)")

    return FockMap(outspace, inspace, spdata)
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalFockMap extension to other APIs ◆
function idmap(fockspace::CrystalFock)
    blocks::Dict = paralleltasks(
        name="idmap",
        tasks=(()->((k, k)=>idmap(getsubspace(fockspace, k))) for k in fockspace|>getcrystal|>brillouinzone),
        count=fockspace|>getcrystal|>vol)|>parallel|>Dict
    # Using the current implementation of `CrystalFockMap`.
    return crystalfockmap(fockspace|>getcrystal, fockspace|>getcrystal, blocks)
end

function idmap(outspace::CrystalFock, inspace::CrystalFock)
    @assert hassamespan(outspace, inspace)
    return idmap(outspace)
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalFockMap display ◆
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
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
