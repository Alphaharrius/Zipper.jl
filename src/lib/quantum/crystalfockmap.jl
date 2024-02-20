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
# ◆ CrystalFockMap arithmetics ◆
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
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalFockMap APIs ◆
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
    spdata = paralleltasks(
        name="FockMap(::CrystalFockMap)",
        tasks=(()->compute(batch) for batch in batches),
        count=getmaxthreads())|>parallel|>sum

    return FockMap(outspace, inspace, spdata)
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
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalFockMap extension to other APIs ◆
Zipper.:crystalspectrum(fockmap::CrystalFockMap)::CrystalSpectrum = crystalspectrum(
    fockmap|>crystalsubmaps, crystal=fockmap|>getoutspace|>getcrystal)

function idmap(fockspace::CrystalFock)
    blocks::Dict = paralleltasks(
        name="idmap",
        tasks=(()->((k, k)=>idmap(getsubspace(fockspace, k))) for k in fockspace|>getcrystal|>brillouinzone),
        count=fockspace|>getcrystal|>vol)|>parallel|>Dict
    return CrystalFockMap(fockspace|>getcrystal, fockspace|>getcrystal, blocks)
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