#\begin:CrystalDenseMap definition
struct CrystalDenseMap <: FockMap{CrystalFock, CrystalFock}
    outspace::CrystalFock
    inspace::CrystalFock
    chunkcount::Integer
    chunksize::Integer
    data::Vector{SparseVector{SparseFockMap}}
    nonzeroids::Set{Tuple}
end

export CrystalDenseMap
#\end

#\begin:FockMap interface implementations
Zipper.getoutspace(fockmap::CrystalDenseMap) = fockmap.outspace
Zipper.getinspace(fockmap::CrystalDenseMap) = fockmap.inspace
#\end

#\begin:CrystalDenseMap internal support utilities
function getchunkinfo(outcrystal::Crystal, incrystal::Crystal)
    chunkcount::Integer = maximum((size(outcrystal)..., size(incrystal)...))*4
    datalength::Integer = vol(outcrystal) * vol(incrystal)
    chunksize::Integer = datalength/chunkcount|>round
    return datalength, chunkcount, chunksize
end

preparedense(chunkcount::Integer, chunksize::Integer) = [
    spzeros(SparseFockMap, chunksize) for _ in 1:chunkcount]
preparedenselocks(chunkcount::Integer) = [ReentrantLock() for _ in 1:chunkcount]

function getdenseindex(outcrystal::Crystal, incrystal::Crystal, outk::Momentum, ink::Momentum)
    return (outcrystal[outk]-1) * vol(incrystal) + incrystal[ink]
end

function getdenseindex(chunkno::Integer, chunkid::Integer, chunksize::Integer)
    return (chunkno-1)*chunksize + chunkid
end

function getchunkindices(denseindex::Integer, chunksize::Integer)
    preindex = denseindex-1
    chunkno::Integer = (denseindex/chunksize|>ceil)
    chunkid::Integer = (preindex%chunksize) + 1
    return chunkno, chunkid
end

function getdenseindices(outcrystal::Crystal, incrystal::Crystal, outk::Momentum)
    start = (outcrystal[outk]-1) * vol(incrystal)
    return (i for i in start+1:start+vol(incrystal))
end

function getdensemomentums(outcrystal::Crystal, incrystal::Crystal, denseindex::Integer)
    outindex = ceil(denseindex/vol(incrystal))|>Integer
    inindex = (denseindex-1)%vol(incrystal)+1
    return outcrystal[outindex], incrystal[inindex]
end
#\end

#\begin:CrystalDenseMap arithmetics
function Base.:+(left::Zipper.CrystalDenseMap, right::Zipper.CrystalDenseMap)
    @assert hassamespan(left|>getoutspace, right|>getoutspace)
    @assert hassamespan(left|>getinspace, right|>getinspace)

    results = paralleltasks(
        name="+(::CrystalDenseMap, ::CrystalDenseMap)",
        tasks=(()->(n, l+r) for (n, (l, r)) in zip(left.data, right.data)|>enumerate),
        count=left.chunkcount)|>parallel|>collect

    data = [el for (_, el) in sort(results, by=first)]

    nonzeroids = Set(id for ids in (left.nonzeroids, right.nonzeroids) for id in ids)

    return CrystalDenseMap(
        left|>getoutspace, left|>getinspace, left.chunkcount, left.chunksize, data, nonzeroids)
end

# TODO: This is really bad in performance...
function Base.:*(left::CrystalDenseMap, right::CrystalDenseMap)
    loutspace = left|>getoutspace
    linspace = left|>getinspace
    routspace = right|>getoutspace
    rinspace = right|>getinspace
    @assert hassamespan(linspace, routspace)
    
    _, chunkcount, chunksize = getchunkinfo(loutspace|>getcrystal, rinspace|>getcrystal)
    data = preparedense(chunkcount, chunksize)
    locks = preparedenselocks(chunkcount)

    function process(chunkno, chunkid)
        # Get all the non-zero blocks from the right for this block on the left.
        denseindex = getdenseindex(chunkno, chunkid, left.chunksize)
        outkl, inkl = getdensemomentums(loutspace|>getcrystal, linspace|>getcrystal, denseindex)
        rdenseindices = getdenseindices(routspace|>getcrystal, rinspace|>getcrystal, inkl)
        rchunkindices = Iterators.filter(
            v->last(v)∈right.nonzeroids, 
            (did, getchunkindices(did, right.chunksize)) for did in rdenseindices)
        rchunks = (
            (getdensemomentums(routspace|>getcrystal, rinspace|>getcrystal, did)|>last, cno, cid) 
            for (did, (cno, cid)) in rchunkindices)
        lblock = left.data[chunkno][chunkid]
        products = (
            (getdenseindex(loutspace|>getcrystal, rinspace|>getcrystal, outkl, inkr), lblock*right.data[cno][cid]) 
            for (inkr, cno, cid) in rchunks)
        nzids = []
        for (did, block) in products
            cno, cid = getchunkindices(did, chunksize)
            locks[cno]|>lock
            data[cno][cid] += block
            locks[cno]|>unlock
            push!(nzids, (cno, cid))
        end
        return nzids
    end

    nonzerobatches = paralleltasks(
        name="*(::CrystalDenseMap, ::CrystalDenseMap)",
        tasks=(()->process(chunkno, chunkid) for (chunkno, chunkid) in left.nonzeroids),
        count=length(left.nonzeroids))|>parallel
    nonzeroids = Set(v for batch in nonzerobatches for v in batch)
    
    return CrystalDenseMap(loutspace, rinspace, chunkcount, chunksize, data, nonzeroids)
end
#\end

#\begin:CrystalDenseMap APIs
function CrystalDenseMap(outcrystal::Crystal, incrystal::Crystal, blocks)
    _, chunkcount, chunksize = getchunkinfo(outcrystal, incrystal)
    data::Vector{SparseVector} = preparedense(chunkcount, chunksize)
    locks::Vector{ReentrantLock} = preparedenselocks(chunkcount)

    function insert(outk::Momentum, ink::Momentum, block::SparseFockMap)
        denseindex = getdenseindex(outcrystal, incrystal, outk, ink)
        chunkno, chunkid = getchunkindices(denseindex, chunksize)
        locks[chunkno]|>lock
        data[chunkno][chunkid] = block
        locks[chunkno]|>unlock
        return chunkno, chunkid
    end

    nonzeroids = paralleltasks(
        name="crystaldensemap",
        tasks=(()->insert(outk, ink, block) for ((outk, ink), block) in blocks),
        count=length(blocks))|>parallel|>Set

    refchunkno, refchunkid = nonzeroids|>first
    refblock::SparseFockMap = data[refchunkno][refchunkid]
    outspace = getcrystalfock(refblock|>getoutspace|>removeattr(:k), outcrystal)
    inspace = getcrystalfock(refblock|>getinspace|>removeattr(:k), incrystal)
    return CrystalDenseMap(outspace, inspace, chunkcount, chunksize, data, nonzeroids)
end

CrystalDenseMap(fockmap::CrystalFockMap) = CrystalDenseMap(fockmap.outcrystal, fockmap.incrystal, fockmap.blocks)
#\end

#\begin:CrystalDenseMap conversions
function CrystalFockMap(fockmap::CrystalDenseMap)
    outcrystal::Crystal = fockmap|>getoutspace|>getcrystal
    incrystal::Crystal = fockmap|>getinspace|>getcrystal

    function getblock(chunkno, chunkid)
        block = fockmap.data[chunkno][chunkid]
        denseindex = getdenseindex(chunkno, chunkid, fockmap.chunksize)
        outk, ink = getdensemomentums(outcrystal, incrystal, denseindex)
        return (outk, ink)=>block
    end

    blocks = paralleltasks(
        name="CrystalFockMap(::CrystalDenseMap)",
        tasks=(()->getblock(chunkno, chunkid) for (chunkno, chunkid) in fockmap.nonzeroids),
        count=length(fockmap.nonzeroids))|>parallel|>Dict

    return CrystalFockMap(outcrystal, incrystal, blocks)
end
#\end
