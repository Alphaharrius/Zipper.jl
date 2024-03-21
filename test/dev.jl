using LinearAlgebra
using Zipper, Plots, SparseArrays
plotlyjs()

setmaxthreads(Threads.nthreads())

square = euclidean(RealSpace, 2)
point = [1/2, 1/2] ∈ square
spatialsnappingcalibration([point])

kspace = convert(MomentumSpace, square)

unitcell = Subset(point)
crystal = Crystal(unitcell, [256, 256])
reciprocalhashcalibration(crystal.sizes)

m = quantize(unitcell, 1)|>first

t_n = ComplexF64(-1.)
t_nn = ComplexF64(0.0)
bonds::FockMap = bondmap([
    (m, m|>setattr(:r=>[1, 0]∈square))=>t_n,
    (m, m|>setattr(:r=>[0, 1]∈square))=>t_n,
    (m, m|>setattr(:r=>[1, 1]∈square))=>t_nn,
    (m, m|>setattr(:r=>[1, -1]∈square))=>t_nn,])

@info "Computing energy spectrum..."
energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
H = energyspectrum|>CrystalFockMap

# struct CrystalDenseMap <: FockMap{CrystalFock, CrystalFock}
#     outspace::CrystalFock
#     inspace::CrystalFock
#     chunkcount::Integer
#     chunksize::Integer
#     data::Vector{SparseVector{SparseFockMap}}
#     nonzeroids::Vector{Tuple}
# end

# Zipper.getoutspace(fockmap::CrystalDenseMap) = fockmap.outspace
# Zipper.getinspace(fockmap::CrystalDenseMap) = fockmap.inspace

# Base.getindex(crystal::Crystal, i::Integer) = brillouinzone(crystal)[i]

# function getdenseindex(outcrystal::Crystal, incrystal::Crystal, outk::Momentum, ink::Momentum)
#     return (outcrystal[outk]-1) * vol(incrystal) + incrystal[ink]
# end

# function getdenseindex(chunkno::Integer, chunkid::Integer, chunksize::Integer)
#     return (chunkno-1)*chunksize + chunkid
# end

# function getchunkindices(denseindex::Integer, chunksize::Integer)
#     preindex = denseindex-1
#     chunkno::Integer = (denseindex/chunksize|>ceil)
#     chunkid::Integer = (preindex%chunksize) + 1
#     return chunkno, chunkid
# end

# function getdensemomentums(outcrystal::Crystal, incrystal::Crystal, denseindex::Integer)
#     outindex = ceil(denseindex/vol(incrystal))|>Integer
#     inindex = (denseindex-1)%vol(incrystal)+1
#     return outcrystal[outindex], incrystal[inindex]
# end

# struct NullFockMap <: FockMap{NormalFock, NormalFock} end

# Base.zero(::FockMap) = NullFockMap()
# Base.zero(::Type{<:FockMap}) = NullFockMap()
# Base.:+(::NullFockMap, right::FockMap) = right
# Base.:+(left::FockMap, ::NullFockMap) = left
# Base.:+(::NullFockMap, ::NullFockMap) = NullFockMap()
# Base.iszero(::NullFockMap) = true
# Base.iszero(::FockMap) = false

# Base.:show(io::IO, ::NullFockMap) = print(io, "NullFockMap")

# function CrystalDenseMap(outcrystal::Crystal, incrystal::Crystal, blocks)
#     chunkcount::Integer = (size(outcrystal)..., size(incrystal)...)|>maximum
#     datalength::Integer = vol(outcrystal) * vol(incrystal)
#     chunksize::Integer = datalength/chunkcount|>round
#     data::Vector{SparseVector} = [spzeros(SparseFockMap, chunksize) for _ in 1:chunkcount]
#     locks::Vector{ReentrantLock} = [ReentrantLock() for _ in 1:chunkcount]

#     function insert(outk::Momentum, ink::Momentum, block::SparseFockMap)
#         denseindex = getdenseindex(outcrystal, incrystal, outk, ink)
#         chunkno, chunkid = getchunkindices(denseindex, chunksize)
#         locks[chunkno]|>lock
#         data[chunkno][chunkid] = block
#         locks[chunkno]|>unlock
#         return chunkno, chunkid
#     end

#     nonzeroids = paralleltasks(
#         name="crystaldensemap",
#         tasks=(()->insert(outk, ink, block) for ((outk, ink), block) in blocks),
#         count=length(blocks))|>parallel|>collect

#     refchunkno, refchunkid = nonzeroids|>first
#     refblock::SparseFockMap = data[refchunkno][refchunkid]
#     outspace = getcrystalfock(refblock|>getoutspace|>removeattr(:k), outcrystal)
#     inspace = getcrystalfock(refblock|>getinspace|>removeattr(:k), incrystal)
#     return CrystalDenseMap(outspace, inspace, chunkcount, chunksize, data, nonzeroids)
# end

# CrystalDenseMap(fockmap::Zipper.CrystalFockMap) = CrystalDenseMap(fockmap.outcrystal, fockmap.incrystal, fockmap.blocks)

# function Base.:+(left::Zipper.CrystalDenseMap, right::Zipper.CrystalDenseMap)
#     @assert hassamespan(left|>getoutspace, right|>getoutspace)
#     @assert hassamespan(left|>getinspace, right|>getinspace)

#     results = paralleltasks(
#         name="+(::CrystalDenseMap, ::CrystalDenseMap)",
#         tasks=(()->(n, l+r) for (n, (l, r)) in zip(left.data, right.data)|>enumerate),
#         count=left.chunkcount)|>parallel|>collect

#     data = [el for (_, el) in sort(results, by=first)]

#     nonzeroids = Set(id for ids in (left.nonzeroids, right.nonzeroids) for id in ids)|>collect

#     return CrystalDenseMap(
#         left|>getoutspace, left|>getinspace, left.chunkcount, left.chunksize, data, nonzeroids)
# end

# @time H+H
dH = @time CrystalDenseMap(H)

# function getdenseindices(outcrystal::Crystal, incrystal::Crystal, outk::Momentum)
#     start = (outcrystal[outk]-1) * vol(incrystal)
#     return (i for i in start+1:start+vol(incrystal))
# end

# function Base.:*(left::CrystalDenseMap, right::CrystalDenseMap)
#     loutspace = left|>getoutspace
#     linspace = left|>getinspace
#     routspace = right|>getoutspace
#     rinspace = right|>getinspace
#     @assert hassamespan(linspace, routspace)
    
#     _, chunkcount, chunksize = Zipper.getchunkinfo(loutspace|>getcrystal, rinspace|>getcrystal)
#     data = Zipper.preparedense(chunkcount, chunksize)
#     locks = Zipper.preparedenselocks(chunkcount)

#     function process(chunkno, chunkid)
#         # Get all the non-zero blocks from the right for this block on the left.
#         denseindex = Zipper.getdenseindex(chunkno, chunkid, left.chunksize)
#         outkl, inkl = Zipper.getdensemomentums(loutspace|>getcrystal, linspace|>getcrystal, denseindex)
#         rdenseindices = getdenseindices(routspace|>getcrystal, rinspace|>getcrystal, inkl)
#         rchunkindices = Iterators.filter(
#             v->last(v)∈right.nonzeroids, 
#             (did, Zipper.getchunkindices(did, right.chunksize)) for did in rdenseindices)
#         rchunks = (
#             (Zipper.getdensemomentums(routspace|>getcrystal, rinspace|>getcrystal, did)|>last, cno, cid) 
#             for (did, (cno, cid)) in rchunkindices)
#         lblock = left.data[chunkno][chunkid]
#         products = (
#             (Zipper.getdenseindex(loutspace|>getcrystal, rinspace|>getcrystal, outkl, inkr), lblock*right.data[cno][cid]) 
#             for (inkr, cno, cid) in rchunks)
#         nzids = []
#         for (did, block) in products
#             cno, cid = Zipper.getchunkindices(did, chunksize)
#             locks[cno]|>lock
#             data[cno][cid] += block
#             locks[cno]|>unlock
#             push!(nzids, (cno, cid))
#         end
#         return nzids
#     end

#     nonzerobatches = paralleltasks(
#         name="*(::CrystalDenseMap, ::CrystalDenseMap)",
#         tasks=(()->process(chunkno, chunkid) for (chunkno, chunkid) in left.nonzeroids),
#         count=length(left.nonzeroids))|>parallel
#     nonzeroids = Set(v for batch in nonzerobatches for v in batch)
    
#     return CrystalDenseMap(loutspace, rinspace, chunkcount, chunksize, data, nonzeroids)
# end

did = Zipper.getdenseindex(16, 10658, dH.chunksize)
Zipper.getdensemomentums(crystal, crystal, did)

@time H+H
@time dH+dH
@time CrystalDenseMap(H)

H*H|>crystalspectrum|>visualize
showtaskmeter(true)
ret = dH*dH
ret
ret|>CrystalFockMap|>crystalspectrum|>visualize
ret.nonzeroids
[v for batch in ret for v in batch]
dH.nonzeroids
(79, 781717) in dH.nonzeroids
did = Zipper.getdenseindex(79, 781717, dH.chunksize)
Zipper.getdensemomentums(crystal, crystal, 13280257)
dH.nonzeroids

# dd = convert(Zipper.CrystalDense, dH)
# dH_ = convert(CrystalDenseMap, dd)
# dH_|>CrystalFockMap|>crystalspectrum|>visualize
# fiodir("/Users/alphaharrius/ZERData")
# fiosave(dH, name="dH")
# dH1 = fioload("dH")
# dH1|>CrystalFockMap|>crystalspectrum|>visualize
# @time dH+dH

# findnz(dH.data[1])[2]

# using DataFrames

# function dense2dataframe(dense::Vector)
#     function process(cno, data::SparseVector)
#         cids, blocks = findnz(data)
#         CI = []
#         I = []
#         J = []
#         Vr = []
#         Vi = []
#         @info "$cno have $(length(cids)) blocks..."
#         for (cid, block) in zip(cids, blocks)
#             bI, bJ, bV = findnz(block|>rep)
#             length(bI) == 0 && (bI = [1]; bJ = [1]; bV = [0im])
#             append!(CI, [cid for _ in eachindex(bI)])
#             append!(I, bI)
#             append!(J, bJ)
#             append!(Vr, [v|>real for v in bV])
#             append!(Vi, [v|>imag for v in bV])
#         end
#         CN = [cno for _ in eachindex(I)]
#         return DataFrame(CN=CN, CI=CI, I=I, J=J, Vr=Vr, Vi=Vi)
#     end

#     dataframes = paralleltasks(
#         name="dense2dataframe",
#         tasks=(()->process(cno, data) for (cno, data) in enumerate(dense)),
#         count=length(dense))|>parallel

#     return paralleldivideconquer(getconquerer(vcat), dataframes, length(dense), "dense2dataframe")
# end

# function dataframe2dense(
#     dataframe::DataFrame; 
#     outcrystal::Crystal, incrystal::Crystal, outhomefock::NormalFock, inhomefock::NormalFock, chunksize::Integer)

#     grouped = groupby(dataframe, :CN)

#     function process(cno::Integer, group::SubDataFrame)
#         blockframes = groupby(group, :CI)
#         chunk = spzeros(SparseFockMap, chunksize)
#         for frame in blockframes
#             I = frame.I|>Vector
#             J = frame.J|>Vector
#             V = [Complex(vr, vi) for (vr, vi) in zip(frame.Vr, frame.Vi)]
#             data = sparse(I, J, V, outhomefock|>dimension, inhomefock|>dimension)
#             cid = frame.CI[1]
#             denseindex = getdenseindex(cno, cid, chunksize)
#             outk, ink = getdensemomentums(outcrystal, incrystal, denseindex)
#             outspace = outhomefock|>setattr(:k=>outk)|>FockSpace
#             inspace = inhomefock|>setattr(:k=>ink)|>FockSpace
#             chunk[cid] = SparseFockMap(outspace, inspace, data)
#         end
#         return cno, chunk
#     end

#     results = paralleltasks(
#         name="dataframe2dense",
#         tasks=(()->process(group.CN[1], group) for group in grouped),
#         count=length(grouped))|>parallel|>collect

#     return [v for (_, v) in sort(results, by=first)]
# end

ret = dH.data|>dense2dataframe

homefock = H|>getoutspace|>unitcellfock
ret2 = dataframe2dense(ret, outcrystal=crystal, incrystal=crystal, outhomefock=homefock, inhomefock=homefock, chunksize=dH.chunksize)

ddH = CrystalDenseMap(dH|>getoutspace, dH|>getinspace, dH.chunkcount, dH.chunksize, ret2, dH.nonzeroids)
ddH|>CrystalFockMap|>crystalspectrum|>visualize

F = [group for group in groupby(ret, :CN)][1]
groupby(F, :CI)
ret|>first

dH.data[1]+copy(dH.data[1])
[v for v in [v for v in dH.data[1]] + [v for v in dH.data[1]] if v isa SparseFockMap]

ret = dH+dH
V = [ret.data[cno][cid] for (cno, cid) in ret.nonzeroids]
Set(id for ids in (ret.nonzeroids, ret.nonzeroids) for id in ids)|>collect
[v for v in V if v isa NullFockMap]

# function CrystalFockMap(fockmap::CrystalDenseMap)
#     outcrystal::Crystal = fockmap|>getoutspace|>getcrystal
#     incrystal::Crystal = fockmap|>getinspace|>getcrystal

#     function getblock(chunkno, chunkid)
#         block = fockmap.data[chunkno][chunkid]
#         denseindex = getdenseindex(chunkno, chunkid, fockmap.chunksize)
#         outk, ink = getdensemomentums(outcrystal, incrystal, denseindex)
#         return (outk, ink)=>block
#     end

#     blocks = paralleltasks(
#         name="CrystalFockMap(::CrystalDenseMap)",
#         tasks=(()->getblock(chunkno, chunkid) for (chunkno, chunkid) in fockmap.nonzeroids),
#         count=length(fockmap.nonzeroids))|>parallel|>Dict

#     return CrystalFockMap(outcrystal, incrystal, blocks)
# end

(H - CrystalFockMap(dH))|>crystalspectrum|>visualize
CrystalFockMap(dH)|>crystalspectrum|>visualize
H - CrystalFockMap(dH) |>crystalspectrum|>visualize

cno, cid = getchunkindices(84049920, dH.chunksize)
getdensemomentums(crystal, crystal, 84049920)
