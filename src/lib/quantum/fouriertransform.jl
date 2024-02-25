# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Fourier transform APIs ◆
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

function fourier(
    crystalfock::CrystalFock, regionfock::RegionFock,
    unitcellfockmapping::Dict{Mode, Mode} = mapunitcellfock(crystalfock, regionfock))
    
    momentummatrix::SparseMatrixCSC = crystalfock|>getcrystal|>computemomentummatrix

    momentumhomefock = crystalfock|>unitcellfock
    values::Array = zeros(Complex, momentumhomefock|>length, size(momentummatrix, 2), regionfock|>dimension)
    
    function fillvalues(n, homemode, m, inmode)
        if !haskey(unitcellfockmapping, homemode) || unitcellfockmapping[homemode] != inmode|>removeattr(:r)
            return
        end
        offsetvector::Vector = inmode|>getattr(:r)|>euclidean|>vec
        values[n, :, m] = exp.(-1im * momentummatrix' * offsetvector)
    end

    # Since each (n, m) only corresponds to one entry, thus this is thread-safe.
    paralleltasks(
        name="fourier $(crystalfock|>dimension)×$(regionfock|>dimension)",
        tasks=(
            ()->fillvalues(n, homemode, m, inmode) 
            for ((n, homemode), (m, inmode))
            in Iterators.product(momentumhomefock|>enumerate, regionfock|>enumerate)),
        count=dimension(momentumhomefock)*dimension(regionfock))|>parallel

    function getktransform(n, k)
        outspace = getsubspace(crystalfock, k)
        data = values[:, n, :]
        return k=>FockMap(outspace, regionfock, data)
    end

    blocks = paralleltasks(
        name="fourier $(crystalfock|>dimension)×$(regionfock|>dimension)",
        tasks=(()->getktransform(n, k) for (n, k) in crystalfock|>getcrystal|>brillouinzone|>enumerate),
        count=crystalfock|>getcrystal|>vol)|>parallel|>Dict

    return FourierMap(crystalfock|>getcrystal, regionfock, blocks)
end
export fourier
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FourierMap definition ◆
struct FourierMap{T} <: FockMap{CrystalFock, T}
    crystal::Crystal
    inspace::T
    data::Dict{Momentum, <:SparseFockMap}
end

struct InvFourierMap{T} <: FockMap{T, CrystalFock}
    crystal::Crystal
    outspace::T
    data::Dict{Momentum, <:SparseFockMap}
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FourierMap essentials ◆
Base.:convert(::Type{SparseMatrixCSC{ComplexF64, Int64}}, v::FourierMap) = v|>FockMap|>rep
Base.:convert(::Type{SparseMatrixCSC{ComplexF64, Int64}}, v::InvFourierMap) = v|>FockMap|>rep

Base.:adjoint(::Type{FourierMap{T}}) where T = InvFourierMap{T}
Base.:adjoint(::Type{InvFourierMap{T}}) where T = FourierMap{T}

FourierMapType = Union{FourierMap, InvFourierMap}
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FourierMap conversions ◆
function FockMap(fockmap::FourierMap)::SparseFockMap{CrystalFock, RegionFock}
    outspace::CrystalFock = fockmap|>getoutspace
    inspace::RegionFock = fockmap|>getinspace

    function compute(batch)
        data::SparseMatrixCSC{ComplexF64, Int64} = spzeros(outspace|>dimension, inspace|>dimension)
        for (k, block) in batch
            blockoutspace::FockSpace = getsubspace(outspace, k)
            outorder::UnitRange = outspace[blockoutspace|>first]:outspace[blockoutspace|>last]
            data[outorder, :] += block|>rep
        end
        return data
    end

    batchsize::Integer = length(fockmap.data)/getmaxthreads()|>ceil
    batches = Iterators.partition(fockmap.data, batchsize)
    datas = paralleltasks(
        name="FockMap(::FourierMap)",
        tasks=(()->compute(batch) for batch in batches),
        count=getmaxthreads())|>parallel

    spdata = paralleldivideconquer(sumwithprogress, datas, count=getmaxthreads(), desc="FockMap(::FourierMap)")

    return FockMap(outspace, inspace, spdata)
end

function FockMap(fockmap::InvFourierMap)
    outspace::RegionFock = fockmap|>getoutspace
    inspace::CrystalFock = fockmap|>getinspace

    function compute(batch)
        data::SparseMatrixCSC{ComplexF64, Int64} = spzeros(outspace|>dimension, inspace|>dimension)
        for (k, block) in batch
            blockinspace::FockSpace = getsubspace(inspace, k)
            inorder::UnitRange = inspace[blockinspace|>first]:inspace[blockinspace|>last]
            data[:, inorder] += block|>rep
        end
        return data
    end

    batchsize::Integer = length(fockmap.data)/getmaxthreads()|>ceil
    batches = Iterators.partition(fockmap.data, batchsize)
    datas = paralleltasks(
        name="FockMap(::FourierMap)",
        tasks=(()->compute(batch) for batch in batches),
        count=getmaxthreads())|>parallel

    spdata = paralleldivideconquer(sumwithprogress, datas, count=getmaxthreads(), desc="FockMap(::InvFourierMap)")

    return FockMap(outspace, inspace, spdata)
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockMap interface implementations ◆
getinspace(fockmap::FourierMap) = fockmap.inspace

function getoutspace(fockmap::FourierMap)::CrystalFock
    basismodes::Subset{Mode} = fockmap.data|>first|>last|>getoutspace|>removeattr(:k)
    return getcrystalfock(basismodes, fockmap.crystal)
end

function getinspace(fockmap::InvFourierMap)::CrystalFock
    basismodes::Subset{Mode} = fockmap.data|>first|>last|>getinspace|>removeattr(:k)
    return getcrystalfock(basismodes, fockmap.crystal)
end

getoutspace(fockmap::InvFourierMap) = fockmap.outspace
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FourierMap APIs ◆
getadjspace(fockmap::FourierMap) = fockmap|>getinspace
getadjspace(fockmap::InvFourierMap) = fockmap|>getoutspace
export getadjspace
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FourierMap indexing ◆
Base.:getindex(::FourierMap, a::Any, b::Any) = error("::FourierMap[$(a|>typeof), $(b|>typeof)] is not supported!")
Base.:getindex(fockmap::FourierMap, k::Momentum, ::Colon) = fockmap.data[k]
Base.:getindex(fockmap::InvFourierMap, ::Colon, k::Momentum) = fockmap.data[k]
Base.:getindex(fockmap::FourierMap, subspace::MomentumFock, ::Colon) = fockmap[subspace|>getmomentum, :]
Base.:getindex(fockmap::InvFourierMap, ::Colon, subspace::MomentumFock) = fockmap[:, subspace|>getmomentum]
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FourierMap arithmetics ◆
""" Perform multiplication of `F' * M`. """
function Base.:*(fouriermap::InvFourierMap, fockmap::CrystalFockMap)::InvFourierMap
    rightblocks::Dict = outspacesubmaps(fockmap)

    compute(k::Momentum, block::SparseFockMap) = Dict(rk=>block*rblock for (rk, rblock) in rightblocks[k])

    batches = paralleltasks(
        name="InvFourierMap * CrystalFockMap",
        tasks=(()->compute(k, block) for (k, block) in fouriermap.data),
        count=length(fouriermap.data))|>parallel
    
    function merger(batch)
        merged = Dict()
        for v in batch, (k, block) in v
            if haskey(merged, k)
                merged[k] += block
            else
                merged[k] = block
            end
            updatedivideconquer()
        end
        return merged
    end

    blocks::Dict{Momentum, SparseFockMap} = paralleldivideconquer(
        merger, batches, count=length(fouriermap.data), desc="InvFourierMap * CrystalFockMap")

    return InvFourierMap(fouriermap.crystal, fouriermap|>getoutspace, blocks)
end

Base.:*(::CrystalFockMap, ::FourierMap)::FourierMap = notimplemented()

function Base.:*(fockmap::FourierMapType, v::Number)
    blocks = paralleltasks(
        name="FourierMap * Number",
        tasks=(()->k=>block*v for (k, block) in fockmap.data),
        count=fockmap.data|>length)|>parallel|>Dict
    return typeof(fockmap)(fockmap.crystal, fockmap|>getadjspace, blocks)
end

Base.:*(v::Number, fockmap::FourierMapType) = fockmap * v

function Base.:*(left::InvFourierMap, right::FourierMap)
    multiplied = paralleltasks(
        name="InvFourierMap * FourierMap",
        tasks=(()->block*right.data[k] for (k, block) in left.data),
        count=left.data|>length)|>parallel
    return paralleldivideconquer(
        sumwithprogress, multiplied, count=left.data|>length, desc="InvFourierMap * FourierMap")
end

Base.:/(fockmap::FourierMapType, v::Number) = fockmap * (1/v)

function Base.:*(left::FourierMap{T}, right::SparseFockMap{T, O}) where {T, O}
    multiplied = paralleltasks(
        name="FourierMap * SparseFockMap",
        tasks=(()->k=>block*right for (k, block) in left.data),
        count=left.data|>length)|>parallel|>Dict
    return FourierMap(left.crystal, right|>getoutspace, multiplied)
end

Base.:*(::FourierMap, ::InvFourierMap) = notimplemented()

""" Diagonal composition. """
function Base.broadcasted(::typeof(*), left::FourierMap, right::InvFourierMap)
    @assert(left.crystal == right.crystal)
    multiplied = paralleltasks(
        name="FourierMap * InvFourierMap",
        tasks=(()->(k, k)=>block*right.data[k] for (k, block) in left.data),
        count=left.data|>length)|>parallel|>Dict
    return CrystalFockMap(left.crystal, left.crystal, multiplied)
end

function Base.:transpose(fockmap::FourierMapType)
    compute(k, block) = k=>transpose(block)
    blocks = paralleltasks(
        name="transpose(::$(fockmap|>typeof))",
        tasks=(()->compute(k, block) for (k, block) in fockmap.data),
        count=fockmap.data|>length)|>parallel|>Dict
    return (fockmap|>typeof)'(fockmap.crystal, fockmap|>getadjspace, blocks)
end

function Base.:adjoint(fockmap::FourierMapType)::InvFourierMap
    compute(k, block) = k=>block'
    blocks = paralleltasks(
        name="adjoint(::$(fockmap|>typeof))",
        tasks=(()->compute(k, block) for (k, block) in fockmap.data),
        count=fockmap.data|>length)|>parallel|>Dict
    return (fockmap|>typeof)'(fockmap.crystal, fockmap|>getadjspace, blocks)
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
