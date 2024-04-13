function writesparse(sparse::SparseMatrixCSC; filename::String)
    filepath = joinpath(fiodir(), filename)
    I, J, V = findnz(sparse)
    Vr = [real(v) for v in V]
    Vi = [imag(v) for v in V]
    dataframe = DataFrame(I=I, J=J, Vr=Vr, Vi=Vi)
    return CSV.write(filepath, dataframe)
end

function readsparse(filename::String)
    filepath = joinpath(fiodir(), filename)
    !isfile(filepath) && error("The required sparse datafile $filepath does not exist!")
    dataframe = CSV.read(filepath, DataFrame)
    I = dataframe.I|>Vector
    J = dataframe.J|>Vector
    Vr = dataframe.Vr
    Vi = dataframe.Vi
    V = [Complex(vr, vi) for (vr, vi) in zip(Vr, Vi)]
    return SparseArrays.sparse(I, J, V)
end

# Defining a more convenient way of serializing a Matrix.
fiolower((RealSpace, :rep), m -> [m[:, i] for i in axes(m, 2)])
fiolower((MomentumSpace, :rep), m -> [m[:, i] for i in axes(m, 2)])
# JSON is parsing Complex into dictionary which is not desired, 
# so we need to convert it as a Vector beforehand.
fiolower((BasisFunction, :rep), v -> [[c|>real, c|>imag] for c in v])
# Since the representation of the Subset fully describes the Subset, we can ignore 
# the ordering which can always be generated from the representation.
fiolower((Subset, :orderings), v -> 0)
# The key of korderings is a Momentum, which is not json serializable, 
# so we need to convert it into a vector of key-value pairs.
fiolower((CrystalFock, :korderings), d -> [[k, i] for (k, i) in d])
# Storing matrix data in JSON format is not efficient, so we have to store it separately 
# in a CSV file and store the filename in the JSON file.
fiolower((SparseFockMap, :rep), function (sparse::SparseMatrixCSC)
    filename = "$(fiotargetname()).csv"
    filepath = writesparse(sparse, filename=filename)
    @debug "sparse -> $filepath"
    return Dict(:filename=>filename)
end)
# The blocks of CrystalFockMap is a Dict keyed by (::Momentum, ::Momentum), we have to lower 
# it into a vector to make everything serializable.
fiolower((CrystalFockMap, :blocks), d -> [[ok, ik, v] for ((ok, ik), v) in d])
# The rules attribute is a Matrix form.
fiolower((BoundaryCondition, :rules), m -> [v for v in eachcol(m)])

# The Dict keys from deserialization are type of String, convert them back to Symbol.
fioconstructor(Mode, d -> Dict(Symbol(k)=>v for (k, v) in d)|>Mode)
# The orderings of the Subset is fully described by its representation.
fioconstructor(Subset, (elements, _) -> Subset(elements))
# NormalFock does not have a default constructor, so we have to define it manually.
fioconstructor(NormalFock, function (reflected, rep)
    actualreflected = reflected isa String ? Nothing : reflected
    return NormalFock(rep, reflected=actualreflected)
end)

# Merge the basis vectors back into a single matrix.
fioparser((RealSpace, :rep), v -> hcat(v...))
fioparser((MomentumSpace, :rep), v -> hcat(v...))
# Since we store Complex into a Vector, we have to convert it back to Complex.
fioparser((BasisFunction, :rep), v -> [Complex(el[1], el[2]) for el in v])
# The keys of ordering is serialized into a Dict, we have to deserialize it back into a Momentum.
fioparser((CrystalFock, :korderings), v -> Dict((fioparse(k), i) for (k, i) in v))
# The rep of SparseFockMap is stored in a CSV file, we have to read it back.
fioparser((SparseFockMap, :rep), function (d::Dict)
    filename = d["filename"]
    @debug "$filepath -> sparse"
    return readsparse(filename)
end)
# The keys and values of blocks are serialized into Dict, we have to deserialize them back.
fioparser((CrystalFockMap, :blocks), v -> Dict((fioparse(el[1]), fioparse(el[2]))=>fioparse(el[3]) for el in v))
# Reconstruct the rules Matrix from Vector{Vector}.
fioparser((BoundaryCondition, :rules), v -> hcat(v...))

struct CrystalFockMapBlob <: Element{Any}
    storagemap::FockMap
    nzkeys::Vector
end
export CrystalFockMapBlob

function Base.:convert(::Type{CrystalFockMapBlob}, map::CrystalFockMap)
    nzkeys = [key for key in keys(map.blocks)]
    return CrystalFockMapBlob(map|>FockMap, nzkeys)
end

function Base.:convert(::Type{CrystalFockMap}, blob::CrystalFockMapBlob)
    outhomefock = blob.storagemap|>getoutspace|>unitcellfock
    inhomefock = blob.storagemap|>getinspace|>unitcellfock
    function extract(ok, ik)
        outspace = outhomefock|>setattr(:k=>ok)|>FockSpace
        inspace = inhomefock|>setattr(:k=>ik)|>FockSpace
        return (ok, ik)=>blob.storagemap[outspace, inspace]
    end
    blocks = paralleltasks(
        name="CrystalFockMapBlob -> CrystalFockMap",
        tasks=(()->extract(ok, ik) for (ok, ik) in blob.nzkeys),
        count=length(blob.nzkeys))|>parallel|>Dict
    return CrystalFockMap(blob.storagemap.outspace|>getcrystal, blob.storagemap.inspace|>getcrystal, blocks)
end

fiostoragetype(CrystalFockMap, CrystalFockMapBlob)

struct CrystalDense <: Element{Any}
    outspace::CrystalFock
    inspace::CrystalFock
    chunkcount::Integer
    chunksize::Integer
    dataframe::DataFrame
    nonzeroids::Set{Tuple}
end

function dense2dataframe(dense::Vector)
    function process(cno, data::SparseVector)
        cids, blocks = findnz(data)
        CI = []
        I = []
        J = []
        Vr = []
        Vi = []
        @debug "$cno have $(length(cids)) blocks..."
        for (cid, block) in zip(cids, blocks)
            bI, bJ, bV = findnz(block|>rep)
            length(bI) == 0 && (bI = [1]; bJ = [1]; bV = [0im])
            append!(CI, [cid for _ in eachindex(bI)])
            append!(I, bI)
            append!(J, bJ)
            append!(Vr, [v|>real for v in bV])
            append!(Vi, [v|>imag for v in bV])
        end
        CN = [cno for _ in eachindex(I)]
        return DataFrame(CN=CN, CI=CI, I=I, J=J, Vr=Vr, Vi=Vi)
    end

    dataframes = paralleltasks(
        name="dense2dataframe",
        tasks=(()->process(cno, data) for (cno, data) in enumerate(dense)),
        count=length(dense))|>parallel

    return paralleldivideconquer(getconquerer(vcat), dataframes, length(dense), "dense2dataframe")
end

function dataframe2dense(
    dataframe::DataFrame; 
    outcrystal::Crystal, incrystal::Crystal, outhomefock::NormalFock, inhomefock::NormalFock, chunksize::Integer)

    grouped = groupby(dataframe, :CN)

    function process(cno::Integer, group::SubDataFrame)
        blockframes = groupby(group, :CI)
        chunk = spzeros(SparseFockMap, chunksize)
        for frame in blockframes
            I = frame.I|>Vector
            J = frame.J|>Vector
            V = [Complex(vr, vi) for (vr, vi) in zip(frame.Vr, frame.Vi)]
            data = sparse(I, J, V, outhomefock|>dimension, inhomefock|>dimension)
            cid = frame.CI[1]
            denseindex = getdenseindex(cno, cid, chunksize)
            outk, ink = getdensemomentums(outcrystal, incrystal, denseindex)
            outspace = outhomefock|>setattr(:k=>outk)|>FockSpace
            inspace = inhomefock|>setattr(:k=>ink)|>FockSpace
            chunk[cid] = SparseFockMap(outspace, inspace, data)
        end
        return cno, chunk
    end

    results = paralleltasks(
        name="dataframe2dense",
        tasks=(()->process(group.CN[1], group) for group in grouped),
        count=length(grouped))|>parallel|>collect

    return [v for (_, v) in sort(results, by=first)]
end

function Base.:convert(::Type{CrystalDense}, fockmap::CrystalDenseMap)
    dataframe::DataFrame = dense2dataframe(fockmap.data)
    return CrystalDense(
        fockmap|>getoutspace, fockmap|>getinspace, 
        fockmap.chunkcount, fockmap.chunksize, dataframe, fockmap.nonzeroids)
end

function Base.:convert(::Type{CrystalDenseMap}, dense::CrystalDense)
    data = dataframe2dense(
        dense.dataframe, 
        outcrystal=dense.outspace|>getcrystal, 
        incrystal=dense.inspace|>getcrystal, 
        outhomefock=dense.outspace|>unitcellfock, 
        inhomefock=dense.inspace|>unitcellfock, 
        chunksize=dense.chunksize)
    return CrystalDenseMap(
        dense.outspace, dense.inspace, 
        dense.chunkcount, dense.chunksize, 
        data, dense.nonzeroids)
end

fiostoragetype(CrystalDenseMap, CrystalDense)

fiolower((CrystalDense, :dataframe), function (dataframe::DataFrame)
    filename = "$(fiotargetname()).csv"
    filepath = joinpath(fiodir(), filename)
    CSV.write(filepath, dataframe)
    @debug "dataframe -> $filepath"
    return Dict(:filename=>filename)
end)

# The dataframe is stored in a CSV file.
fioparser((CrystalDense, :dataframe), function (d::Dict)
    filename = d["filename"]
    filepath = joinpath(fiodir(), filename)
    dataframe = CSV.read(filepath, DataFrame)
    @debug "$filepath -> dataframe"
    return dataframe
end)

# The dataframe is stored in a CSV file.
fioparser((CrystalDense, :nonzeroids), function (vec::Vector)
    return Set(tuple(el...) for el in vec)
end)
