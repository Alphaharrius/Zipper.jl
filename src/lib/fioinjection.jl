function writesparse(sparse::SparseMatrixCSC; filename::String)
    filepath = joinpath(fiodir(), filename)
    I, J, V = findnz(sparse)
    Vr = [real(v) for v in V]
    Vi = [imag(v) for v in V]
    dataframe = DataFrame(I=I, J=J, Vr=Vr, Vi=Vi)
    @info "Sparse data is written to $filepath"
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
    return sparse(I, J, V)
end

# Defining a more convenient way of serializing a Matrix.
fiolower((RealSpace, :rep), m -> [m[:, i] for i in axes(m, 2)])
fiolower((MomentumSpace, :rep), m -> [m[:, i] for i in axes(m, 2)])
# JSON is parsing Complex into dictionary which is not desired, 
# so we need to convert it as a Vector beforehand.
fiolower((BasisFunction, :rep), v -> [[c|>real, c|>imag] for c in v])
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

# The Dict keys from deserialization are type of String, convert them back to Symbol.
fioconstructor(Mode, d -> Dict(Symbol(k)=>v for (k, v) in d)|>Mode)
# SparseFock does not have a default constructor, so we have to define it manually.
fioconstructor(SparseFock, function (reflected, rep, ordering)
    actualreflected = reflected isa String ? Nothing : reflected
    return SparseFock(rep, Dict(m=>i for (m, i) in ordering), reflected=actualreflected)
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
    readsparse(filename)
    @debug "$filepath -> sparse"
    return readsparse(filename)
end)
# The keys and values of blocks are serialized into Dict, we have to deserialize them back.
fioparser((CrystalFockMap, :blocks), v -> Dict((fioparse(el[1]), fioparse(el[2]))=>fioparse(el[3]) for el in v))

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
