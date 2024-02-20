abstract type FockMap{A <: FockSpace, B <: FockSpace} <: Element{SparseMatrixCSC{ComplexF64, Int64}} end
export FockMap

"""
    SparseFockMap(outspace::FockSpace, inspace::FockSpace, rep::SparseMatrixCSC{ComplexF64, Int64})
    SparseFockMap(outspace::FockSpace, inspace::FockSpace, rep::AbstractArray{<:Number})
    SparseFockMap(outspace::FockSpace, inspace::FockSpace, mapping::Dict{Tuple{Mode, Mode}, ComplexF64})
    SparseFockMap(fockmap::FockMap; outspace::FockSpace = fockmap|>getoutspace, inspace::FockSpace = fockmap|>getinspace)

Represents an mapping between two fockspaces with the same span of the underlying Hilbert space.

### Input
- `outspace` The output `FockSpace` of this map. If this object is the multiplier, then this will be the `outspace` of the resulting `FockMap`;
             if this object is the factor, then this must have the same span as the `inspace` of the multiplier.
- `inspace`  The input `FockSpace` of this map. If this object is the multiplier, then must have the same span as the `outspace` of the multiplier;
             if this object is the factor, then this will be the `inspace` of the resulting `FockMap`.
- `rep`      A complex sparse matrix represents the 2-point maps between the elements of the `inspace` & `outspace`.
- `mapping`  The values of the map have to be specified for a distinct 2-point pair, keyed by the pair in `Tuple`.

### Examples
- `SparseFockMap(outspace::FockSpace, inspace::FockSpace, rep::SparseMatrixCSC{ComplexF64, Int64})` is used when every ingredients are precomputed.
- `SparseFockMap(outspace::FockSpace, inspace::FockSpace, rep::AbstractArray{<:Number})` is used when the `rep` is an arbitary array like object.
- `SparseFockMap(outspace::FockSpace, inspace::FockSpace, mapping::Dict{Tuple{Mode, Mode}, ComplexF64})` is used when the values of the map have to be specified
  for a distinct 2-point pair.
- `SparseFockMap(fockmap::FockMap; outspace::FockSpace = fockmap|>getoutspace, inspace::FockSpace = fockmap|>getinspace)` is for using new `outspace` and `inspace`.
"""
struct SparseFockMap{A <: FockSpace, B <: FockSpace} <: FockMap{A, B}
    outspace::A
    inspace::B
    rep::SparseMatrixCSC{ComplexF64, Int64}

    SparseFockMap(outspace::FockSpace, inspace::FockSpace, rep::SparseMatrixCSC{ComplexF64, Int64}) = new{outspace |> typeof, inspace |> typeof}(outspace, inspace, rep)
    SparseFockMap(outspace::FockSpace, inspace::FockSpace, rep::AbstractArray{<:Number}) = new{outspace |> typeof, inspace |> typeof}(outspace, inspace, SparseMatrixCSC{ComplexF64, Int64}(rep))

    function SparseFockMap(outspace::FockSpace, inspace::FockSpace, mapping::Dict{Tuple{Mode, Mode}, T})::SparseFockMap where {T <: Complex}
        rep::SparseMatrixCSC{ComplexF64, Int64} = spzeros(dimension(outspace), dimension(inspace))
        for ((out_mode::Mode, in_mode::Mode), value::ComplexF64) in mapping
            rep[outspace[out_mode], inspace[in_mode]] = value
        end
        return SparseFockMap(outspace, inspace, rep)
    end

    SparseFockMap(fockmap::FockMap; outspace::FockSpace = fockmap|>getoutspace, inspace::FockSpace = fockmap|>getinspace, performpermute::Bool = true) = (
        SparseFockMap(outspace, inspace, performpermute ? permute(fockmap, outspace=outspace, inspace=inspace) |> rep : fockmap |> rep))
end
export SparseFockMap

"""
Using `FockMap` as the entry point of instantiating a default `SparseFockMap` object.
"""
# Using SparseFockMap as the default implementation of FockMap
FockMap(outspace::FockSpace, inspace::FockSpace, rep::SparseMatrixCSC{ComplexF64, Int64}) = SparseFockMap(outspace, inspace, rep)
FockMap(outspace::FockSpace, inspace::FockSpace, rep::AbstractArray{<:Number}) = SparseFockMap(outspace, inspace, rep)
FockMap(outspace::FockSpace, inspace::FockSpace, mapping::Dict{Tuple{Mode, Mode}})::SparseFockMap= SparseFockMap(outspace, inspace, mapping)
FockMap(fockmap::FockMap; outspace::FockSpace = fockmap|>getoutspace, inspace::FockSpace = fockmap|>getinspace, performpermute::Bool = true) = SparseFockMap(fockmap, outspace=outspace, inspace=inspace, performpermute=performpermute)

getoutspace(fockmap::SparseFockMap)::FockSpace = fockmap.outspace
export getoutspace

getinspace(fockmap::SparseFockMap)::FockSpace = fockmap.inspace
export getinspace

Base.:size(fockmap::FockMap)::Tuple{Int64, Int64} = (dimension(fockmap|>getoutspace), dimension(fockmap|>getinspace))

""" Shorthand to update the `inspace` & `outspace` of the `FockMap`. """
Base.:*(fockspace::FockSpace, fockmap::FockMap)::FockMap = FockMap(fockmap, outspace=fockspace)
Base.:*(fockmap::FockMap, fockspace::FockSpace)::FockMap = FockMap(fockmap, inspace=fockspace)

Base.:show(io::IO, fockmap::FockMap) = print(io, string("$(fockmap|>getinspace) => $(fockmap|>getoutspace)"))

Base.:convert(::Type{SparseMatrixCSC{ComplexF64, Int64}}, source::SparseFockMap) = source.rep

"""
    extractindices(fockmap::FockMap, indices)

This function extracts a subset of indices from a `FockMap` object and returns a new `FockMap` object that only includes these indices. 

### Input
- `fockmap::FockMap`: The input `FockMap` object.
- `indices`: The indices to be extracted.

### Output
- A new `FockMap` object that only includes the specified indices.
"""
function extractindices(fockmap::FockMap, indices)
    function getindex(index)::Tuple{Integer, Integer}
        frommode, tomode = index
        return (fockmap|>getoutspace)[tomode], (fockmap|>getinspace)[frommode]
    end

    extracted::SparseMatrixCSC = spzeros(Complex, fockmap|>getoutspace|>dimension, fockmap|>getinspace|>dimension)
    for index in indices
        fromindex, toindex = getindex(index)
        extracted[fromindex, toindex] = (fockmap|>rep)[fromindex, toindex]
    end
    return FockMap(fockmap|>getoutspace, fockmap|>getinspace, extracted)
end
export extractindices

LinearAlgebra.:tr(fockmap::FockMap) = fockmap|>rep|>tr

LinearAlgebra.:det(fockmap::FockMap) = fockmap|>rep|>det

# ===================================================================================================================================================
# Added to support getindex of FockMap objects.
Base.:getindex(fockmap::FockMap, row, col) = restrict(fockmap, (fockmap |> getoutspace)[row] |> FockSpace, (fockmap |> getinspace)[col] |> FockSpace)
Base.:getindex(fockmap::FockMap, row::Mode, col::Mode) = restrict(fockmap, row |> FockSpace, col |> FockSpace)
Base.:getindex(fockmap::FockMap, ::Colon, col) = restrict(fockmap, fockmap |> getoutspace, (fockmap |> getinspace)[col] |> FockSpace)
Base.:getindex(fockmap::FockMap, ::Colon, col::Mode) = restrict(fockmap, fockmap |> getoutspace, col |> FockSpace)
Base.:getindex(fockmap::FockMap, row, ::Colon) = restrict(fockmap, (fockmap |> getoutspace)[row] |> FockSpace, fockmap |> getinspace)
Base.:getindex(fockmap::FockMap, row::Mode, ::Colon) = restrict(fockmap, row |> FockSpace, fockmap |> getinspace)
Base.:getindex(fockmap::FockMap, ::Colon, ::Colon) = fockmap
Base.:getindex(fockmap::FockMap, rowspace::FockSpace, colspace::FockSpace) = restrict(fockmap, rowspace, colspace)
Base.:getindex(fockmap::FockMap, ::Colon, colspace::FockSpace) = columns(fockmap, colspace)
Base.:getindex(fockmap::FockMap, rowspace::FockSpace, ::Colon) = rows(fockmap, rowspace)
# ===================================================================================================================================================

getmodepair(fockmap::FockMap, coords::CartesianIndex)::Pair{Mode, Mode} = (fockmap|>getinspace)[coords[2]] => (fockmap|>getoutspace)[coords[1]]
export getmodepair

"""
    idmap(outspace::FockSpace, inspace::FockSpace)::FockMap

Create an injective `FockMap` from `inspace` to `outspace` of the same dimension, and the mapping pairs are determined by the order of both spaces,
i.e. the n-th element of the `inspace` will be mapped to the n-th element of the outspace.
"""
function idmap(outspace::FockSpace, inspace::FockSpace)::FockMap
    @assert(dimension(outspace) == dimension(inspace))
    FockMap(outspace, inspace, SparseMatrixCSC(Matrix{Float64}(I(dimension(outspace)))))
end
export idmap

idmap(fockspace::FockSpace) = idmap(fockspace, fockspace)

"""
    onesmap(outspace::FockSpace, inspace::FockSpace)::FockMap

Generate a `FockMap` full of `1`s from `inspace` to `outspace`.
"""
onesmap(outspace::FockSpace, inspace::FockSpace)::FockMap = FockMap(outspace, inspace, spzeros(dimension(outspace), dimension(inspace)) .+ 1)
export onesmap

"""
    zerosmap(outspace::FockSpace, inspace::FockSpace)::FockMap

Generate a `FockMap` full of `0`s from `inspace` to `outspace`.
"""
zerosmap(outspace::FockSpace, inspace::FockSpace)::FockMap = FockMap(outspace, inspace, spzeros(dimension(outspace), dimension(inspace)))
export zerosmap

"""
    colmap(inmode::Mode, rowdata::Vector{Pair{Mode, ComplexF64}})::FockMap

Create a column `FockMap` with the specified complex value entries associated with the modes which forms the outspace, with a dimension `1` inspace.

### Input
- `inmode`  The single `Mode` that will form the inspace of the fock map.
- `rowdata` The `Mode` to `ComplexF64` pairs which the modes forms the outspace of the fock map and the complex values as the entries.

### Output
A column `FockMap` with the outspace ordered by the order of `rowdata`.
"""
function colmap(inmode::Mode, rowdata::Vector{Pair{Mode, ComplexF64}})::FockMap
    outfock::FockSpace = FockSpace(Subset(p.first for p in rowdata))
    spmat::SparseMatrixCSC{ComplexF64, Int64} = spzeros(dimension(outfock), 1)
    foreach(n -> spmat[n, 1] += rowdata[n].second, 1:dimension(outfock))
    return FockMap(outfock, FockSpace(Subset(inmode)), spmat)
end
export colmap

"""
    columns(fockmap::FockMap, restrictspace::FockSpace)::FockMap

Restrict the `inspace` of the `fockmap` by a sub-fockspace `restrictspace`.
"""
function columns(fockmap::FockMap, restrictspace::FockSpace)::FockMap
    restrictindices::Vector{Integer} = [(fockmap|>getinspace)[mode] for mode in orderedmodes(restrictspace)]
    spmat::SparseMatrixCSC{ComplexF64, Int64} = sparse(view(rep(fockmap), :, restrictindices))
    return FockMap(fockmap|>getoutspace, restrictspace, spmat)
end
export columns

"""
    rows(fockmap::FockMap, restrictspace::FockSpace)::FockMap

Restrict the `outspace` of the `fockmap` by a sub-fockspace `restrictspace`.
"""
function rows(fockmap::FockMap, restrictspace::FockSpace)::FockMap
    restrictindices::Vector{Integer} = [(fockmap|>getoutspace)[mode] for mode in orderedmodes(restrictspace)]
    spmat::SparseMatrixCSC{ComplexF64, Int64} = sparse(view(rep(fockmap), restrictindices, :))
    return FockMap(restrictspace, fockmap|>getinspace, spmat)
end
export rows

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

"""
    restrict(fockmap::FockMap, outspace::FockSpace, inspace::FockSpace)::FockMap

Restrict the `outspace` & `inspace` of the `fockmap` by a sub-fockspaces `outspace` & `inspace` respectively.
"""
function restrict(fockmap::FockMap, outspace::FockSpace, inspace::FockSpace)::FockMap
    outindices::Vector{Integer} = [(fockmap|>getoutspace)[mode] for mode in orderedmodes(outspace)]
    inindices::Vector{Integer} = [(fockmap|>getinspace)[mode] for mode in orderedmodes(inspace)]
    spmat::SparseMatrixCSC{ComplexF64, Int64} = sparse(view(rep(fockmap), outindices, inindices))
    return FockMap(outspace, inspace, spmat)
end
export restrict

"""
    permute(source::FockMap; outspace::FockSpace=source|>getoutspace, inspace::FockSpace=source|>getinspace)::FockMap

Permute the columns and rows of the representation of the `source` `FockMap` by `outspace` & `inspace` respectively.

### Input
- `source`   The target `FockMap` to be permuted.
- `outspace` A `FockSpace` with the same span as the `outspace` of the `source` `FockMap`.
- `inspace`  A `FockSpace` with the same span as the `inspace` of the `source` `FockMap`.
"""
function permute(source::FockMap; outspace::FockSpace=source|>getoutspace, inspace::FockSpace=source|>getinspace)::FockMap
    row_rule::Vector{Int64} = orderingrule(source|>getoutspace, outspace)
    col_rule::Vector{Int64} = orderingrule(source|>getinspace, inspace)
    return FockMap(outspace, inspace, SparseArrays.permute(rep(source), row_rule, col_rule))
end
export permute

""" Shorthand for creating a function with single parameter `FockMap` to perform `permute`. """
permute(; outspace::FockSpace, inspace::FockSpace)::Function = fockmap::FockMap -> Zipper.permute(fockmap, outspace=outspace, inspace=inspace)

Base.:-(target::FockMap)::FockMap = FockMap(target|>getoutspace, target|>getinspace, -rep(target))
Base.:+(a::FockMap, b::FockMap)::FockMap = fockadd(a, b)
Base.:-(a::FockMap, b::FockMap)::FockMap = a + (-b)

function Base.:*(a::FockMap, b::FockMap)::FockMap
    @assert(hassamespan(a|>getinspace, b|>getoutspace)) # Even if the fockspaces are different, composition works as long as they have same span.
    return FockMap(a|>getoutspace, b|>getinspace, rep(a) * rep(permute(b, outspace=a|>getinspace, inspace=b|>getinspace)))
end

Base.:*(fockmap::FockMap, number::Number)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, rep(fockmap) * number)
Base.:*(number::Number, fockmap::FockMap)::FockMap = fockmap * number

Base.:/(fockmap::FockMap, number::Number)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, rep(fockmap) / number)

Base.:transpose(source::FockMap)::FockMap = FockMap(source|>getinspace, source|>getoutspace, transpose(rep(source)))
""" Corresponds to the Hermitian adjoint. """
Base.:adjoint(source::FockMap)::FockMap = FockMap(source|>getinspace, source|>getoutspace, rep(source)')

"""
Tensor product between two `FockMap` objects, the outspace and inspace of the result `FockMap` is the tensor product 
between the `FockSpace` of the parameter `primary` and `secondary`.
"""
function Base.kron(primary::FockMap, secondary::FockMap)
    outspace::FockSpace = kron(primary|>getoutspace, secondary|>getoutspace)
    inspace::FockSpace = kron(primary|>getinspace, secondary|>getinspace)
    data::SparseMatrixCSC = kron(primary|>rep, secondary|>rep)
    FockMap(outspace, inspace, data)
end

LinearAlgebra.:norm(fockmap::FockMap)::Number = norm(fockmap |> rep)
LinearAlgebra.:normalize(fockmap::FockMap)::FockMap = FockMap(fockmap |> getoutspace, fockmap |> getinspace, fockmap |> rep |> normalize)
Base.:abs(fockmap::FockMap)::FockMap = FockMap(fockmap |> getoutspace, fockmap |> getinspace, map(abs, fockmap |> rep))

""" Shorthand for retrieving the eigenvectors from the `eigspech` function. """
eigvecsh(hermitian::FockMap, attrs::Pair{Symbol}...)::FockMap = eigspech(hermitian, attrs...) |> geteigenvectors
export eigvecsh

""" Shorthand for retrieving the eigenvalues from the `eigspech` function. """
eigvalsh(hermitian::FockMap, attrs::Pair{Symbol}...)::Dict{Mode, Real} = eigspech(hermitian, attrs...) |> geteigenvalues
export eigvalsh

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
Perform addition of two `FockMap`s, if they carries `inspace` and `outspace` of same span, then this is just a sum of representation values; if they carries
non-overlapping `inspace` and `outspace`, this corresponds to a direct sum; if they have overlapping but different span of `inspace` or `outspace`, the result will
be a `FockMap` with `outspace` and `inspace` with 3 subspaces of (!intersect x 2, intersect x 1), with the corresponding representation summed. Noted that this
function is accessed via `a::FockMap + b::FockMap`.
"""
function fockadd(a::FockMap, b::FockMap)::FockMap
    outspacesamespan::Bool = hassamespan(a|>getoutspace, b|>getoutspace)
    inspacesamespan::Bool = hassamespan(a|>getinspace, b|>getinspace)

    if outspacesamespan && inspacesamespan
        return fockaddsamespan(a, b)
    end

    # FockMap: a   FockMap: b
    # -----------  -----------
    # | 11 | 12 |  | 22 | 23 |
    # -----------  -----------
    # | 21 | 22 |  | 32 | 33 |
    # -----------  -----------

    # Addition
    # ----------------
    # | 11 | 12 | -- |
    # ----------------
    # | 21 | 22 | 23 |
    # ----------------
    # | -- | 32 | 33 |
    # ----------------

    outspace2::FockSpace = intersect(a|>getoutspace, b|>getoutspace)
    inspace2::FockSpace = intersect(a|>getinspace, b|>getinspace)

    if outspace2 |> dimension == 0 && inspace2 |> dimension == 0
        return directsum(a, b)
    end

    outspace1::FockSpace = (a|>getoutspace) - (b|>getoutspace)
    outspace3::FockSpace = (b|>getoutspace) - (a|>getoutspace)

    inspace1::FockSpace = (a|>getinspace) - (b|>getinspace)
    inspace3::FockSpace = (b|>getinspace) - (a|>getinspace)

    outspace::FockSpace = union(outspace1, outspace2, outspace3) # Orthogonal
    inspace::FockSpace = union(inspace1, inspace2, inspace3) # Orthogonal

    data::SparseMatrixCSC{ComplexF64, Int64} = spzeros(outspace |> dimension, inspace |> dimension)

    function internaladdition(source::FockMap, os::FockSpace, is::FockSpace)
        if (os |> dimension) * (is |> dimension) == 0
            return
        end
        outmodes::Subset{Mode} = os |> orderedmodes
        inmodes::Subset{Mode} = is |> orderedmodes

        restricted::FockMap = restrict(source, os, is)
        data[outspace[outmodes|>first]:outspace[outmodes|>last], inspace[inmodes|>first]:inspace[inmodes|>last]] += (restricted|>rep)
    end

    internaladdition(a, outspace1, inspace1)
    internaladdition(a, outspace1, inspace2)
    internaladdition(a, outspace2, inspace1)
    internaladdition(a, outspace2, inspace2)
    internaladdition(b, outspace2, inspace2)
    internaladdition(b, outspace2, inspace3)
    internaladdition(b, outspace3, inspace2)
    internaladdition(b, outspace3, inspace3)

    # Prioritize the fockspaces of `a`.
    returnoutspace::FockSpace = hassamespan(outspace, a|>getoutspace) ? a|>getoutspace : hassamespan(outspace, b|>getoutspace) ? b|>getoutspace : outspace
    returninspace::FockSpace = hassamespan(inspace, a|>getinspace) ? a|>getinspace : hassamespan(inspace, b|>getinspace) ? b|>getinspace : inspace

    result::FockMap = FockMap(outspace, inspace, data)
    return FockMap(result, outspace=returnoutspace, inspace=returninspace)
end

""" Addition of two `FockMap` objects with the same `inspace` and `outspace`. """
function fockaddsamespan(a::FockMap, b::FockMap)::FockMap
    data::SparseMatrixCSC{ComplexF64, Int64} = (
        (a |> rep) + (Zipper.permute(b, outspace=a|>getoutspace, inspace=a|>getinspace) |> rep))
    return FockMap(a|>getoutspace, a|>getinspace, data)
end

"""
    directsum(a::FockMap, b::FockMap)::FockMap

Given two `FockMap` objects with orthogonal span of `inspace` and `outspace`, perform direct sum of the two.

### Output
The direct summed `FockMap`, with both `inspace` and `outspace` as the union of the spaces from both `FockMap` objects,
each forming a subspace within the new `inspace` and `outspace`.
"""
function directsum(a::FockMap, b::FockMap)::FockMap
    outspace::FockSpace = (a|>getoutspace) + (b|>getoutspace)
    inspace::FockSpace = (a|>getinspace) + (b|>getinspace)
    data::SparseMatrixCSC{ComplexF64, Int64} = spzeros(outspace |> dimension, inspace |> dimension)
    data[1:(a|>getoutspace |> dimension), 1:(a|>getinspace |> dimension)] += a |> rep
    data[(a|>getoutspace |> dimension) + 1:end, (a|>getinspace |> dimension) + 1:end] += b |> rep
    return FockMap(outspace, inspace, data)
end
export directsum

"""
    directsum(fockmaps)::FockMap

Given a collection of `FockMap` objects, and perform direct sum of the `FockMap` objects.

### Input
- `fockmaps` An iterable of `FockMap` objects.
"""
function directsum(fockmaps)::FockMap
    outspace::FockSpace = Iterators.map(fockmap -> fockmap|>getoutspace, fockmaps) |> fockspaceunion
    inspace::FockSpace = Iterators.map(fockmap -> fockmap|>getinspace, fockmaps) |> fockspaceunion
    data::SparseMatrixCSC{ComplexF64, Int64} = spzeros(outspace |> dimension, inspace |> dimension)
    function filldata(fockmap::FockMap)
        # The procedure beneath assumes that the fockspace elements are concatenated in order during union operations.
        outmodes::Subset{Mode} = fockmap|>getoutspace |> orderedmodes
        outrange::UnitRange = outspace[outmodes|>first]:outspace[outmodes|>last]
        inmodes::Subset{Mode} = fockmap|>getinspace |> orderedmodes
        inrange::UnitRange = inspace[inmodes|>first]:inspace[inmodes|>last]
        data[outrange, inrange] += fockmap |> rep
    end
    foreach(filldata, fockmaps)
    return FockMap(outspace, inspace, data)
end

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

    CrystalSpectrum(
        crystal::Crystal,
        eigenmodes::Dict{Momentum, Subset{Mode}},
        eigenvalues::Dict{Mode, Number},
        eigenvectors::Dict{Momentum, FockMap}) = new{crystal|>dimension}(crystal, eigenmodes, eigenvalues, eigenvectors)
end
export CrystalSpectrum

Base.:show(io::IO, spectrum::CrystalSpectrum) = print(io, string("$(spectrum |> typeof)(entries=$(spectrum.eigenvalues |> length))"))

""" Shorthand to retrieve the unitcell fockspace from a `CrystalSpectrum`. """
function unitcellfock(spectrum::CrystalSpectrum)
    sourcefock::FockSpace = spectrum |> geteigenvectors |> first |> last |> getoutspace
    return sourcefock|>orderedmodes|>removeattr(:k)|>FockSpace
end

"""
    crystalspectrum(momentumfockmaps; crystal::Crystal)::CrystalSpectrum

Given a collection of Hermitian `FockMap` objects each associated with a specific `Momentum` from the brillouin zone, pack into a
`CrystalSpectrum` object.
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
    crystalspectrum(fockmap::FockMap)::CrystalSpectrum

Given a Hermitian `FockMap` with `inspace` and `outspace` of type `CrystalFock` of same span, pack into a `CrystalSpectrum` object.
"""
crystalspectrum(fockmap::FockMap)::CrystalSpectrum = crystalspectrum(fockmap|>crystalsubmaps, crystal=fockmap|>getinspace|>getcrystal)

"""
    linespectrum(spectrum::CrystalSpectrum)::CrystalSpectrum{1}

Since some spectrum might be defined in a crystal that is implicitly on 1-D space, such as a line in the 2D plane, inorder to visualize
the spectrum properly instead of using the plot in the original dimension, we can embed the spectrum into a 1-D space for plotting.
"""
function linespectrum(spectrum::CrystalSpectrum)::CrystalSpectrum{1}
    crystal::Crystal = spectrum|>getcrystal
    nontrivialbasis = Iterators.filter(v -> (v[1]) == 1, zip(crystal|>size, crystal|>getspace|>getbasisvectors))|>collect
    @assert(nontrivialbasis|>collect|>length == 1, "More than 1 non-trivial basis in this crystal.")
    embeddedbasis::Real = nontrivialbasis[1][2]|>norm
    embeddedspace::RealSpace = RealSpace([embeddedbasis][:, :])
    unitcelllength = crystal|>getunitcell|>length
    dummyunitcell::Region = Subset([r / unitcelllength] ∈ embeddedspace for r in (0:unitcelllength - 1))
    embeddedcrystal::Crystal = Crystal(dummyunitcell, [nontrivialbasis[1][1]])
    return CrystalSpectrum(embeddedcrystal, spectrum|>geteigenmodes, spectrum|>geteigenvalues, spectrum|>geteigenvectors)
end
export linespectrum

""" Get the associated crystal of a `CrystalSpectrum` object. """
Zipper.:getcrystal(spectrum::CrystalSpectrum)::Crystal = spectrum.crystal

"""
    geteigenmodes(spectrum::CrystalSpectrum)::Dict{Momentum, Subset{Mode}}

Get the eigenmodes of a `CrystalSpectrum` object, indexed by the `Momentum` of the brillouin zone of the crystal.

### Output
A dictionary with keys of `Momentum` and values of `Subset{Mode}` which contains the momentum indexed unitcell modes within the crystal.
"""
geteigenmodes(spectrum::CrystalSpectrum)::Dict{Momentum, Subset{Mode}} = spectrum.eigenmodes
export geteigenmodes

"""
    geteigenvalues(spectrum::CrystalSpectrum)::Dict{Mode, Number}

Get the eigenvalues of a `CrystalSpectrum` object, indexed by all the momentum indexed `Mode` objects of the `CrystalFock`.
"""
geteigenvalues(spectrum::CrystalSpectrum)::Dict{Mode, Number} = spectrum.eigenvalues
export geteigenvalues

"""
    geteigenvectors(spectrum::CrystalSpectrum)::Dict{Momentum, FockMap}

Get the eigenvectors associated with each indexing `Momentum` within the brillouin zone, each with `inspace` corresponds to the
returned `Subset{Mode}` of the same indexing `Momentum` from `geteigenmodes(spectrum)`.
"""
geteigenvectors(spectrum::CrystalSpectrum)::Dict{Momentum, FockMap} = spectrum.eigenvectors
export geteigenvectors

"""
    Fockmap(spectrum::CrystalSpectrum)::FockMap

Pack the `CrystalSpectrum` into a `FockMap` with `outspace` and `inspace` of type `CrystalFock` of same span, the unitcell fockspace
of the packed `FockMap` should corresponds directly to the `outspace` of individual engenvectors in the `CrystalSpectrum`.
"""
function FockMap(crystalspectrum::CrystalSpectrum)::FockMap
    function momentumfockmap(k::Momentum)
        modes::Subset{Mode} = crystalspectrum.eigenmodes[k]
        eigenfock::FockSpace = modes |> FockSpace
        diagonal::FockMap = FockMap(eigenfock, eigenfock, Dict((m, m) => crystalspectrum.eigenvalues[m] |> ComplexF64 for m in modes))
        return crystalspectrum.eigenvectors[k] * diagonal * crystalspectrum.eigenvectors[k]'
    end
    fockmap::FockMap = directsum(k |> momentumfockmap for (k, _) in crystalspectrum.eigenmodes)
    crystalfock::FockSpace = FockSpace(fockmap|>getinspace, reflected=crystalspectrum.crystal)
    return FockMap(fockmap, inspace=crystalfock, outspace=crystalfock, performpermute=false)
end

""" Packaging the result computed from eigenvalue decomposition. """
struct EigenSpectrum
    eigenvalues::Dict{Mode, Number}
    eigenvectors::FockMap
end
export EigenSpectrum

function geteigenmodes(spectrum::EigenSpectrum)::Subset{Mode}
    return spectrum.eigenvectors|>getinspace |> orderedmodes
end

geteigenvalues(spectrum::EigenSpectrum)::Dict{Mode, Number} = spectrum.eigenvalues

geteigenvectors(spectrum::EigenSpectrum)::FockMap = spectrum.eigenvectors

Base.:show(io::IO, spectrum::EigenSpectrum) = print(io, string("$(spectrum |> typeof)(entries=$(spectrum.eigenvalues |> length))"))

"""
    eigspech(hermitian::FockMap, attrs::Pair{Symbol}...; groupingthreshold::Real = 1e-7)::EigenSpectrum

Perform Hermitian eigenvalue decomposition to find the eigenvalues and eigenvectors simultaneously, the corresponding eigenmodes will have
the attributes of `:eigenindex` which corresponds to the degenerate group associated to a eigenvalue, and in ascending order of the eigenvalues;
`:flavor` which indicates the individual degrees of freedom within the degenerate group; along with the attributes supplied by `attrs`.

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

""" Internal method that generates the `eigenmode => eigenvalue` pairs. """
function digesteigenvalues(H::Type, V::Type, vals, groupingthreshold::Real, attrs::Pair{Symbol}...)::Base.Iterators.Flatten
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

"""
    eigspec(fockmap::FockMap, attrs::Pair{Symbol}...; groupingthreshold::Real = 1e-7)::EigenSpectrum

Perform eigenvalue decomposition to find the eigenvalues and eigenvectors simultaneously, the corresponding eigenmodes will have
the attributes of `:eigenindex` which corresponds to the degenerate group associated to a eigenvalue; `:flavor` which indicates the
individual degrees of freedom within the degenerate group; along with the attributes supplied by `attrs`.

### Input
- `fockmap`             The source of the decomposition.
- `attrs`               Attributes to be inserted to the generated eigenmodes.
- `groupingthreshold`   The threshold for grouping degenerated eigenvalues, since the eigenvalues are complex numbers, this threshold
                        will be applied to the real and imaginary parts separately.
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

Given a spectrum, attempt to group the eigenmodes based on their corresponding eigenvalues with a eigenvalue grouping threshold.

### Output
A generator yielding `Pair{Number, Subset{Mode}}` objects, with the eigenvalues as keys and the corresponding eigenmodes as values.
"""
function groupbyeigenvalues(spectrum; groupingthreshold::Number = 1e-7)::Base.Generator
    denominator::Integer = (1 / groupingthreshold) |> round |> Integer
    actualvalues::Dict{Rational, Number} = Dict(hashablereal(v, denominator) => v for (_, v) in spectrum |> geteigenvalues)
    items::Base.Generator = (hashablereal(v, denominator) => m for (m, v) in spectrum |> geteigenvalues)
    groups::Dict{Rational, Vector{Mode}} = foldl(items; init=Dict{Rational, Vector{Mode}}()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end
    sortedrationals::Vector{Rational} = sort([(groups |> keys)...])
    return (actualvalues[r] => groups[r] |> Subset for r in sortedrationals)
end
export groupbyeigenvalues

function LinearAlgebra.log(fockmap::FockMap)::FockMap
    mat::SparseMatrixCSC = fockmap |> rep |> Matrix |> log |> SparseMatrixCSC
    return FockMap(fockmap|>getoutspace, fockmap|>getinspace, mat)
end

Base.iszero(fockmap::FockMap)::Bool = iszero(fockmap |> rep)

function LinearAlgebra.:svd(fockmap::FockMap)::Tuple{FockMap, Base.Generator, FockMap}
    leftmodes::Subset{Mode} = Subset(Mode([:svdindex => n]) for n in 1:dimension(fockmap|>getinspace))
    rightmodes::Subset{Mode} = Subset(Mode([:svdindex => n]) for n in 1:dimension(fockmap|>getoutspace))
    U, Σ, Vt = fockmap |> rep |> Matrix |> svd
    svdvalues::Base.Generator = (
        (Mode([:svdindex => n]), Mode([:svdindex => n])) => Σ[n] for n in 1:min(leftmodes |> length, rightmodes |> length))
    return FockMap(fockmap|>getoutspace, leftmodes |> FockSpace, U), svdvalues, FockMap(rightmodes |> FockSpace, fockmap|>getinspace, Vt')
end

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
