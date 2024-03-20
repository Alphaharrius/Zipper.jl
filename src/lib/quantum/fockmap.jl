# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockMap definition ◆
abstract type FockMap{A <: FockSpace, B <: FockSpace} <: Element{SparseMatrixCSC{ComplexF64, Int64}} end
export FockMap
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockMap interfaces ◆

# This allows rep(fockmap) to work.
Base.:convert(::Type{SparseMatrixCSC{ComplexF64, Int64}}, v::FockMap) = typeof(v)|>notimplemented

getoutspace(::Any) = notimplemented()
getoutspace(v::FockMap) = typeof(v)|>notimplemented
export getoutspace

getinspace(::Any) = notimplemented()
getinspace(v::FockMap) = typeof(v)|>notimplemented
export getinspace
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockMap logicals ◆
Base.iszero(fockmap::FockMap)::Bool = iszero(fockmap |> rep)
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockMap indexing ◆
Base.:getindex(fockmap::FockMap, row, col) = restrict(
    fockmap, (fockmap|>getoutspace)[row]|>FockSpace, (fockmap|>getinspace)[col]|>FockSpace)
Base.:getindex(fockmap::FockMap, row::Mode, col::Mode) = restrict(fockmap, row|>FockSpace, col|>FockSpace)
Base.:getindex(fockmap::FockMap, ::Colon, col) = restrict(
    fockmap, fockmap|>getoutspace, (fockmap|>getinspace)[col]|>FockSpace)
Base.:getindex(fockmap::FockMap, ::Colon, col::Mode) = restrict(fockmap, fockmap|>getoutspace, col|>FockSpace)
Base.:getindex(fockmap::FockMap, row, ::Colon) = restrict(
    fockmap, (fockmap|>getoutspace)[row]|>FockSpace, fockmap|>getinspace)
Base.:getindex(fockmap::FockMap, row::Mode, ::Colon) = restrict(fockmap, row|>FockSpace, fockmap|>getinspace)
Base.:getindex(fockmap::FockMap, ::Colon, ::Colon) = fockmap
Base.:getindex(fockmap::FockMap, rowspace::FockSpace, colspace::FockSpace) = restrict(fockmap, rowspace, colspace)
Base.:getindex(fockmap::FockMap, ::Colon, colspace::FockSpace) = columns(fockmap, colspace)
Base.:getindex(fockmap::FockMap, rowspace::FockSpace, ::Colon) = rows(fockmap, rowspace)
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockMap arithmetics ◆
Base.:+(v::Number, fockmap::FockMap) = idmap(fockmap|>getoutspace, fockmap|>getinspace)*v + fockmap
Base.:-(v::Number, fockmap::FockMap) = v + (-fockmap)

""" Shorthand to update the `outspace` of the `FockMap`. """
Base.:*(fockspace::FockSpace, fockmap::FockMap)::FockMap = FockMap(fockmap, outspace=fockspace)
""" Shorthand to update the `inspace` of the `FockMap`. """
Base.:*(fockmap::FockMap, fockspace::FockSpace)::FockMap = FockMap(fockmap, inspace=fockspace)

LinearAlgebra.:tr(fockmap::FockMap) = fockmap|>rep|>tr
LinearAlgebra.:det(fockmap::FockMap) = fockmap|>rep|>det

""" Negate the `FockMap`. """
Base.:-(v::FockMap)::FockMap = FockMap(v|>getoutspace, v|>getinspace, -rep(v))
""" Adding two `FockMap` aligned with their corresponding `FockSpace`. """
Base.:+(a::FockMap, b::FockMap)::FockMap = fockadd(a, b)
""" Subtracting two `FockMap` aligned with their corresponding `FockSpace`. """
Base.:-(a::FockMap, b::FockMap)::FockMap = a + (-b)

""" Matrix product of two `FockMap` aligned with their corresponding `FockSpace`. """
function Base.:*(a::FockMap, b::FockMap)::FockMap
    # Even if the fockspaces are different, composition works as long as they have same span.
    @assert(hassamespan(a|>getinspace, b|>getoutspace))
    return FockMap(a|>getoutspace, b|>getinspace, rep(a) * rep(permute(b, outspace=a|>getinspace, inspace=b|>getinspace)))
end
    
""" Scalar product of a `FockMap` with a number. """
Base.:*(fockmap::FockMap, number::Number)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, rep(fockmap)*number)
""" Scalar product of a `FockMap` with a number. """
Base.:*(number::Number, fockmap::FockMap)::FockMap = fockmap*number
""" Scalar division of a `FockMap` by a number. """
Base.:/(fockmap::FockMap, number::Number)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, rep(fockmap)/number)

""" Transposing a `FockMap`. """
Base.:transpose(source::FockMap)::FockMap = FockMap(source|>getinspace, source|>getoutspace, transpose(rep(source)))
""" Hermitian adjoint of a `FockMap`. """
Base.:adjoint(source::FockMap)::FockMap = FockMap(source|>getinspace, source|>getoutspace, rep(source)')

""" Get the matrix norm of a `FockMap`. """
LinearAlgebra.:norm(fockmap::FockMap)::Number = norm(fockmap|>rep)
""" Normalize a `FockMap`. """
LinearAlgebra.:normalize(fockmap::FockMap)::FockMap = FockMap(
    fockmap|>getoutspace, fockmap|>getinspace, fockmap|>rep|>normalize)
""" Find the absolute of a `FockMap` and store the result in the real part each index of its representation. """
Base.:abs(fockmap::FockMap)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, map(abs, fockmap|>rep))

""" Find the logarithm of a `FockMap` and store the result in the real part each index of its representation. """
function LinearAlgebra.log(fockmap::FockMap)::FockMap
    mat::SparseMatrixCSC = fockmap |> rep |> Matrix |> log |> SparseMatrixCSC
    return FockMap(fockmap|>getoutspace, fockmap|>getinspace, mat)
end

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
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockMap APIs ◆
Base.:size(fockmap::FockMap)::Tuple{Int64, Int64} = (dimension(fockmap|>getoutspace), dimension(fockmap|>getinspace))

"""
    extractindices(fockmap::FockMap, indices)

This function extracts a subset of indices from a `FockMap` object and returns 
a new `FockMap` object that only includes these indices. 

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

    extracted::SparseMatrixCSC = spzeros(
        Complex, fockmap|>getoutspace|>dimension, fockmap|>getinspace|>dimension)
    for index in indices
        fromindex, toindex = getindex(index)
        extracted[fromindex, toindex] = (fockmap|>rep)[fromindex, toindex]
    end
    return FockMap(fockmap|>getoutspace, fockmap|>getinspace, extracted)
end
export extractindices

getmodepair(fockmap::FockMap, coords::CartesianIndex)::Pair{Mode, Mode} = (
    (fockmap|>getinspace)[coords[2]] => (fockmap|>getoutspace)[coords[1]])
export getmodepair

"""
    idmap(outspace::FockSpace, inspace::FockSpace)::FockMap

Create an injective `FockMap` from `inspace` to `outspace` of the same dimension, and the 
mapping pairs are determined by the order of both spaces, i.e. the n-th element of the `inspace` 
will be mapped to the n-th element of the outspace.
"""
function idmap(outspace::FockSpace, inspace::FockSpace)::FockMap
    @assert(dimension(outspace) == dimension(inspace))
    FockMap(outspace, inspace, SparseMatrixCSC(Matrix{Float64}(I(dimension(outspace)))))
end
export idmap

""" Generate an identity map from a `FockSpace` to itself. """
idmap(fockspace::FockSpace) = idmap(fockspace, fockspace)

"""
    onesmap(outspace::FockSpace, inspace::FockSpace)::FockMap

Generate a `FockMap` full of `1`s from `inspace` to `outspace`.
"""
onesmap(outspace::FockSpace, inspace::FockSpace)::FockMap = FockMap(
    outspace, inspace, spzeros(dimension(outspace), dimension(inspace)) .+ 1)
export onesmap

"""
    zerosmap(outspace::FockSpace, inspace::FockSpace)::FockMap

Generate a `FockMap` full of `0`s from `inspace` to `outspace`.
"""
zerosmap(outspace::FockSpace, inspace::FockSpace)::FockMap = FockMap(
    outspace, inspace, spzeros(dimension(outspace), dimension(inspace)))
export zerosmap

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
    rowrule::Vector{Int64} = orderingrule(source|>getoutspace, outspace)
    colrule::Vector{Int64} = orderingrule(source|>getinspace, inspace)
    return FockMap(outspace, inspace, SparseArrays.permute(rep(source), rowrule, colrule))
end
export permute

""" Shorthand for creating a function with single parameter `FockMap` to perform `permute`. """
permute(; outspace::FockSpace, inspace::FockSpace)::Function = (
    fockmap::FockMap -> Zipper.permute(fockmap, outspace=outspace, inspace=inspace))

"""
Perform addition of two `FockMap`s, if they carries `inspace` and `outspace` of same span, 
then this is just a sum of representation values; if they carries non-overlapping `inspace` 
and `outspace`, this corresponds to a direct sum; if they have overlapping but different span 
of `inspace` or `outspace`, the result will be a `FockMap` with `outspace` and `inspace` with 
3 subspaces of (!intersect x 2, intersect x 1), with the corresponding representation summed. 
Noted that this function is accessed via `a::FockMap + b::FockMap`.
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
        outrange::UnitRange = outspace[outmodes|>first]:outspace[outmodes|>last]
        inrange::UnitRange = inspace[inmodes|>first]:inspace[inmodes|>last]
        data[outrange, inrange] += (restricted|>rep)
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
    returnoutspace::FockSpace = hassamespan(outspace, a|>getoutspace) ? 
        a|>getoutspace : hassamespan(outspace, b|>getoutspace) ? b|>getoutspace : outspace
    returninspace::FockSpace = hassamespan(inspace, a|>getinspace) ? 
        a|>getinspace : hassamespan(inspace, b|>getinspace) ? b|>getinspace : inspace

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
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockMap numerical APIs. ◆
"""
    makezero(fockmap::FockMap, eps::Number = 1e-7)::FockMap

Round all numerical zeros within the `FockMap` to actual zero with a given tolerance `eps`.
"""
makezero(fockmap::FockMap, eps::Number = 1e-7)::FockMap = FockMap(
    fockmap|>getoutspace, fockmap|>getinspace, map(v -> abs(v|>real) < eps && abs(v|>imag) < eps ? 0im : v, fockmap|>rep))
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
    columnspec(fockmap::FockMap)::Vector{Pair{Mode, ComplexF64}}

Extract the column values as a `Mode` to `ComplexF64` pair from a `N×1` `FockMap`, 
this is used when you need to visualize the column spectrum.
"""
function columnspec(fockmap::FockMap)::Vector{Pair{Mode, ComplexF64}}
    @assert(dimension(fockmap|>getinspace) == 1)
    mat::SparseMatrixCSC{ComplexF64, Int64} = rep(fockmap)
    return [outmode => mat[(fockmap|>getoutspace)[outmode], 1] for outmode in orderedmodes(fockmap|>getoutspace)]
end

function fockmapper(mapper::Function, inspace::FockSpace)
    outspace = inspace|>mapmodes(mapper)|>FockSpace
    return idmap(outspace, inspace)
end
export fockmapper
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockMap commutation APIs ◆
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
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockMap display ◆
Base.:show(io::IO, fockmap::FockMap) = print(io, string("$(fockmap|>getinspace) => $(fockmap|>getoutspace)"))
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
