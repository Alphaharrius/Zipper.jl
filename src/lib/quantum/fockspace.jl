# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockSpace definition ◆
abstract type FockSpace <: AbstractSpace{Subset} end
export FockSpace
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockSpace interfaces ◆
# --------------------
# Every implementation of FockSpace should implement the following interface 
# methods to ensure that the FockSpace APIs are compatible with the subtype.

# Overloads that makes rep(::FockSpace) works.
Base.:convert(::Type{Subset}, ::FockSpace) = error("Not implemented!")

Zipper.:dimension(::FockSpace) = error("Not implemented!")

orderedmodes(::Any) = error("Not implemented!")
export orderedmodes

getmodes(::FockSpace)::Set{Mode} = error("Not implemented!")
export getmodes

subspacecount(::FockSpace) = error("Not implemented!")
export subspacecount

unitcellfock(::FockSpace) = error("Not implemented!")
export unitcellfock

Base.:getindex(::FockSpace, v) = error("Not implemented!")

Base.:getindex(::FockSpace, range::UnitRange) = error("Not implemented!")

Base.:in(::Mode, ::FockSpace)::Bool = error("Not implemented!")
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockSpace conversions ◆
Base.:convert(::Type{FockSpace}, subset::Subset) = FockSpace(subset)
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockSpace iterators implementation ◆
Base.:iterate(fockspace::FockSpace, i...) = iterate(fockspace|>orderedmodes, i...)
Base.:length(fockspace::FockSpace) = fockspace|>dimension
Base.:lastindex(fockspace::FockSpace) = fockspace|>dimension
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockSpace APIs ◆
"""
    orderedmodes(fockspace::FockSpace)::Subset{Mode}

Returns an ordered `Subset` of modes of `fockspace`.
"""
orderedmodes(fockspace::FockSpace)::Subset{Mode} = Iterators.flatten(rep(fockspace))
export orderedmodes

""" Check whether `sub` is a subspace of `super`. """
issubspace(super::FockSpace, sub::FockSpace)::Bool = issubset(sub|>orderedmodes, super|>orderedmodes)
export issubspace

"""
    fockspaceunion(fockspaces)::FockSpace

Union of fockspaces, the resulting fockspace will have the same span as the union of the 
underlying modes of the fockspaces.

### Input
- `fockspaces` An iterable of fockspaces to be unioned.
"""
fockspaceunion(focks)::FockSpace = subsetunion(fock|>rep for fock in focks)|>FockSpace
export fockspaceunion

"""
    orderingrule(fromspace::FockSpace, tospace::FockSpace)::Vector{Int64}

This method is used when you need to compose two elements with the composing port fockspaces 
having the same span but different orderings.

### Input
- `fromspace` The fockspace as the source of the permutation.
- `tospace` The fockspace as the target of the permutation.

### output
The returned vector that at each mode position `n` of `tospace`, contains the order index 
a mode in `fromspace` to be moved/permuted into slot `n`, the mode ordering of the permuted 
list of modes will matches the ordering of `tospace`.
"""
function orderingrule(fromspace::FockSpace, tospace::FockSpace)::Vector{Int64}
    @assert(hassamespan(fromspace, tospace))
    return [fromspace[mode] for mode in orderedmodes(tospace)]
end
export orderingrule

"""
    sparsegrouping(fockspace::FockSpace, byattrs::Symbol...)::FockSpace

Grouping the modes within the `fockspace` by the attributes specified by `byattrs` and 
return all groups as an iterator of `FockSpace` objects.
"""
function sparsegrouping(fockspace::FockSpace, byattrs::Symbol...)
    getidentifier(mode)::Vector = [mode |> getattr(attr) for attr in byattrs]

    identifiedmodes::Base.Generator = (m |> getidentifier => m for m in fockspace)
    groups::Dict{Vector, Vector} = foldl(identifiedmodes; init=Dict{Vector, Vector}()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end

    return (group|>FockSpace for group in groups|>values)
end
export sparsegrouping

"""
Shorthand for `sparsegrouping` which generates a function, that takes a `FockSpace` 
and returns the sparse grouped fockspace.
"""
sparsegrouping(byattrs::Symbol...)::Function = fockspace -> sparsegrouping(fockspace, byattrs...)
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Quantization APIs ◆
"""
    quantize(position::Offset, flavorcount::Integer)

Quantizes a position into a `FockSpace` object.

### Input
- `position` The position to be quantized. This is an `Offset` object that represents a point in a quantum system.
- `flavorcount` The number of flavors (types) of particles each site.

### Output
- A `FockSpace`` that represents the quantum state of the position in the Fock space.
"""
function quantize(position::Offset, flavorcount::Integer)
    b::Offset = position|>basispoint
    r::Offset = position - b
    return FockSpace((Mode(:r=>r, :b=>b, :flavor=>flavor) for flavor in 1:flavorcount), reflected=r)
end
export quantize

"""
    quantize(region::Region, flavorcount::Integer)::RegionFock

Quantizes a region into a `RegionFock`.

### Input
- `region` The region to be quantized. This is a `Region` object that represents a part of a quantum system.
- `flavorcount` The number of flavors (types) of particles in the system.

### Output
- `RegionFock`: A `RegionFock` object that represents the quantum state of the region in the Fock space.
"""
quantize(region::Region, flavorcount::Integer)::RegionFock = FockSpace(
    fockspaceunion(quantize(position, flavorcount) for position in region),
    reflected=region)
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockSpace logicals ◆
Base.:union(fockspaces::FockSpace...)::FockSpace = fockspaceunion(fockspaces)

Base.:intersect(a::FockSpace, b::FockSpace) = FockSpace(intersect(a |> orderedmodes, b |> orderedmodes))
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockSpace arithmetics ◆
Base.:+(a::FockSpace, b::FockSpace)::FockSpace = union(a, b)

Base.:-(a::FockSpace, b::FockSpace)::FockSpace = FockSpace(setdiff(a |> orderedmodes , b |> orderedmodes))

""" Tensor product between two `FockSpace` objects. """
Base.kron(primary::FockSpace, secondary::FockSpace) = fockspaceunion(
    FockSpace(merge(pmode, smode) for smode in secondary) for pmode in primary)
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockSpace implementation of spatial APIs ◆
"""
    hassamespan(a::FockSpace, b::FockSpace)::Bool

Check if fockspaces of `a` and `b` shares the same set of underlying modes regardless of partitions.
"""
Zipper.:hassamespan(a::FockSpace, b::FockSpace)::Bool = getmodes(a) == getmodes(b)
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ FockSpace display ◆
""" Displays the fock type, subspace count and dimension information of a `FockSpace`. """
function Base.:show(io::IO, fockspace::FockSpace)
    T::Type = typeof(fockspace)
    sub = subspacecount(fockspace)
    dim = dimension(fockspace)
    print(io, string("$T(sub=$sub, dim=$dim)"))
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
