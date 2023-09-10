if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Transformations) include("transformations.jl") end

module Quantum

using LinearAlgebra, SparseArrays, OrderedCollections, Base.Iterators
using ..Spaces, ..Geometries, ..Transformations

export quantized, transformed, symmetrized
export Mode, FockSpace, FockMap
export hasattr, getattr, setattr, removeattr, setorbital, getorbital, quantize, flavorcount, spanoffset
export dimension, crystalof, crystalsubspaces, commonattr, unitcellfock, subspaces, subspacecount, flattensubspaces
export order, orderedmodes, orderingrule, modes, hassamespan, sparsefock, crystalfock, issparse, fockspaceunion
export columns, rows, colsubmaps, rowsubmaps, restrict, eigvecsh, eigvalsh, eigh, fourier, fockaddsamespan, directsum, idmap, onesmap, colmap, columnspec

"""
    Mode(attrs::Dict{Symbol})
    Mode(generator::Base.Generator)
    Mode(input::Base.Iterators.Flatten)
    Mode(datas::Vector{Pair{Symbol, T}}) where {T}

Represents an element in a `FockSpace`, and uniquely identifies a physical mode.

### Attributes to put in `attrs`
- `:offset` stores a `Point` which is the offset in lattice unit, of this mode relative to the associated basis mode.
- `:pos` stores a `Point` which is the unit cell offset, this is associated to the attribute `:flavor`.
- `:flavor` stores an `Integer` that identifies a fermionic freedom at a lattice site.
- `:orbital` Defines which orbital this mode transforms like under a symmetry, for example for 𝐶₃ symmetry and 𝑠 like orbital `Dict(:c3 => :s)`.

### Input
- `attrs` The attributes which uniquely identifies the `Mode` object.
"""
struct Mode <: AbstractSubset{Mode}
    attrs::Dict{Symbol}

    Mode(attrs::Dict{Symbol}) = new(Dict((:orbital => swave, attrs...)))

    Mode(generator::Base.Generator) = Mode(Dict(generator))
    Mode(input::Base.Iterators.Flatten) = Mode(Dict(input))

    Mode(attrs::Vector{Pair{Symbol, T}}) where {T} = Mode(Dict(attrs...))
end

""" Display the number of attributes that identifies this `Mode`. """
Base.:show(io::IO, mode::Mode) = print(io, string("$(typeof(mode))$(tuple(keys(mode.attrs)...))"))

""" Allow for offset the `:offset` attribute of a mode. """
Base.:+(mode::Mode, offset::Point)::Mode = setattr(mode, :offset => getattr(mode, :offset) + offset)
Base.:-(mode::Mode, offset::Point)::Mode = mode + (-offset)

"""
    spaceof(mode::Mode)

The space of a `Mode` comes from the physical quantities its defined on, such as `:offset` and `:pos`, if none of those are defined,
it will be `euclidean(RealSpace, 1)` as the mode and it's siblings can always be parameterized by a scalar.

### Output
The space of the attribute `:offset` if available, fall back to `:pos` if otherwise; returns a Euclidean space of dimension `1` if both `:offset` & `:pos` is missing.
"""
function Spaces.spaceof(mode::Mode)
    # :offset have a higher priority in determining the space of the mode.
    if hasattr(mode, :offset) return spaceof(getattr(mode, :offset)) end
    if hasattr(mode, :pos) return spaceof(getattr(mode, :pos)) end
    # If the mode does not based on any physical position or quantities for associating with any space, then it will simply lives
    # in a 1D euclidean space as the mode and it's siblings can always be parameterized by a scalar.
    return euclidean(RealSpace, 1)
end

Spaces.Subset(modes::Mode...) = Subset(m for m in modes)

"""
    Spaces.pos(mode::Mode)::Point

Get the actual position of the mode, this method only works when `:offset` and `:pos` are defined in the same space.
"""
Spaces.pos(mode::Mode)::Point = convert(Point, mode)

"""
    hasattr(mode::Mode, key::Symbol)::Bool

Check if the `mode` has the attribute identified by `key`.
"""
hasattr(mode::Mode, key::Symbol)::Bool = haskey(mode.attrs, key)

"""
    getattr(mode::Mode, key::Symbol)

Retrieve the attribute value identified by `key` from `mode`.
"""
getattr(mode::Mode, key::Symbol) = mode.attrs[key]

"""
    getattr(key::Symbol)

Shorthand of `getattr(v, key::Symbol)` with the pipe operator `|>`.

### Examples
The line `mode |> getattr(:flavor)` is equal to `getattr(mode, :flavor)`.
"""
getattr(key::Symbol) = v -> getattr(v, key)

"""
    removeattr(mode::Mode, keys::Symbol...)::Mode

Create a **copy** of `mode` **without** the attributes identified by `keys`.

### Examples
- To remove the attribute of `:offset` and `:pos`, we use `removeattr(mode, :offset, :pos)`.
"""
removeattr(mode::Mode, keys::Symbol...)::Mode = Mode(Dict(filter(p -> !(p.first ∈ keys), mode.attrs)))

"""
    removeattr(keys::Symbol...)

Shorthand of `removeattr(v, keys::Symbol...)::Mode` with the pipe operator `|>`.

### Examples
The line `mode |> removeattr(:flavor, :index)` is equal to `removeattr(mode, :flavor, :index)`.
"""
removeattr(keys::Symbol...) = v -> removeattr(v, keys...)

"""
    removeattr(modes::Subset{Mode}, keys::Symbol...)::Subset{Mode}

Create a **copy** of every `Mode` of `modes` **without** the attributes identified by `keys`, the resulting `Subset` might not have the
same length as the input `modes` as some `Mode` might be **condensed** into a single one after some unique identifier attributes is removed.

### Examples
To remove the attribute of `:offset` and `:pos`, we use `removeattr(modes, :offset, :pos)`.
"""
removeattr(modes::Subset{Mode}, keys::Symbol...)::Subset{Mode} = Subset(OrderedSet{Mode}(removeattr(mode, keys...) for mode in modes))

"""
    setattr(mode::Mode, attrs::Pair{Symbol}...)::Mode

Create a **copy** of `mode` with the new attributes identified by `keys` added to the current attributes, if the attribute exists in `mode`, the current
record will be overwritten.

### Examples
- To add `:offset` and `:flavor`, we use `newmode = setattr(mode, :offset => point, :flavor => 1)`
"""
setattr(mode::Mode, attrs::Pair{Symbol}...)::Mode = Mode(Dict(mode.attrs..., attrs...))

"""
    setattr(attrs::Pair{Symbol}...)

Shorthand of `setattr(v, attrs::Pair{Symbol}...)::Mode` with the pipe operator `|>`.

### Examples
The line `mode |> setattr(:flavor => 1)` is equal to `setattr(mode, :flavor => 1)`.
"""
setattr(attrs::Pair{Symbol}...) = v -> setattr(v, attrs...)

"""
    setattr(subset::Subset{Mode}, attrs::Pair{Symbol}...)::Subset{Mode}

Create a **copy** of `modes` with the new attributes identified by `keys` added to the current attributes for each mode in `modes`, if the attribute
exists in the modes, the current record will be overwritten.

### Examples
- To add `:offset` and `:flavor`, we use `newmodes = setattr(modes, :offset => point, :flavor => 1)`
"""
setattr(subset::Subset{Mode}, attrs::Pair{Symbol}...)::Subset{Mode} = Subset(setattr(mode, attrs...) for mode in subset)

getorbital(mode::Mode, default::BasisFunction)::BasisFunction = hasattr(mode, :orbital) ? getattr(mode, :orbital) : default

getorbital(default::BasisFunction = swave) = mode -> getorbital(mode, default)

setorbital(mode::Mode, basis::BasisFunction)::Mode = setattr(mode, :orbital => basis)

setorbital(basis::BasisFunction) = mode -> setorbital(mode, basis)

"""
    spanoffset(basismodes::Subset{Mode}, points::Subset{<: Point})::Subset{Mode}

Given a set of `basismodes`, and the generator `points`, span the basis modes to the generator `points` with attribute `:offset`, the primary ordering will be
the ordering of `points`, then follows the ordering of `basismodes`.
"""
spanoffset(basismodes::Subset{Mode}, points::Subset{<: Point})::Subset{Mode} = Subset(setattr(mode, :offset => point) for point in points for mode in basismodes)

"""
    flavorcount(basismodes::Subset{Mode})::Integer

Determines the number of flavor with the information provided within `basismodes`, if the provided parameter is not actually the basismodes, the return value
is just the number of distinct attribute `:flavor` within the given parameter set.
"""
flavorcount(basismodes::Subset{Mode})::Integer = trunc(Integer, length(basismodes)) / length(removeattr(basismodes, :flavor))

"""
    sparsefock(basismodes::Subset{Mode}, points::Subset{<: Point})::FockSpace

Given a set of `basismodes`, and the generator `points`, span the basis modes to the generator `points` with attribute `:offset` and form a `FockSpace`. Not to
be mistakened with `spanoffset`, this method will partition the modes by the generator points.
Noted that the ordering of the partitions will follow the ordering of `points`, and the ordering within each partition will follow the ordering of `basismodes`.
"""
function sparsefock(basismodes::Subset{Mode}, points::Subset{<: Point})::FockSpace
    partitions::Vector{Subset{Mode}} = [setattr(basismodes, :offset => point) for point in points]
    modes::Subset{Mode} = spanoffset(basismodes, points)
    orderings::Dict{Mode, Integer} = Dict(mode => index for (index, mode) in enumerate(modes))
    return FockSpace(Subset(partitions::Vector{Subset{Mode}}), orderings)
end

"""
    crystalfock(basismodes::Subset{Mode}, crystal::Crystal)::FockSpace

A short hand to build the crystal fockspace, which is the fockspace containing all modes spanned from `basismodes` by the brillouin zone of the `crystal`.
"""
crystalfock(basismodes::Subset{Mode}, crystal::Crystal)::FockSpace = FockSpace(sparsefock(basismodes, brillouinzone(crystal)), reflected=crystal)

""" By this conversion, one can obtain the actual position of the mode, this method only works when `:offset` and `:pos` are defined in the same space. """
Base.:convert(::Type{Point}, source::Mode)::Point = getattr(source, :offset) + getattr(source, :pos)

Base.:(==)(a::Mode, b::Mode)::Bool = a.attrs == b.attrs

# ==================================================
# Mode can be used as keys in dictionaries and sets.
Base.:hash(mode::Mode)::UInt = hash(mode.attrs)
Base.:isequal(a::Mode, b::Mode) = a == b
# ==================================================

"""
    FockSpace(subsets::Subset{Subset{Mode}}, ordering::Dict{Mode, <: Integer}; reflected=Nothing)
    FockSpace(subset::Subset{Mode}; reflected=Nothing)
    FockSpace(fockspace::FockSpace; reflected=Nothing)
    FockSpace(mode::Mode) = FockSpace(Subset(mode); reflected=Nothing)

A collection of `Modes` or `Mode` partitions in the case of sparse fockspace, implicit ordering of underlying modes is assumed.

The structure of fockspace is assume to hold partitions, a fockspace that **does not** have partition is assumed to have a single partition containing
all underlying modes. The design can be seen with the type of the representation `Subset{Subset{Mode}}`, which the higher order `Subset` holds partitions
represented by `Subset{Mode}`.

The `reflected` attribute is used to store the object this fockspace is reflected to, such as `Crystal` for a crystal fockspace, by default it will be `Nothing`.

### Examples
- `FockSpace(subsets::Subset{Subset{Mode}}, ordering::Dict{Mode, <: Integer}; reflected=Nothing)` is used when all the components of the `FockSpace`
  is already constructed prior instantiation.
- `FockSpace(subset::Subset{Mode}; reflected=Nothing)` is the normal use case to convert a `Subset{Mode}` into `FockSpace`.
- `FockSpace(fockspace::FockSpace; reflected=Nothing)` is used to set the `reflected` attribute.
- `FockSpace(mode::Mode; reflected=Nothing)` is used to create a fockspace with a single mode.
"""
struct FockSpace{T} <: AbstractSpace{Subset{Subset{Mode}}}
    reflected::T
    rep::Subset{Subset{Mode}}
    ordering::Dict{Mode, Integer}

    FockSpace(subsets::Subset{Subset{Mode}}, ordering::Dict{Mode, T}; reflected=Nothing) where {T <: Integer} = new{typeof(reflected)}(reflected, subsets, ordering)
    FockSpace(subset::Subset{Mode}; reflected=Nothing) = FockSpace(
        Subset(subset),
        Dict(mode => order for (order, mode) in enumerate(subset)),
        reflected=reflected)
    FockSpace(fockspace::FockSpace; reflected=Nothing) = FockSpace(rep(fockspace), fockspace.ordering, reflected=reflected)
    FockSpace(mode::Mode; reflected=Nothing) = FockSpace(Subset(mode), reflected=reflected)
end

CrystalFock = FockSpace{Crystal}
export CrystalFock

""" Displays the fock type, subspace count and dimension information of a `FockSpace`. """
Base.:show(io::IO, fockspace::FockSpace) = print(io, string("$(typeof(fockspace))(sub=$(fockspace |> subspacecount), dim=$(fockspace |> dimension))"))

function fockspaceunion(fockspaces)::FockSpace
    subspaces::Subset{Subset{Mode}} = subsetunion(fockspace |> rep for fockspace in fockspaces if fockspace |> dimension != 0)
    modes::Subset{Mode} = subspaces |> Iterators.flatten
    ordering::Dict{Mode, Integer} = Dict(mode => index for (index, mode) in enumerate(modes))
    return FockSpace(subspaces::Subset{Subset{Mode}}, ordering)
end

Base.:union(fockspaces::FockSpace...)::FockSpace = fockspaceunion(fockspaces)

Base.:+(a::FockSpace, b::FockSpace)::FockSpace = union(a, b)

Base.:-(a::FockSpace, b::FockSpace)::FockSpace = FockSpace(setdiff(a |> orderedmodes , b |> orderedmodes))

Base.:intersect(a::FockSpace, b::FockSpace) = FockSpace(intersect(a |> orderedmodes, b |> orderedmodes))

Base.:iterate(fock_space::FockSpace, i...) = iterate(flatten(rep(fock_space)), i...)

"""
    Spaces.:dimension(fockspace::FockSpace)

Returns the number of unique member modes within the `fockspace`, each of those represents a vector from the Hilbert space.
"""
Spaces.:dimension(fockspace::FockSpace) = length(fockspace.ordering) # This is a short cut to retrieve the length, which is the dimension.

"""
    crystalof(crystalfock::FockSpace{Crystal})::Crystal

Shorthand for retrieving the `Crystal` of a `FockSpace{Crystal}`.
"""
crystalof(crystalfock::FockSpace{Crystal})::Crystal = crystalfock.reflected

crystalsubsets(crystalfock::FockSpace{Crystal})::Dict{Momentum, Subset{Mode}} = Dict(commonattr(subspace, :offset) => subspace for subspace in crystalfock |> rep)
export crystalsubsets

"""
    crystalsubspaces(crystalfock::FockSpace{Crystal})::Dict{Momentum, FockSpace}

Retrieve mappings from the crystal momentums to the corresponding fockspaces.
"""
crystalsubspaces(crystalfock::FockSpace{Crystal})::Dict{Momentum, FockSpace} = Dict(commonattr(subspace, :offset) => subspace |> FockSpace for subspace in crystalfock |> rep)

"""
    commonattr(subset::Subset{Mode}, key::Symbol)

Retrieve the common attribute associated to `key` of all the child modes of the `fockspace`, and throws assertion error if the attribute
is not unique within the `fockspace`.
"""
function commonattr(subset::Subset{Mode}, key::Symbol)
    set::Set = Set()
    foreach(m -> push!(set, getattr(m, key)), subset)
    @assert(length(set) == 1, "The modes in this fockspace does not share the same attr `$(key)`!")
    return first(set)
end

"""
    unitcellfock(crystalfock::FockSpace{Crystal})::FockSpace

Retrieve the unit cell fockspace of the system from a `FockSpace{Crystal}`, positioned at the origin of the parent `AffineSpace`.
"""
function unitcellfock(crystalfock::FockSpace{Crystal})::FockSpace
    firstpartition::Subset{Mode} = crystalfock |> rep |> first
    originpoint::Point = firstpartition |> first |> getattr(:pos) |> spaceof |> origin
    return FockSpace(firstpartition |> setattr(:offset => originpoint))
end

"""
    subspaces(fockspace::FockSpace)::Base.Generator

Retrieve the sub-fockspaces of the `fockspace`.
"""
subspaces(fockspace::FockSpace)::Base.Generator = (FockSpace(partition) for partition in rep(fockspace))

"""
    subspacecount(fockspace::FockSpace)::Integer

Get the number of sub-fockspaces of this `fockspace`.
"""
subspacecount(fockspace::FockSpace)::Integer = fockspace |> rep |> length

"""
    flattensubspaces(fockspace::FockSpace)::FockSpace

Merge all subspaces within the `fockspace`.
"""
flattensubspaces(fockspace::FockSpace)::FockSpace = FockSpace(Subset(mode for mode in orderedmodes(fockspace)))

"""
    order(fockspace::FockSpace, mode::Mode)::Int64

Since the fockspace have implicit ordering, this function returns the order index of the `mode` in the `fockspace`.
"""
order(fockspace::FockSpace, mode::Mode)::Int64 = fockspace.ordering[mode]

"""
    modes(fockspace::FockSpace)::Set{Mode}

Returns an unordered set of modes of `fockspace`, this is a more efficient way to retrieve the underlying modes for a fockspace with
more than one partitions.
"""
modes(fockspace::FockSpace)::Set{Mode} = Set(keys(fockspace.ordering)) # This is the most efficient way to get all distinct modes.


"""
    orderedmodes(fockspace::FockSpace)::Subset{Mode}

Returns an ordered `Subset` of modes of `fockspace`.
"""
orderedmodes(fockspace::FockSpace)::Subset{Mode} = flatten(rep(fockspace))

"""
    orderingrule(fromspace::FockSpace, tospace::FockSpace)::Vector{Int64}

This method is used when you need to compose two elements with the composing port fockspaces having the same span but different orderings.

### Input
- `fromspace` The fockspace as the source of the permutation.
- `tospace` The fockspace as the target of the permutation.

### output
The returned vector that at each mode position `n` of `tospace`, contains the order index a mode in `fromspace` to be moved/permuted into slot `n`,
the mode ordering of the permuted list of modes will matches the ordering of `tospace`.
"""
function orderingrule(fromspace::FockSpace, tospace::FockSpace)::Vector{Int64}
    @assert(hassamespan(fromspace, tospace))
    return [fromspace.ordering[mode] for mode in orderedmodes(tospace)]
end

"""
    hassamespan(a::FockSpace, b::FockSpace)::Bool

Check if fockspaces of `a` and `b` shares the same set of underlying modes regardless of partitions.
"""
hassamespan(a::FockSpace, b::FockSpace)::Bool = modes(a) == modes(b)

"""
    issparse(fockspace::FockSpace)::Bool

Check if `fockspace` is a sparse fockspace.
"""
issparse(fockspace::FockSpace)::Bool = subspacecount(fockspace) > 1

""" Check if fockspaces of `a` and `b` has the exact same structure. """
Base.:(==)(a::FockSpace, b::FockSpace)::Bool = rep(a) == rep(b)

# ======================================================================================
# Overloads that makes rep(::FockSpace) works.
Base.:convert(::Type{Subset}, source::FockSpace) = source.rep
Base.:convert(::Type{Subset{Subset}}, source::FockSpace) = convert(Subset, source)
Base.:convert(::Type{Subset{Subset{Mode}}}, source::FockSpace) = convert(Subset, source)
# ======================================================================================

Base.:convert(::Type{FockSpace}, source::Subset{Mode}) = FockSpace(source) # Added for completeness.

"""
    quantize(index::Integer, identifier::Symbol, point::Point, flavor::Integer)::Mode

Quantizing a mode from a given `Point`.

### Input
- `identifier` The identifying atttibute key which the `point` object will be linked to.
- `point` The `Point` as the physical attribute or object to be quantized.
- `flavor` The flavor index of the `Mode`, don't mistaken this with the flavor count.

### Output
The quantized `Mode` object.
"""
function quantize(identifier::Symbol, point::Point, flavor::Integer)::Mode
    @assert(identifier == :offset || identifier == :pos)
    home::Point = origin(spaceof(point))
    # Since there are only one of the attribute :offset or :pos will take the point, the left over shall take the origin.
    couple::Pair{Symbol, Point} = identifier == :offset ? :pos => home : :offset => home
    # The new mode will take a group of q:$(name).
    return Mode([identifier => point, :flavor => flavor, couple])
end

"""
    quantize(identifier::Symbol, subset::Subset{Position}, count::Integer)::Subset{Mode}

Quantizing a set of mode from a given set of `Point`.

### Input
- `identifier` The identifying atttibute key which the `point` object will be linked to.
- `subset` The set of `Point` provided as the physical attributes or objects to be quantized.
- `count` The flavor count of the quantization, if it is greater than `1`, it means the given site defined by a `Point` in `subset` has more
  than one fermionic degree of freedom.

### Output
The quantized set of `Mode` objects.
"""
quantize(identifier::Symbol, subset::Subset{Position}, count::Int64)::Subset{Mode} = (
    Subset(quantize(identifier, point, flavor) for point in subset for flavor in 1:count))

"""
    FockMap(outspace::FockSpace, inspace::FockSpace, rep::SparseMatrixCSC{ComplexF64, Int64})
    FockMap(outspace::FockSpace, inspace::FockSpace, rep::AbstractArray{<:Number})
    FockMap(outspace::FockSpace, inspace::FockSpace, mapping::Dict{Tuple{Mode, Mode}, ComplexF64})
    FockMap(fockmap::FockMap; outspace::FockSpace = fockmap.outspace, inspace::FockSpace = fockmap.inspace)

Represents an mapping between two fockspaces with the same span of the underlying Hilbert space.

### Input
- `outspace` The output `FockSpace` of this map. If this object is the multiplier, then this will be the `outspace` of the resulting `FockMap`;
             if this object is the factor, then this must have the same span as the `inspace` of the multiplier.
- `inspace`  The input `FockSpace` of this map. If this object is the multiplier, then must have the same span as the `outspace` of the multiplier;
             if this object is the factor, then this will be the `inspace` of the resulting `FockMap`.
- `rep`      A complex sparse matrix represents the 2-point maps between the elements of the `inspace` & `outspace`.
- `mapping`  The values of the map have to be specified for a distinct 2-point pair, keyed by the pair in `Tuple`.

### Examples
- `FockMap(outspace::FockSpace, inspace::FockSpace, rep::SparseMatrixCSC{ComplexF64, Int64})` is used when every ingredients are precomputed.
- `FockMap(outspace::FockSpace, inspace::FockSpace, rep::AbstractArray{<:Number})` is used when the `rep` is an arbitary array like object.
- `FockMap(outspace::FockSpace, inspace::FockSpace, mapping::Dict{Tuple{Mode, Mode}, ComplexF64})` is used when the values of the map have to be specified
  for a distinct 2-point pair.
- `FockMap(fockmap::FockMap; outspace::FockSpace = fockmap.outspace, inspace::FockSpace = fockmap.inspace)` is for using new `outspace` and `inspace`.
"""
struct FockMap <: Element{SparseMatrixCSC{ComplexF64, Int64}}
    outspace::FockSpace
    inspace::FockSpace
    rep::SparseMatrixCSC{ComplexF64, Int64}

    FockMap(outspace::FockSpace{<: Any}, inspace::FockSpace{<: Any}, rep::SparseMatrixCSC{ComplexF64, Int64}) = new(outspace, inspace, rep)
    FockMap(outspace::FockSpace{<: Any}, inspace::FockSpace{<: Any}, rep::AbstractArray{<:Number}) = new(outspace, inspace, SparseMatrixCSC{ComplexF64, Int64}(rep))

    function FockMap(outspace::FockSpace{<: Any}, inspace::FockSpace{<: Any}, mapping::Dict{Tuple{Mode, Mode}, ComplexF64})::FockMap
        rep::SparseMatrixCSC{ComplexF64, Int64} = spzeros(dimension(outspace), dimension(inspace))
        for ((out_mode::Mode, in_mode::Mode), value::ComplexF64) in mapping
            rep[outspace.ordering[out_mode], inspace.ordering[in_mode]] = value
        end
        return new(outspace, inspace, rep)
    end

    FockMap(fockmap::FockMap; outspace::FockSpace{<: Any} = fockmap.outspace, inspace::FockSpace{<: Any} = fockmap.inspace, performpermute::Bool = true) = (
        FockMap(outspace, inspace, performpermute ? permute(fockmap, outspace=outspace, inspace=inspace) |> rep : fockmap |> rep))
end

Base.:show(io::IO, fockmap::FockMap) = print(io, string("$(typeof(fockmap))(in=$(fockmap.inspace), out=$(fockmap.outspace))"))

Base.:convert(::Type{SparseMatrixCSC{ComplexF64, Int64}}, source::FockMap) = source.rep

"""
    idmap(outspace::FockSpace, inspace::FockSpace)::FockMap

Create an injective `FockMap` from `inspace` to `outspace` of the same dimension, and the mapping pairs are determined by the order of both spaces,
i.e. the n-th element of the `inspace` will be mapped to the n-th element of the outspace.
"""
function idmap(outspace::FockSpace, inspace::FockSpace)::FockMap
    @assert(dimension(outspace) == dimension(inspace))
    FockMap(outspace, inspace, SparseMatrixCSC(Matrix{Float64}(I(dimension(outspace)))))
end

"""
    onesmap(outspace::FockSpace, inspace::FockSpace)::FockMap

Generate a `FockMap` full of `1`s from `inspace` to `outspace`.
"""
onesmap(outspace::FockSpace, inspace::FockSpace)::FockMap = FockMap(outspace, inspace, spzeros(dimension(outspace), dimension(inspace)) .+ 1)

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

"""
    columns(fockmap::FockMap, restrictspace::FockSpace)::FockMap

Restrict the `inspace` of the `fockmap` by a sub-fockspace `restrictspace`.
"""
function columns(fockmap::FockMap, restrictspace::FockSpace)::FockMap
    restrictindices::Vector{Integer} = [order(fockmap.inspace, mode) for mode in orderedmodes(restrictspace)]
    spmat::SparseMatrixCSC{ComplexF64, Int64} = sparse(view(rep(fockmap), :, restrictindices))
    return FockMap(fockmap.outspace, restrictspace, spmat)
end

"""
    rows(fockmap::FockMap, restrictspace::FockSpace)::FockMap

Restrict the `outspace` of the `fockmap` by a sub-fockspace `restrictspace`.
"""
function rows(fockmap::FockMap, restrictspace::FockSpace)::FockMap
    restrictindices::Vector{Integer} = [order(fockmap.outspace, mode) for mode in orderedmodes(restrictspace)]
    spmat::SparseMatrixCSC{ComplexF64, Int64} = sparse(view(rep(fockmap), restrictindices, :))
    return FockMap(restrictspace, fockmap.inspace, spmat)
end

"""
    rowsubmaps(fockmap::FockMap)::Base.Generator

Partition the `fockmap` into smaller `fockmaps` by the subspaces of its `outspace`.
"""
rowsubmaps(fockmap::FockMap)::Base.Generator = (rows(fockmap, subspace) for subspace in fockmap.outspace |> subspaces)

"""
    colsubmaps(fockmap::FockMap)::Base.Generator

Partition the `fockmap` into smaller `fockmaps` by the subspaces of its `inspace`.
"""
colsubmaps(fockmap::FockMap)::Base.Generator = (columns(fockmap, subspace) for subspace in fockmap.inspace |> subspaces)

"""
    restrict(fockmap::FockMap, outspace::FockSpace, inspace::FockSpace)::FockMap

Restrict the `outspace` & `inspace` of the `fockmap` by a sub-fockspaces `outspace` & `inspace` respectively.
"""
function restrict(fockmap::FockMap, outspace::FockSpace, inspace::FockSpace)::FockMap
    outindices::Vector{Integer} = [order(fockmap.outspace, mode) for mode in orderedmodes(outspace)]
    inindices::Vector{Integer} = [order(fockmap.inspace, mode) for mode in orderedmodes(inspace)]
    spmat::SparseMatrixCSC{ComplexF64, Int64} = sparse(view(rep(fockmap), outindices, inindices))
    return FockMap(outspace, inspace, spmat)
end

"""
    permute(source::FockMap; outspace::FockSpace=source.outspace, inspace::FockSpace=source.inspace)::FockMap

Permute the columns and rows of the representation of the `source` `FockMap` by `outspace` & `inspace` respectively.

### Input
- `source`   The target `FockMap` to be permuted.
- `outspace` A `FockSpace` with the same span as the `outspace` of the `source` `FockMap`.
- `inspace`  A `FockSpace` with the same span as the `inspace` of the `source` `FockMap`.
"""
function permute(source::FockMap; outspace::FockSpace=source.outspace, inspace::FockSpace=source.inspace)::FockMap
    row_rule::Vector{Int64} = orderingrule(source.outspace, outspace)
    col_rule::Vector{Int64} = orderingrule(source.inspace, inspace)
    return FockMap(outspace, inspace, SparseArrays.permute(rep(source), row_rule, col_rule))
end

Base.:-(target::FockMap)::FockMap = FockMap(target.outspace, target.inspace, -rep(target))

Base.:+(a::FockMap, b::FockMap)::FockMap = fockadd(a, b)

Base.:-(a::FockMap, b::FockMap)::FockMap = a + (-b)

function Base.:*(a::FockMap, b::FockMap)::FockMap
    @assert(hassamespan(a.inspace, b.outspace)) # Even if the fockspaces are different, composition works as long as they have same span.
    return FockMap(a.outspace, b.inspace, rep(a) * rep(permute(b, outspace=a.inspace, inspace=b.inspace)))
end

Base.:*(fockmap::FockMap, number::Number)::FockMap = FockMap(fockmap.outspace, fockmap.inspace, rep(fockmap) * number)
Base.:*(number::Number, fockmap::FockMap)::FockMap = fockmap * number

Base.:/(fockmap::FockMap, number::Number)::FockMap = FockMap(fockmap.outspace, fockmap.inspace, rep(fockmap) / number)

Base.:transpose(source::FockMap)::FockMap = FockMap(source.inspace, source.outspace, transpose(rep(source)))
Base.:adjoint(source::FockMap)::FockMap = FockMap(source.inspace, source.outspace, rep(source)')

"""
    eigmodes(fockmap::FockMap, attrs::Pair{Symbol}...)::Subset{Mode}

According to the provided Hermitian `FockMap`, generate a set of eigenmodes that corresponds one to one to the eigenvectors and eigenvalues when
performing eigenvalue decomposition.

### Input
- `fockmap` The source of the generation.
- `attrs`   Attributes to be inserted to the generated eigenmodes, please see `Mode` for more detains.
"""
function eigmodes(fockmap::FockMap, attrs::Pair{Symbol}...)::Subset{Mode}
    @assert(hassamespan(fockmap.inspace, fockmap.outspace))
    return Subset(
        Mode([:flavor => index, attrs...]) for index in 1:dimension(fockmap.inspace))
end

"""
    eigvecsh(fockmap::FockMap, attrs::Pair{Symbol}...)::FockMap

Perform eigenvalue decomposition to find the eigenvectors horizontally concatenated into a `FockMap`, ordered by the eigenvalues in ascending order.

### Input
- `fockmap` The source of the decomposition.
- `attrs`   Attributes to be inserted to the generated eigenmodes, please see `Mode` for more detains.

### Output
A `FockMap` with `outspace` of `fockmap` and the inspace of the associated eigenmodes, containing the horizontally concatenated eigenvectors.
"""
function eigvecsh(fockmap::FockMap, attrs::Pair{Symbol}...)::FockMap
    @assert(fockmap.inspace == fockmap.outspace)
    evecs::Matrix{ComplexF64} = eigvecs(Hermitian(Matrix(rep(fockmap))))
    return FockMap(fockmap.outspace, FockSpace(eigmodes(fockmap, attrs...)), SparseMatrixCSC{ComplexF64, Int64}(evecs))
end

"""
    eigvalsh(fockmap::FockMap, attrs::Pair{Symbol}...)::Vector{Pair{Mode, Float64}}

Perform eigenvalue decomposition to find the eigenvalues as pairs of eigenmode to eigenvalue ordered by the eigenmode attribute `index`, ordered by the
eigenvalues in ascending order.

### Input
- `fockmap` The source of the decomposition.
- `attrs`   Attributes to be inserted to the generated eigenmodes, please see `Mode` for more detains.
"""
function eigvalsh(fockmap::FockMap, attrs::Pair{Symbol}...)::Vector{Pair{Mode, Float64}}
    @assert(fockmap.inspace == fockmap.outspace)
    evs::Vector{Number} = eigvals(Hermitian(Matrix(rep(fockmap))))
    return [first(tup) => last(tup) for tup in Iterators.zip(eigmodes(fockmap, attrs...), evs)]
end

function LinearAlgebra.:eigvals(fockmap::FockMap, attrs::Pair{Symbol}...)::Base.Generator
    eigenvalues::Vector{Number} = fockmap |> rep |> Matrix |> eigvals
    return (pair.first => pair.second for pair in Iterators.zip(eigmodes(fockmap, attrs...), eigenvalues))
end

function LinearAlgebra.:eigvecs(fockmap::FockMap, attrs::Pair{Symbol}...)::FockMap
    eigenvectors::Matrix = fockmap |> rep |> Matrix |> eigvecs
    return FockMap(fockmap.outspace, FockSpace(eigmodes(fockmap, attrs...)), eigenvectors)
end

function LinearAlgebra.:eigen(fockmap::FockMap, attrs::Pair{Symbol}...)::Tuple{Base.Generator, FockMap}
    eigenvalues, eigenvectors = fockmap |> rep |> Matrix |> eigen
    modes::Subset{Mode} = eigmodes(fockmap, attrs...)
    return (
        (pair |> first => pair |> last for pair in Iterators.zip(eigmodes(fockmap, attrs...), eigenvalues)),
        FockMap(fockmap.outspace, FockSpace(eigmodes(fockmap, attrs...)), eigenvectors))
end

"""
    eigh(fockmap::FockMap, attrs::Pair{Symbol}...)::Tuple{Vector{Pair{Mode, Float64}}, FockMap}

Perform eigenvalue decomposition to find the eigenvalues and eigenvectors simultaneously.

### Input
- `fockmap` The source of the decomposition.
- `attrs`   Attributes to be inserted to the generated eigenmodes, please see `Mode` for more detains.

### Output
The first returned value will be pairs of eigenmode to eigenvalue ordered by the eigenmode attribute `index`; the second returned value is a `FockMap`
with `outspace` of `fockmap` and the inspace of the associated eigenmodes, containing the horizontally concatenated eigenvectors. The eigenmodes are
ordered by the eigenvalues in ascending order.
"""
function eigh(fockmap::FockMap, attrs::Pair{Symbol}...)::Tuple{Vector{Pair{Mode, Float64}}, FockMap}
    decomposed::Eigen = eigen(Hermitian(Matrix(rep(fockmap))))
    modes::Subset{Mode} = eigmodes(fockmap, attrs...)
    evals::Vector{Pair{Mode, Float64}} = [m => v for (m, v) in Iterators.zip(modes, decomposed.values)]
    evecs::FockMap = FockMap(fockmap.outspace, FockSpace(modes), decomposed.vectors)
    return evals, evecs
end

"""
    fourier(momentums::Subset{Momentum}, inspace::FockSpace)::FockMap

Create a `FockMap` corresponds to a Fourier transform of a `FockSpace`. This map assumed orthogonality of modes with the attribute `:offset` dropped, which means
entries corresponds to different fermionic site within the same translational invariant unit cell will be default to `0 + 0im`.

### Input
- `momentums` The momentums which spans the output reciprocal subspace, for `spaceof(momentum) isa MomentumSpace`.
- `inspace` The input space, all contituent modes must have the attribute `:offset` defined or result in errors. The `inspace` should include all possible basis modes
  so that they can be identified and span the momentum `FockSpace`.

### Output
The `FockMap` represents this specific Fourier transform, with sizes `(N, M)` which `N = length(momentums) * count` for `count` is the number of fermionic site
within the translational invariant unit cell, supplied by `inmodes`; `M = length(inmodes)`. The `inspace` of the returned `FockMap` will equates to `FockSpace(inmodes)`;
the `outspace` is the product of `momentums` and the supplied fermionic sites.
"""
function fourier(momentums::Subset{Momentum}, inspace::FockSpace)::FockMap
    ∑𝑘::Matrix{Float64} = hcat([𝑘 |> euclidean |> pos for 𝑘 in momentums]...)
    inmodes::Subset{Mode} = orderedmodes(inspace)
    basismodes::Subset{Mode} = removeattr(inmodes, :offset)
    outspace::FockSpace = sparsefock(basismodes, momentums)
    return fourier(outspace, inspace, ∑𝑘, basismodes)
end

"""
    fourier(outspace::FockSpace, inspace::FockSpace)::FockMap

Create a `FockMap` corresponds to a Fourier transform of a physical fockspace `inspace` to Fourier fockspace `outspace`. This map assumed orthogonality of modes with the
attribute `:offset` dropped, which means entries corresponds to different fermionic site within the same translational invariant unit cell will be default to `0 + 0im`.

### Input
- `outspace` The output space of the Fourier transform, which is the momentum `FockSpace`, most likely with type `FockSpace{Crystal}`, noted that the basis modes in
  each momentum subspace should matches with the possible basis modes within `inspace`.
- `inspace`  The input space, all contituent modes must have the attribute `:offset` defined or result in errors.

### Output
The `FockMap` represents this specific Fourier transform, with sizes `(N, M)` which `N = dimension(outspace); M = dimension(inspace)`.
"""
function fourier(outspace::FockSpace, inspace::FockSpace)::FockMap
    ∑𝑘::Matrix{Float64} = hcat([getattr(first(partition), :offset) |> euclidean |> pos for partition in rep(outspace)]...)
    basismodes::Subset{Mode} = removeattr(outspace |> rep |> first, :offset) # Assumed the similarity in structure for each partitions.
    return fourier(outspace, inspace, ∑𝑘, basismodes)
end

""" Internal use only. """
function fourier(outspace::FockSpace, inspace::FockSpace, momentummatrix::Matrix{Float64}, basismodes::Subset{Mode})::FockMap
    values::Array{ComplexF64} = zeros(ComplexF64, length(basismodes), size(momentummatrix, 2), dimension(inspace))
    for ((n, basismode), (m, inmode)) in Iterators.product(enumerate(basismodes), enumerate(orderedmodes(inspace)))
        if removeattr(inmode, :offset) != basismode continue end
        offset::Point = inmode |> getattr(:offset) |> euclidean
        values[n, :, m] = exp.(-1im * momentummatrix' * pos(offset))
    end
    spmat::SparseMatrixCSC = SparseMatrixCSC(reshape(values, (length(basismodes) * size(momentummatrix, 2), dimension(inspace))))
    return FockMap(outspace, inspace, spmat)
end

"""
Perform addition of two `FockMap`s, if they carries `inspace` and `outspace` of same span, then this is just a sum of representation values; if they carries
non-overlapping `inspace` and `outspace`, this corresponds to a direct sum; if they have overlapping but different span of `inspace` or `outspace`, the result will
be a `FockMap` with `outspace` and `inspace` with 3 subspaces of (!intersect x 2, intersect x 1), with the corresponding representation summed. Noted that this
function is accessed via `a::FockMap + b::FockMap`.
"""
function fockadd(a::FockMap, b::FockMap)::FockMap
    outspacesamespan::Bool = hassamespan(a.outspace, b.outspace)
    inspacesamespan::Bool = hassamespan(a.inspace, b.inspace)

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

    outspace2::FockSpace = intersect(a.outspace, b.outspace)
    inspace2::FockSpace = intersect(a.inspace, b.inspace)

    if outspace2 |> dimension == 0 && inspace2 |> dimension == 0
        return directsum(a, b)
    end

    outspace1::FockSpace = a.outspace - b.outspace
    outspace3::FockSpace = b.outspace - a.outspace

    inspace1::FockSpace = a.inspace - b.inspace
    inspace3::FockSpace = b.inspace - a.inspace

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
        data[order(outspace, outmodes |> first):order(outspace, outmodes |> last), order(inspace, inmodes |> first):order(inspace, inmodes |> last)] += (restricted |> rep)
    end

    internaladdition(a, outspace1, inspace1)
    internaladdition(a, outspace1, inspace2)
    internaladdition(a, outspace2, inspace1)
    internaladdition(a, outspace2, inspace2)
    internaladdition(b, outspace2, inspace2)
    internaladdition(b, outspace2, inspace3)
    internaladdition(b, outspace3, inspace2)
    internaladdition(b, outspace3, inspace3)

    return FockMap(outspace, inspace, data)
end

function fockaddsamespan(a::FockMap, b::FockMap)::FockMap
    data::SparseMatrixCSC{ComplexF64, Int64} = (
        (a |> rep) + (Quantum.permute(b, outspace=a.outspace, inspace=a.inspace) |> rep))
    return FockMap(a.outspace, a.inspace, data)
end

function directsum(a::FockMap, b::FockMap)::FockMap
    outspace::FockSpace = a.outspace + b.outspace
    inspace::FockSpace = a.inspace + b.inspace
    data::SparseMatrixCSC{ComplexF64, Int64} = spzeros(outspace |> dimension, inspace |> dimension)
    data[1:(a.outspace |> dimension), 1:(a.inspace |> dimension)] += a |> rep
    data[(a.outspace |> dimension) + 1:end, (a.inspace |> dimension) + 1:end] += b |> rep
    return FockMap(outspace, inspace, data)
end

function directsum(fockmaps)::FockMap
    outspace::FockSpace = Iterators.map(fockmap -> fockmap.outspace, fockmaps) |> fockspaceunion
    inspace::FockSpace = Iterators.map(fockmap -> fockmap.inspace, fockmaps) |> fockspaceunion
    data::SparseMatrixCSC{ComplexF64, Int64} = spzeros(outspace |> dimension, inspace |> dimension)
    function filldata(fockmap::FockMap)
        # The procedure beneath assumes that the fockspace elements are concatenated in order during union operations.
        outmodes::Subset{Mode} = fockmap.outspace |> orderedmodes
        outrange::UnitRange = order(outspace, outmodes |> first):order(outspace, outmodes |> last)
        inmodes::Subset{Mode} = fockmap.inspace |> orderedmodes
        inrange::UnitRange = order(inspace, inmodes |> first):order(inspace, inmodes |> last)
        data[outrange, inrange] += fockmap |> rep
    end
    foreach(filldata, fockmaps)
    return FockMap(outspace, inspace, data)
end

function crystalsubmaps(fockmap::FockMap)
    if !(fockmap.inspace isa FockSpace{Crystal} && fockmap.outspace isa FockSpace{Crystal})
        @error("The in/out spaces of the fock map must be crystal fock-spaces!")
    end
    if !hassamespan(fockmap.inspace, fockmap.outspace)
        @error("Not a Hermitian!")
    end
    return (k => restrict(fockmap, fockspace, fockspace) for (k, fockspace) in fockmap.inspace |> crystalsubspaces)
end
export crystalsubmaps

struct CrystalSpectrum
    crystal::Crystal
    eigenmodes::Dict{Momentum, Subset{Mode}}
    eigenvalues::Dict{Mode, Number}
    eigenvectors::Dict{Momentum, FockMap}
end
export CrystalSpectrum

Base.:show(io::IO, spectrum::CrystalSpectrum) = print(io, string("$(spectrum |> typeof)(entries=$(spectrum.eigenvalues |> length))"))

function crystalspectrum(fockmap::FockMap)::CrystalSpectrum
    blocks = fockmap |> crystalsubmaps
    crystalemodes::Dict{Momentum, Subset{Mode}} = Dict()
    crystalevals::Dict{Mode, Number} = Dict()
    crystalevecs::Dict{Momentum, FockMap} = Dict()
    for (k, submap) in blocks
        eigenvalues, eigenvectors = eigh(submap, :offset => k)
        crystalemodes[k] = Subset(m for (m, _) in eigenvalues)
        crystalevecs[k] = eigenvectors
        for (m, v) in eigenvalues
            crystalevals[m] = v
        end
    end
    return CrystalSpectrum(fockmap.inspace |> crystalof, crystalemodes, crystalevals, crystalevecs)
end
export crystalspectrum

struct EigenSpectrum
    eigenvalues::Dict{Mode, Number}
    eigenvectors::FockMap
end
export EigenSpectrum

Base.:show(io::IO, spectrum::EigenSpectrum) = print(io, string("$(spectrum |> typeof)(entries=$(spectrum.eigenvalues |> length))"))

function eigenspectrum(hermitian::FockMap, attrs::Pair{Symbol}...)::EigenSpectrum
    evals, eigenvectors = eigh(hermitian, attrs...)
    eigenvalues::Dict{Mode, Number} = Dict(evals)
    return EigenSpectrum(eigenvalues, eigenvectors)
end
export eigenspectrum

function groupbyeigenvalues(spectrum::EigenSpectrum; groupingthreshold::Number = 1e-2)
    denominator::Integer = (1 / groupingthreshold) |> round |> Integer
    actualvalues::Dict{Rational, Number} = Dict(hashablereal(v, denominator) => v for (_, v) in spectrum.eigenvalues)
    items::Base.Generator = (hashablereal(v, denominator) => m for (m, v) in spectrum.eigenvalues)
    groups::Dict{Rational, Vector{Mode}} = foldl(items; init=Dict{Rational, Vector{Mode}}()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end
    sortedrationals::Vector{Rational} = sort([(groups |> keys)...])
    return (actualvalues[r] => groups[r] |> Subset for r in sortedrationals)
end
export groupbyeigenvalues

function LinearAlgebra.log(fockmap::FockMap)::FockMap
    mat::SparseMatrixCSC = fockmap |> rep |> Matrix |> log |> SparseMatrixCSC
    return FockMap(fockmap.outspace, fockmap.inspace, mat)
end

Base.iszero(fockmap::FockMap)::Bool = iszero(fockmap |> rep)

function LinearAlgebra.:svd(fockmap::FockMap)::Tuple{FockMap, Base.Generator, FockMap}
    leftmodes::Subset{Mode} = Subset(Mode([:svdindex => n]) for n in 1:dimension(fockmap.inspace))
    rightmodes::Subset{Mode} = Subset(Mode([:svdindex => n]) for n in 1:dimension(fockmap.outspace))
    U, Σ, Vt = fockmap |> rep |> Matrix |> svd
    svdvalues::Base.Generator = (
        (Mode([:svdindex => n]), Mode([:svdindex => n])) => Σ[n] for n in 1:min(leftmodes |> length, rightmodes |> length))
    return FockMap(fockmap.outspace, leftmodes |> FockSpace, U), svdvalues, FockMap(rightmodes |> FockSpace, fockmap.inspace, Vt')
end

"""
    columnspec(fockmap::FockMap)::Vector{Pair{Mode, ComplexF64}}

Extract the column values as a `Mode` to `ComplexF64` pair from a `N×1` `FockMap`, this is used when you need to visualize the column spectrum.
"""
function columnspec(fockmap::FockMap)::Vector{Pair{Mode, ComplexF64}}
    @assert(dimension(fockmap.inspace) == 1)
    mat::SparseMatrixCSC{ComplexF64, Int64} = rep(fockmap)
    return [outmode => mat[order(fockmap.outspace, outmode), 1] for outmode in orderedmodes(fockmap.outspace)]
end

end
