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
    Mode(attrs::Pair{Symbol, <:Any}...) = Mode(Dict(attrs...))
end
export Mode

""" Stringify the `Mode` object. """
Base.:string(mode::Mode)::String = join((k => v for (k, v) in mode |> getattrs), " ")

""" Display the number of attributes that identifies this `Mode`. """
Base.:show(io::IO, mode::Mode) = print(io, string("$(typeof(mode))$(tuple(keys(mode.attrs)...))"))

""" Allow for offset the `:offset` attribute of a mode. """
Base.:+(mode::Mode, offset::Point)::Mode = setattr(mode, :offset => getattr(mode, :offset) + offset)
Base.:-(mode::Mode, offset::Point)::Mode = mode + (-offset)

"""
    getspace(mode::Mode)

The space of a `Mode` comes from the physical quantities its defined on, such as `:offset` and `:pos`, if none of those are defined,
it will be `euclidean(RealSpace, 1)` as the mode and it's siblings can always be parameterized by a scalar.

### Output
The space of the attribute `:offset` if available, fall back to `:pos` if otherwise; returns a Euclidean space of dimension `1` if both `:offset` & `:pos` is missing.
"""
function Zipper.:getspace(mode::Mode)
    # :offset have a higher priority in determining the space of the mode.
    if hasattr(mode, :offset) return getspace(getattr(mode, :offset)) end
    if hasattr(mode, :pos) return getspace(getattr(mode, :pos)) end
    # If the mode does not based on any physical position or quantities for associating with any space, then it will simply lives
    # in a 1D euclidean space as the mode and it's siblings can always be parameterized by a scalar.
    return euclidean(RealSpace, 1)
end

""" Shorthand to create a one element `Subset{Mode}`. """
Zipper.:Subset(modes::Mode...) = Subset(m for m in modes)

"""
    getpos(mode::Mode)::Point

Get the actual position of the mode, this method only works when `:offset` and `:pos` are defined in the same space.
"""
Zipper.:getpos(mode::Mode)::Point = convert(Point, mode)

"""
    hasattr(mode::Mode, key::Symbol)::Bool

Check if the `mode` has the attribute identified by `key`.
"""
hasattr(mode::Mode, key::Symbol)::Bool = haskey(mode.attrs, key)
export hasattr

"""
    getattr(mode::Mode, key::Symbol)

Retrieve the attribute value identified by `key` from `mode`.
"""
getattr(mode::Mode, key::Symbol) = mode.attrs[key]
export getattr

"""
    getattr(key::Symbol)

Shorthand of `getattr(v, key::Symbol)` with the pipe operator `|>`.

### Examples
The line `mode |> getattr(:flavor)` is equal to `getattr(mode, :flavor)`.
"""
getattr(key::Symbol) = v -> getattr(v, key)

""" A shorthand for getting the mode attributes for visualization in the REPL. """
getattrs(mode::Mode) = mode.attrs
export getattrs

"""
    removeattr(mode::Mode, keys::Symbol...)::Mode

Create a **copy** of `mode` **without** the attributes identified by `keys`.

### Examples
- To remove the attribute of `:offset` and `:pos`, we use `removeattr(mode, :offset, :pos)`.
"""
removeattr(mode::Mode, keys::Symbol...)::Mode = Mode(Dict(filter(p -> !(p.first ∈ keys), mode.attrs)))
export removeattr

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
export setattr

"""
    setattr(attrs::Pair{Symbol}...)

Shorthand of `setattr(v, attrs::Pair{Symbol}...)::Mode` with the pipe operator `|>`.

### Examples
The line `mode |> setattr(:flavor => 1)` is equal to `setattr(mode, :flavor => 1)`.
"""
setattr(attrs::Pair{Symbol}...) = v -> setattr(v, attrs...)

""" A shorthand for setting an attribute to a mode. """
Base.:&(mode::Mode, attr::Pair{Symbol}) = setattr(mode, attr)

"""
    setattr(subset::Subset{Mode}, attrs::Pair{Symbol}...)::Subset{Mode}

Create a **copy** of `modes` with the new attributes identified by `keys` added to the current attributes for each mode in `modes`, if the attribute
exists in the modes, the current record will be overwritten.

### Examples
- To add `:offset` and `:flavor`, we use `newmodes = setattr(modes, :offset => point, :flavor => 1)`
"""
setattr(subset::Subset{Mode}, attrs::Pair{Symbol}...)::Subset{Mode} = Subset(setattr(mode, attrs...) for mode in subset)

"""
    commonattrs(modes)::Base.Generator

Retrieve all the attributes within the group of modes that has the same values and returned as a generator of `Symbol`.
"""
function commonattrs(modes)::Base.Generator
    attrs::Base.Generator = (Set(k => v for (k, v) in mode |> getattrs) for mode in modes)
    return (k for (k, _) in intersect(attrs...))
end
export commonattrs

"""
    indexmodes(modes)::Subset{Mode}

Retrieve the representative modes of this group of modes with their corresponding common attributes stripped.
"""
function indexmodes(modes)::Subset{Mode}
    cleanattrs::Base.Generator = modes |> commonattrs
    return Subset(mode for mode in modes) |> removeattr(cleanattrs...) 
end
export indexmodes

getorbital(mode::Mode, default::BasisFunction)::BasisFunction = hasattr(mode, :orbital) ? getattr(mode, :orbital) : default
export getorbital

getorbital(default::BasisFunction = swave) = mode -> getorbital(mode, default)

setorbital(mode::Mode, basis::BasisFunction)::Mode = setattr(mode, :orbital => basis)
export setorbital

setorbital(basis::BasisFunction) = mode -> setorbital(mode, basis)

"""
    spanoffset(basismodes::Subset{Mode}, points::Subset{<: Point})::Subset{Mode}

Given a set of `basismodes`, and the generator `points`, span the basis modes to the generator `points` with attribute `:offset`, the primary ordering will be
the ordering of `points`, then follows the ordering of `basismodes`.
"""
spanoffset(basismodes::Subset{Mode}, points::Subset{<: Point})::Subset{Mode} = Subset(setattr(mode, :offset => point) for point in points for mode in basismodes)
export spanoffset

""" By this conversion, one can obtain the actual position of the mode, this method only works when `:offset` and `:pos` are defined in the same space. """
Base.:convert(::Type{Point}, source::Mode)::Point = getattr(source, :offset) + getattr(source, :pos)

""" Two `Mode` objects are equivalent if they held the same informations. """
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
    
    FockSpace(input::Vector{Mode}; reflected=Nothing) = FockSpace(Subset(input), reflected=reflected)
    FockSpace(input::Base.Generator; reflected=Nothing) = FockSpace(Subset(input), reflected=reflected)
    FockSpace(input::Base.Iterators.Flatten; reflected=Nothing) = FockSpace(Subset(input), reflected=reflected)
end
export FockSpace

""" To retrieve the reflected property this `FockSpace` is reflected to. """
getreflected(fockspace::FockSpace) = fockspace.reflected
export getreflected

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
export sparsefock

"""
    crystalfock(basismodes::Subset{Mode}, crystal::Crystal)::FockSpace

A short hand to build the crystal fockspace, which is the fockspace containing all modes spanned from `basismodes` by the brillouin zone of the `crystal`.
"""
crystalfock(basismodes::Subset{Mode}, crystal::Crystal)::CrystalFock = FockSpace(sparsefock(basismodes, brillouinzone(crystal)), reflected=crystal)
export crystalfock

""" To create a regional `FockSpace`."""
function FockSpace{Region}(input)
    region::Region = Subset(m |> getpos for m in input)
    return FockSpace(input, reflected=region)
end

""" Get the reflected region of the regional `FockSpace`. """
Zipper.:getregion(regionfock::FockSpace{Region})::Region = regionfock |> getreflected

""" Shorthand alias for `FockSpace{Crystal}`. """
CrystalFock = FockSpace{Crystal}
export CrystalFock

""" Displays the fock type, subspace count and dimension information of a `FockSpace`. """
Base.:show(io::IO, fockspace::FockSpace) = print(io, string("$(typeof(fockspace))(sub=$(fockspace |> subspacecount), dim=$(fockspace |> dimension))"))
Base.:show(io::IO, ::Type{FockSpace{Region}}) = print(io, string("RegionFock"))
Base.:show(io::IO, ::Type{CrystalFock}) = print(io, string("CrystalFock"))

""" Check whether a `Mode` is in the `FockSpace`. """
Base.:in(mode::Mode, fockspace::FockSpace)::Bool = haskey(fockspace.ordering, mode)

""" Allow the retrieval of `Mode` at a given `order` index with syntax `fockspace[order]`. """
Base.:getindex(fockspace::FockSpace, order::Integer)::Mode = (fockspace |> orderedmodes)[order]
Base.:getindex(fockspace::FockSpace, range::UnitRange)::FockSpace = (fockspace |> orderedmodes)[range] |> FockSpace
Base.:getindex(fockspace::FockSpace, mode::Mode)::Mode = mode ∈ fockspace ? mode : error("The mode is not in the fockspace!")

"""
    fockspaceunion(fockspaces)::FockSpace

Union of fockspaces, the resulting fockspace will have the same span as the union of the underlying modes of the fockspaces.

### Input
- `fockspaces` An iterable of fockspaces to be unioned.
"""
function fockspaceunion(fockspaces)::FockSpace
    subspaces::Subset{Subset{Mode}} = subsetunion(fockspace |> rep for fockspace in fockspaces if fockspace |> dimension != 0)
    modes::Subset{Mode} = subspaces |> Iterators.flatten
    ordering::Dict{Mode, Integer} = Dict(mode => index for (index, mode) in enumerate(modes))
    return FockSpace(subspaces::Subset{Subset{Mode}}, ordering)
end
export fockspaceunion

Base.:union(fockspaces::FockSpace...)::FockSpace = fockspaceunion(fockspaces)

Base.:+(a::FockSpace, b::FockSpace)::FockSpace = union(a, b)

Base.:-(a::FockSpace, b::FockSpace)::FockSpace = FockSpace(setdiff(a |> orderedmodes , b |> orderedmodes))

Base.:intersect(a::FockSpace, b::FockSpace) = FockSpace(intersect(a |> orderedmodes, b |> orderedmodes))

Base.:iterate(fockspace::FockSpace, i...) = iterate(fockspace |> orderedmodes, i...)
Base.:length(fockspace::FockSpace) = fockspace |> dimension

"""
    dimension(fockspace::FockSpace)

Returns the number of unique member modes within the `fockspace`, each of those represents a vector from the Hilbert space.
"""
Zipper.:dimension(fockspace::FockSpace) = length(fockspace.ordering) # This is a short cut to retrieve the length, which is the dimension.

"""
    getcrystal(crystalfock::CrystalFock)::Crystal

Shorthand for retrieving the `Crystal` of a `CrystalFock`.
"""
Zipper.:getcrystal(crystalfock::CrystalFock)::Crystal = crystalfock |> getreflected

"""
    crystalsubsets(crystalfock::FockSpace{Crystal})::Dict{Momentum, Subset{Mode}}

    Retrieve mappings from the crystal momentums to the corresponding `Subset{Mode}`.
"""
crystalsubsets(crystalfock::FockSpace{Crystal})::Dict{Momentum, Subset{Mode}} = Dict(commonattr(subspace, :offset) => subspace for subspace in crystalfock |> rep)
export crystalsubsets

"""
    crystalsubspaces(crystalfock::FockSpace{Crystal})::Base.Generator

Retrieve mappings from the crystal momentums to the corresponding fockspaces.
"""
function crystalsubspaces(crystalfock::FockSpace{Crystal})::Base.Generator
    function subsettofockspace(subset::Subset{Mode})::Pair{Momentum, FockSpace}
        k::Momentum = commonattr(subset, :offset)
        return k => FockSpace(subset, reflected=k)
    end
    return (subspace |> subsettofockspace for subspace in crystalfock |> rep)
end
export crystalsubspaces

"""
    getmomentum(subspace::FockSpace{Momentum})::Momentum

Shorthand for retrieving the `Momentum` of a `FockSpace{Momentum}`.
"""
getmomentum(subspace::FockSpace{Momentum})::Momentum = subspace.reflected
export getmomentum

"""
    commonattr(modes, key::Symbol)

Retrieve the common attribute associated to `key` of all the child modes of the `fockspace`, and throws assertion error if the attribute
is not unique within the `fockspace`.
"""
function commonattr(modes, key::Symbol)
    set::Set = Set()
    foreach(m -> push!(set, getattr(m, key)), modes)
    @assert(length(set) == 1, "The modes in this fockspace does not share the same attr `$(key)`!")
    return first(set)
end
export commonattr

""" Shorthand for accessing the attribute information of all modes within the `FockSpace` for visualization. """
modeattrs(fockspace::FockSpace)::OrderedSet{Dict} = OrderedSet(mode |> getattrs for mode in fockspace |> orderedmodes)
modeattrs(subset::Subset{Mode})::OrderedSet{Dict} = OrderedSet(mode |> getattrs for mode in subset)
export modeattrs

"""
    unitcellfock(crystalfock::FockSpace{Crystal})::FockSpace

Retrieve the unit cell fockspace of the system from a `FockSpace{Crystal}`, positioned at the origin of the parent `AffineSpace`.
"""
function unitcellfock(crystalfock::FockSpace{Crystal})::FockSpace{Region}
    firstpartition::Subset{Mode} = crystalfock |> rep |> first
    originpoint::Point = firstpartition |> first |> getattr(:pos) |> getspace |> getorigin
    return FockSpace{Region}(firstpartition |> setattr(:offset => originpoint))
end
export unitcellfock

"""
    subspaces(fockspace::FockSpace)::Base.Generator

Retrieve the sub-fockspaces of the `fockspace`.
"""
subspaces(fockspace::FockSpace)::Base.Generator = (FockSpace(partition) for partition in rep(fockspace))
export subspaces

"""
    subspacecount(fockspace::FockSpace)::Integer

Get the number of sub-fockspaces of this `fockspace`.
"""
subspacecount(fockspace::FockSpace)::Integer = fockspace |> rep |> length
export subspacecount

"""
    flattensubspaces(fockspace::FockSpace)::FockSpace

Merge all subspaces within the `fockspace`.
"""
flattensubspaces(fockspace::FockSpace)::FockSpace = FockSpace(Subset(mode for mode in orderedmodes(fockspace)))
export flattensubspaces

"""
    fockorder(fockspace::FockSpace, mode::Mode)::Int64

Since the fockspace have implicit ordering, this function returns the order index of the `mode` in the `fockspace`.
"""
fockorder(fockspace::FockSpace, mode::Mode)::Int64 = fockspace.ordering[mode]
export fockorder

"""
    getmodes(fockspace::FockSpace)::Set{Mode}

Returns an unordered set of modes of `fockspace`, this is a more efficient way to retrieve the underlying modes for a fockspace with
more than one partitions.
"""
getmodes(fockspace::FockSpace)::Set{Mode} = Set(keys(fockspace.ordering)) # This is the most efficient way to get all distinct modes.
export getmodes

"""
    orderedmodes(fockspace::FockSpace)::Subset{Mode}

Returns an ordered `Subset` of modes of `fockspace`.
"""
orderedmodes(fockspace::FockSpace)::Subset{Mode} = flatten(rep(fockspace))
export orderedmodes

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
export orderingrule

"""
    hassamespan(a::FockSpace, b::FockSpace)::Bool

Check if fockspaces of `a` and `b` shares the same set of underlying modes regardless of partitions.
"""
Zipper.:hassamespan(a::FockSpace, b::FockSpace)::Bool = getmodes(a) == getmodes(b)

"""
    issparse(fockspace::FockSpace)::Bool

Check if `fockspace` is a sparse fockspace.
"""
issparse(fockspace::FockSpace)::Bool = subspacecount(fockspace) > 1
export issparse

"""
    sparsegrouping(fockspace::FockSpace, byattrs::Symbol...)::FockSpace

Grouping the modes within the `fockspace` by the attributes specified by `byattrs`, the resulting fockspace will have a sparse structure
with the subspaces containing modes with the same specified attributes.
"""
function sparsegrouping(fockspace::FockSpace, byattrs::Symbol...)::FockSpace
    getidentifier(mode)::Vector = [mode |> getattr(attr) for attr in byattrs]

    identifiedmodes::Base.Generator = (m |> getidentifier => m for m in fockspace)
    groups::Dict{Vector, Vector} = foldl(identifiedmodes; init=Dict{Vector, Vector}()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end

    return fockspaceunion(group |> FockSpace for group in groups |> values)
end
export sparsegrouping

""" Shorthand for `sparsegrouping` which generates a function, that takes a `FockSpace` and returns the sparse grouped fockspace. """
sparsegrouping(byattrs::Symbol...)::Function = fockspace -> sparsegrouping(fockspace, byattrs...)

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
    home::Point = getorigin(getspace(point))
    # Since there are only one of the attribute :offset or :pos will take the point, the left over shall take the origin.
    couple::Pair{Symbol, Point} = identifier == :offset ? :pos => home : :offset => home
    # The new mode will take a group of q:$(name).
    return Mode([identifier => point, :flavor => flavor, couple])
end
export quantize

"""
    quantize(identifier::Symbol, subset::Subset{Offset}, count::Integer)::Subset{Mode}

Quantizing a set of mode from a given set of `Point`.

### Input
- `identifier` The identifying atttibute key which the `point` object will be linked to.
- `subset` The set of `Point` provided as the physical attributes or objects to be quantized.
- `count` The flavor count of the quantization, if it is greater than `1`, it means the given site defined by a `Point` in `subset` has more
  than one fermionic degree of freedom.

### Output
The quantized set of `Mode` objects.
"""
quantize(identifier::Symbol, subset::Subset{Offset}, count::Int64)::Subset{Mode} = (
    Subset(quantize(identifier, point, flavor) for point in subset for flavor in 1:count))

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

    SparseFockMap(outspace::FockSpace{<: Any}, inspace::FockSpace{<: Any}, rep::SparseMatrixCSC{ComplexF64, Int64}) = new{outspace |> typeof, inspace |> typeof}(outspace, inspace, rep)
    SparseFockMap(outspace::FockSpace{<: Any}, inspace::FockSpace{<: Any}, rep::AbstractArray{<:Number}) = new{outspace |> typeof, inspace |> typeof}(outspace, inspace, SparseMatrixCSC{ComplexF64, Int64}(rep))

    function SparseFockMap(outspace::FockSpace{<: Any}, inspace::FockSpace{<: Any}, mapping::Dict{Tuple{Mode, Mode}, T})::SparseFockMap where {T <: Complex}
        rep::SparseMatrixCSC{ComplexF64, Int64} = spzeros(dimension(outspace), dimension(inspace))
        for ((out_mode::Mode, in_mode::Mode), value::ComplexF64) in mapping
            rep[outspace.ordering[out_mode], inspace.ordering[in_mode]] = value
        end
        return SparseFockMap(outspace, inspace, rep)
    end

    SparseFockMap(fockmap::FockMap; outspace::FockSpace{<: Any} = fockmap|>getoutspace, inspace::FockSpace{<: Any} = fockmap|>getinspace, performpermute::Bool = true) = (
        SparseFockMap(outspace, inspace, performpermute ? permute(fockmap, outspace=outspace, inspace=inspace) |> rep : fockmap |> rep))
end
export SparseFockMap

"""
Using `FockMap` as the entry point of instantiating a default `SparseFockMap` object.
"""
# Using SparseFockMap as the default implementation of FockMap
FockMap(outspace::FockSpace{<: Any}, inspace::FockSpace{<: Any}, rep::SparseMatrixCSC{ComplexF64, Int64}) = SparseFockMap(outspace, inspace, rep)
FockMap(outspace::FockSpace{<: Any}, inspace::FockSpace{<: Any}, rep::AbstractArray{<:Number}) = SparseFockMap(outspace, inspace, rep)
FockMap(outspace::FockSpace{<: Any}, inspace::FockSpace{<: Any}, mapping::Dict{Tuple{Mode, Mode}})::SparseFockMap= SparseFockMap(outspace, inspace, mapping)
FockMap(fockmap::FockMap; outspace::FockSpace{<: Any} = fockmap|>getoutspace, inspace::FockSpace{<: Any} = fockmap|>getinspace, performpermute::Bool = true) = SparseFockMap(fockmap, outspace=outspace, inspace=inspace, performpermute=performpermute)

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

""" Allows the retrieval of mapping data via `fockmap[outmode, inmode]`. """
function Base.:getindex(fockmap::FockMap, outmode::Mode, inmode::Mode)::Complex
    inspaceorder::Integer = fockorder(fockmap|>getinspace, inmode)
    outspaceorder::Integer = fockorder(fockmap|>getoutspace, outmode)
    return (fockmap |> rep)[CartesianIndex(outspaceorder, inspaceorder)]
end

LinearAlgebra.:tr(fockmap::FockMap) = fockmap|>rep|>tr

LinearAlgebra.:det(fockmap::FockMap) = fockmap|>rep|>det

# ===================================================================================================================================================
# Added to support getindex of FockMap objects.
Base.:getindex(fockmap::FockMap, row, col) = restrict(fockmap, (fockmap |> getoutspace)[row] |> FockSpace, (fockmap |> getinspace)[col] |> FockSpace)
Base.:getindex(fockmap::FockMap, ::Colon, col) = restrict(fockmap, fockmap |> getoutspace, (fockmap |> getinspace)[col] |> FockSpace)
Base.:getindex(fockmap::FockMap, row, ::Colon) = restrict(fockmap, (fockmap |> getoutspace)[row] |> FockSpace, fockmap |> getinspace)
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
    restrictindices::Vector{Integer} = [fockorder(fockmap|>getinspace, mode) for mode in orderedmodes(restrictspace)]
    spmat::SparseMatrixCSC{ComplexF64, Int64} = sparse(view(rep(fockmap), :, restrictindices))
    return FockMap(fockmap|>getoutspace, restrictspace, spmat)
end
export columns

"""
    rows(fockmap::FockMap, restrictspace::FockSpace)::FockMap

Restrict the `outspace` of the `fockmap` by a sub-fockspace `restrictspace`.
"""
function rows(fockmap::FockMap, restrictspace::FockSpace)::FockMap
    restrictindices::Vector{Integer} = [fockorder(fockmap|>getoutspace, mode) for mode in orderedmodes(restrictspace)]
    spmat::SparseMatrixCSC{ComplexF64, Int64} = sparse(view(rep(fockmap), restrictindices, :))
    return FockMap(restrictspace, fockmap|>getinspace, spmat)
end
export rows

"""
    rowsubmaps(fockmap::FockMap)::Base.Generator

Partition the `fockmap` into smaller `fockmaps` by the subspaces of its `outspace`.
"""
rowsubmaps(fockmap::FockMap)::Base.Generator = (rows(fockmap, subspace) for subspace in fockmap|>getoutspace |> subspaces)
export rowsubmaps

"""
    colsubmaps(fockmap::FockMap)::Base.Generator

Partition the `fockmap` into smaller `fockmaps` by the subspaces of its `inspace`.
"""
colsubmaps(fockmap::FockMap)::Base.Generator = (columns(fockmap, subspace) for subspace in fockmap|>getinspace |> subspaces)
export colsubmaps

"""
    restrict(fockmap::FockMap, outspace::FockSpace, inspace::FockSpace)::FockMap

Restrict the `outspace` & `inspace` of the `fockmap` by a sub-fockspaces `outspace` & `inspace` respectively.
"""
function restrict(fockmap::FockMap, outspace::FockSpace, inspace::FockSpace)::FockMap
    outindices::Vector{Integer} = [fockorder(fockmap|>getoutspace, mode) for mode in orderedmodes(outspace)]
    inindices::Vector{Integer} = [fockorder(fockmap|>getinspace, mode) for mode in orderedmodes(inspace)]
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

LinearAlgebra.:norm(fockmap::FockMap)::Number = norm(fockmap |> rep)
LinearAlgebra.:normalize(fockmap::FockMap)::FockMap = FockMap(fockmap |> getoutspace, fockmap |> getinspace, fockmap |> rep |> normalize)
Base.:abs(fockmap::FockMap)::FockMap = FockMap(fockmap |> getoutspace, fockmap |> getinspace, map(abs, fockmap |> rep))

""" Shorthand for retrieving the eigenvectors from the `eigspech` function. """
eigvecsh(hermitian::FockMap, attrs::Pair{Symbol}...)::FockMap = eigspech(hermitian, attrs...) |> geteigenvectors
export eigvecsh

""" Shorthand for retrieving the eigenvalues from the `eigspech` function. """
eigvalsh(hermitian::FockMap, attrs::Pair{Symbol}...)::Dict{Mode, Real} = eigspech(hermitian, attrs...) |> geteigenvalues
export eigvalsh

"""
    fourier(momentums::Subset{Momentum}, inspace::FockSpace)::FockMap

Create a `FockMap` corresponds to a Fourier transform of a `FockSpace`. This map assumed orthogonality of modes with the attribute `:offset` dropped, which means
entries corresponds to different fermionic site within the same translational invariant unit cell will be default to `0 + 0im`.

### Input
- `momentums` The momentums which spans the output reciprocal subspace, for `getspace(momentum) isa MomentumSpace`.
- `inspace` The input space, all contituent modes must have the attribute `:offset` defined or result in errors. The `inspace` should include all possible basis modes
  so that they can be identified and span the momentum `FockSpace`.

### Output
The `FockMap` represents this specific Fourier transform, with sizes `(N, M)` which `N = length(momentums) * count` for `count` is the number of fermionic site
within the translational invariant unit cell, supplied by `inmodes`; `M = length(inmodes)`. The `inspace` of the returned `FockMap` will equates to `FockSpace(inmodes)`;
the `outspace` is the product of `momentums` and the supplied fermionic sites.
"""
function fourier(momentums::Subset{Momentum}, inspace::FockSpace)::FockMap
    ∑𝑘::Matrix{Float64} = hcat([𝑘 |> euclidean |> vec for 𝑘 in momentums]...)
    inmodes::Subset{Mode} = orderedmodes(inspace)
    basismodes::Subset{Mode} = removeattr(inmodes, :offset)
    outspace::FockSpace = sparsefock(basismodes, momentums)
    return fourier(outspace, inspace, ∑𝑘, basismodes)
end
export fourier

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
    ∑𝑘::Matrix{Float64} = hcat([getattr(first(partition), :offset) |> euclidean |> vec for partition in rep(outspace)]...)
    basismodes::Subset{Mode} = removeattr(outspace |> rep |> first, :offset) # Assumed the similarity in structure for each partitions.
    return fourier(outspace, inspace, ∑𝑘, basismodes)
end

""" Internal use only. """
function fourier(outspace::FockSpace, inspace::FockSpace, momentummatrix::Matrix{Float64}, basismodes::Subset{Mode})::FockMap
    values::Array{ComplexF64} = zeros(ComplexF64, length(basismodes), size(momentummatrix, 2), dimension(inspace))
    for ((n, basismode), (m, inmode)) in Iterators.product(enumerate(basismodes), enumerate(orderedmodes(inspace)))
        if removeattr(inmode, :offset) != basismode continue end
        offset::Point = inmode |> getattr(:offset) |> euclidean
        values[n, :, m] = exp.(-1im * momentummatrix' * vec(offset))
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
        data[fockorder(outspace, outmodes |> first):fockorder(outspace, outmodes |> last), fockorder(inspace, inmodes |> first):fockorder(inspace, inmodes |> last)] += (restricted |> rep)
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
        outrange::UnitRange = fockorder(outspace, outmodes |> first):fockorder(outspace, outmodes |> last)
        inmodes::Subset{Mode} = fockmap|>getinspace |> orderedmodes
        inrange::UnitRange = fockorder(inspace, inmodes |> first):fockorder(inspace, inmodes |> last)
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
        eigenvectors::Dict{Momentum, FockMap}) = new{crystal |> dimension}(crystal, eigenmodes, eigenvalues, eigenvectors)
end
export CrystalSpectrum

Base.:show(io::IO, spectrum::CrystalSpectrum) = print(io, string("$(spectrum |> typeof)(entries=$(spectrum.eigenvalues |> length))"))

""" Shorthand to retrieve the unitcell fockspace from a `CrystalSpectrum`. """
function unitcellfock(spectrum::CrystalSpectrum)::FockSpace{Region}
    sourcefock::FockSpace = spectrum |> geteigenvectors |> first |> last |> getoutspace
    originposition::Offset = sourcefock |> first |> getattr(:pos) |> getspace |> getorigin
    return sourcefock |> orderedmodes |> setattr(:offset => originposition) |> FockSpace{Region}
end

"""
    crystalspectrum(momentumfockmaps; crystal::Crystal)::CrystalSpectrum

Given a collection of Hermitian `FockMap` objects each associated with a specific `Momentum` from the brillouin zone, pack into a
`CrystalSpectrum` object.
"""
function crystalspectrum(momentumfockmaps; crystal::Crystal)::CrystalSpectrum
    crystaleigenmodes::Dict{Momentum, Subset{Mode}} = Dict()
    crystaleigenvalues::Dict{Mode, Number} = Dict()
    crystaleigenvectors::Dict{Momentum, FockMap} = Dict()
    for (k, fockmap) in momentumfockmaps
        eigenspectrum::EigenSpectrum = eigspech(fockmap, :offset => k)
        crystaleigenmodes[k] = Subset(m for (m, _) in eigenspectrum |> geteigenvalues)
        crystaleigenvectors[k] = eigenspectrum |> geteigenvectors
        for (m, v) in eigenspectrum |> geteigenvalues
            crystaleigenvalues[m] = v
        end
    end
    return CrystalSpectrum(crystal, crystaleigenmodes, crystaleigenvalues, crystaleigenvectors)
end
export crystalspectrum

"""
    crystalspectrum(fockmap::FockMap)::CrystalSpectrum

Given a Hermitian `FockMap` with `inspace` and `outspace` of type `CrystalFock` of same span, pack into a `CrystalSpectrum` object.
"""
crystalspectrum(fockmap::FockMap)::CrystalSpectrum = crystalspectrum(fockmap |> crystalsubmaps, crystal=fockmap|>getinspace |> getcrystal)

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
    return [outmode => mat[fockorder(fockmap|>getoutspace, outmode), 1] for outmode in orderedmodes(fockmap|>getoutspace)]
end

struct RegionState{Dim} <: Element{FockMap}
    spstates::Dict{Mode, FockMap}
end
export RegionState

Base.:show(io::IO, state::RegionState) = print(io, string("$(typeof(state))(count=$(state |> getinspace |> dimension))"))

Base.:convert(::Type{FockMap}, source::RegionState)::FockMap = reduce(+, spstate for (_, spstate) in spstates)

function regionalrestriction(crystalstate::FockMap, regionfock::FockSpace)::RegionState
    eigenmodes::Subset{Mode} = crystalstate |> getinspace |> unitcellfock |> orderedmodes

    function extractregionstate(mode::Mode)
        rightfourier::FockMap = fourier(crystalstate |> getinspace, mode |> FockSpace)
        leftfourier::FockMap = fourier(crystalstate |> getoutspace, regionfock)
        return leftfourier' * crystalstate * rightfourier
    end

    return Dict(mode => mode |> extractregionstate for mode in eigenmodes) |> RegionState{crystalstate |> getoutspace |> getcrystal |> dimension}
end
export regionalrestriction

getinspace(state::RegionState) = FockSpace(m for (m, _) in state.spstates)
getoutspace(state::RegionState) = state.spstates |> first |> last |> getoutspace

Base.:iterate(state::RegionState, i...) = iterate(state.spstates, i...)
Base.:length(state::RegionState) = state.spstates |> length

struct CrystalFockMap <: FockMap{CrystalFock, CrystalFock}
    outcrystal::Crystal
    incrystal::Crystal
    outsubspaces::Dict{Momentum, FockSpace}
    insubspaces::Dict{Momentum, FockSpace}
    outbz::Subset{Momentum} # The outspace brillouin zone also served as the ordering.
    momentummappings::Dict{Momentum, Momentum} # The momentum mapping from the outspace to the inspace.
    subfockmaps::Dict{Momentum, FockMap} # The sub-fockmaps indexed by the momentums in the outspace.
end
export CrystalFockMap

function Zipper.:getoutspace(fockmap::CrystalFockMap)::CrystalFock
    fockspace::FockSpace = fockspaceunion(fockmap.outsubspaces[k] for k in fockmap.outbz)
    return FockSpace(fockspace, reflected=fockmap.outcrystal)
end

function Zipper.:getinspace(fockmap::CrystalFockMap)::CrystalFock
    fockspace::FockSpace = fockspaceunion(fockmap.insubspaces[fockmap.momentummappings[k]] for k in fockmap.outbz)
    return FockSpace(fockspace, reflected=fockmap.incrystal)
end

function CrystalFockMap(fockmap::FockMap)
    @assert(fockmap|>getoutspace isa CrystalFock)
    @assert(fockmap|>getinspace isa CrystalFock)

    outsubspaces::Base.Generator = fockmap|>getoutspace|>crystalsubspaces
    insubspaces::Base.Generator = fockmap|>getinspace|>crystalsubspaces

    outbz::Subset{Momentum} = Subset(k for (k, _) in outsubspaces)
    inbz::Subset{Momentum} = Subset(k for (k, _) in insubspaces)

    momentummappings::Dict{Momentum, Momentum} = Dict(outk => ink for (outk, ink) in Iterators.zip(outbz, inbz))
    subfockmaps::Dict = Dict(k => fockmap[outspace, inspace] for ((k, outspace), inspace) in Iterators.zip(((k, subspace) for (k, subspace) in outsubspaces), (subspace for (_, subspace) in insubspaces)))

    return CrystalFockMap(fockmap|>getoutspace|>getcrystal, fockmap|>getinspace|>getcrystal, outsubspaces|>Dict, insubspaces|>Dict, outbz, momentummappings, subfockmaps)
end

function FockMap(fockmap::CrystalFockMap)
    ret::FockMap = directsum(fockmap.subfockmaps[k] for k in fockmap.outbz)
    return FockMap(ret, inspace=fockmap|>getinspace, outspace=fockmap|>getoutspace)
end

Base.:convert(::Type{SparseMatrixCSC{ComplexF64, Int64}}, source::CrystalFockMap) = source|>FockMap|>rep

function crystalsubmaps(fockmap::CrystalFockMap)::Base.Generator
    @assert(fockmap.incrystal == fockmap.outcrystal)
    return (k => fockmap.subfockmaps[k] for k in fockmap.outbz)
end

function Base.:+(a::CrystalFockMap, b::CrystalFockMap)
    bsubmaps::Dict{Vector{Momentum}, FockMap} = Dict([k, b.momentummappings[k]] => b.subfockmaps[k] for k in b.outbz)
    addedsubmaps::Dict{Momentum, FockMap} = Dict(k => a.subfockmaps[k] + bsubmaps[[k, a.momentummappings[k]]] for k in a.outbz)
    return CrystalFockMap(a.outcrystal, a.incrystal, a.outsubspaces, a.insubspaces, a.outbz, a.momentummappings, addedsubmaps)
end

Base.:+(a::CrystalFockMap, b::FockMap)::CrystalFockMap = a + (b|>CrystalFockMap)
Base.:+(a::FockMap, b::CrystalFockMap)::FockMap = (a|>FockMap) + b

function Base.:*(a::CrystalFockMap, b::CrystalFockMap)
    momentummappings::Dict{Momentum, Momentum} = Dict(k => b.momentummappings[a.momentummappings[k]] for k in a.outbz)
    multipliedsubmaps::Dict{Momentum, FockMap} = Dict(k => a.subfockmaps[k] * b.subfockmaps[a.momentummappings[k]] for k in a.outbz)

    return CrystalFockMap(a.outcrystal, b.incrystal, a.outsubspaces, b.insubspaces, a.outbz, momentummappings, multipliedsubmaps)
end

Base.:*(a::CrystalFockMap, b::FockMap)::CrystalFockMap = a * (b|>CrystalFockMap)
Base.:*(a::FockMap, b::CrystalFockMap)::FockMap = (a|>FockMap) * b

Base.:*(num::Number, f::CrystalFockMap)::CrystalFockMap = CrystalFockMap(f.outcrystal, f.incrystal, f.outsubspaces, f.insubspaces, f.outbz, f.momentummappings, Dict(k => num * submap for (k, submap) in f.subfockmaps))
Base.:*(f::CrystalFockMap, num::Number)::CrystalFockMap = num * f
