if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end

module Quantum

using LinearAlgebra, SparseArrays, OrderedCollections
using ..Spaces, ..Geometries

export quantized, transformed, symmetrized
export ModeGroupType, ModeGroup, Mode, AnyFock, SparseFock, CrystalFock, FockSpace, FockMap
export groupname, hasattr, getattr, identify, unidentify, reidentify, setattr, removeattr, addgroup, quantize, flavorcount, spanoffset
export dimension, order, orderedmodes, orderingrule, modes, hassamespan, sparsefock, crystalfock, issparse
export columns, rows, restrict, eigvecsh, eigvalsh, eigh, fourier, focksum, idmap, columnspec

"""
    ModeGroupType

Classifiers of the type of a `ModeGroup`.
- `quantized`   The mode is being quantized from some physical objects.
- `transformed` The mode is being created by transforming existing set of modes.
- `symmetrized` The mode is being created by symmetry-transforming existing set of modes.
"""
@enum ModeGroupType begin
    quantized
    transformed
    symmetrized
end

"""
    ModeGroup(type::ModeGroupType, name::String)

A structure that holds the information of a specific group of mode.

### Input
- `type` The type of the group.
- `name` The name of the group.
"""
struct ModeGroup
    type::ModeGroupType
    name::String
end

"""
    Mode(attrs::Dict{Symbol})
    Mode(datas::Vector{Pair{Symbol, T}}) where {T}

Represents an element in a `FockSpace`, and uniquely identifies a physical mode.

### Attributes to put in `attrs`
- `:groups` stores a `Vector{ModeGroup}` which `ModeGroup` identifies a group of modes created by some action.
- `:index` stores an `Integer` that identifies the basis index.
- `:offset` stores a `Point` which is the offset in lattice unit, of this mode relative to the associated basis mode.
- `:pos` stores a `Point` which is the unit cell offset, this is associated to the attributes of `:index` and `:flavor`.
- `:flavor` stores an `Integer` that identifies a fermionic freedom at a lattice site.

### Input
- `attrs` The attributes which uniquely identifies the `Mode` object.

"""
struct Mode <: AbstractSubset{Mode}
    attrs::Dict{Symbol}

    Mode(attrs::Dict{Symbol}) = new(attrs)
    Mode(datas::Vector{Pair{Symbol, T}}) where {T} = Mode(Dict(datas...))
end

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
    removeattr(mode::Mode, keys::Symbol...)::Mode

Create a **copy** of `mode` **without** the attributes identified by `keys`.

### Examples
- To remove the attribute of `:offset` and `:pos`, we use `removeattr(mode, :offset, :pos)`.
"""
removeattr(mode::Mode, keys::Symbol...)::Mode = Mode(Dict(filter(p -> !(p.first ∈ keys), mode.attrs)))

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
- To add `:offset` and `index`, we use `newmode = setattr(mode, :offset => point, :index => 1)`
"""
setattr(mode::Mode, attrs::Pair{Symbol}...)::Mode = Mode(Dict(mode.attrs..., attrs...))

"""
    setattr(subset::Subset{Mode}, attrs::Pair{Symbol}...)::Subset{Mode}

Create a **copy** of `modes` with the new attributes identified by `keys` added to the current attributes for each mode in `modes`, if the attribute
exists in the modes, the current record will be overwritten.

### Examples
- To add `:offset` and `index`, we use `newmodes = setattr(modes, :offset => point, :index => 1)`
"""
setattr(subset::Subset{Mode}, attrs::Pair{Symbol}...)::Subset{Mode} = Subset([setattr(mode, attrs...) for mode in subset])

"""
    spanoffset(basismodes::Subset{Mode}, points::Subset{Point})::Subset{Mode}

Given a set of `basismodes`, and the generator `points`, span the basis modes to the generator `points` with attribute `:offset`, the primary ordering will be
the ordering of `points`, then follows the ordering of `basismodes`.
"""
spanoffset(basismodes::Subset{Mode}, points::Subset{Point})::Subset{Mode} = Subset([setattr(mode, :offset => point) for point in points for mode in basismodes])

"""
    flavorcount(basismodes::Subset{Mode})::Integer

Determines the number of flavor with the information provided within `basismodes`, if the provided parameter is not actually the basismodes, the return value
is just the number of distinct attribute `:flavor` within the given parameter set.
"""
flavorcount(basismodes::Subset{Mode})::Integer = trunc(Integer, length(basismodes)) / length(removeattr(basismodes, :flavor))

"""
    sparsefock(basismodes::Subset{Mode}, points::Subset{Point})::FockSpace

Given a set of `basismodes`, and the generator `points`, span the basis modes to the generator `points` with attribute `:offset` and form a `FockSpace`. Not to
be mistakened with `spanoffset`, this method will partition the modes by the generator points, in normal conditions the return type will be `FockSpace{SparseFock}`.
Noted that the ordering of the partitions will follow the ordering of `points`, and the ordering within each partition will follow the ordering of `basismodes`.
"""
function sparsefock(basismodes::Subset{Mode}, points::Subset{Point})::FockSpace
    partitions::Vector{Subset{Mode}} = [setattr(basismodes, :offset => point) for point in points]
    modes::Subset{Mode} = spanbasis(basismodes, points)
    orderings::Dict{Mode, Integer} = Dict(mode => index for (index, mode) in enumerate(modes))
    return FockSpace(Subset(partitions::Vector{Subset{Mode}}), orderings)
end

""" By this conversion, one can obtain the actual position of the mode, this method only works when `:offset` and `:pos` are defined in the same space. """
Base.:convert(::Type{Point}, source::Mode)::Point = getattr(source, :offset) + getattr(source, :pos)

Base.:(==)(a::Mode, b::Mode)::Bool = a.attrs == b.attrs

# =====================================================
# Mode can be used as keys in dictionaries and sets.
Base.:hash(mode::Mode)::UInt = hash(mode.attrs)
Base.:isequal(a::Mode, b::Mode) = a == b
# =====================================================

"""
    FockSpace(subsets::Subset{Subset{Mode}}, ordering::Dict{Mode, <: Integer}; T::Type{<: AnyFock})
    FockSpace(subset::Subset{Mode}; T::Type{<: AnyFock})
    FockSpace(fockspace::FockSpace; T::Type{<: AnyFock} = AnyFock)

A collection of `Modes` or `Mode` partitions in the case of sparse fockspace, implicit ordering of underlying modes is assumed.

The structure of fockspace is assume to hold partitions, a fockspace that **does not** have partition is assumed to have a single partition containing
all underlying modes. The design can be seen with the type of the representation `Subset{Subset{Mode}}`, which the higher order `Subset` holds partitions
represented by `Subset{Mode}`.

The type of the fockspace will be determined by the constructor with the structure of the representation, if there is more than one partition,
`SparseFock` will be selected, else `AnyFock` is the default option.

### Examples
- `FockSpace(subsets::Subset{Subset{Mode}}, ordering::Dict{Mode, <: Integer}; T::Type{<: AnyFock})` is used when all the components of the `FockSpace`
  is already constructed prior instantiation.
- `FockSpace(subset::Subset{Mode}; T::Type{<: AnyFock})` is the normal use case to convert a `Subset{Mode}` into `FockSpace`.
- `FockSpace(fockspace::FockSpace; T::Type{<: AnyFock} = AnyFock)` is used to change the fockspace type to `AnyFock`, `SparseFock` or `CrystalFock`.
"""
struct FockSpace{T} <: AbstractSpace{Subset{Subset{Mode}}}
    rep::Subset{Subset{Mode}}
    ordering::Dict{Mode, Integer}

    FockSpace(subsets::Subset{Subset{Mode}}, ordering::Dict{Mode, <: Integer};
        T::Type{<: AnyFock} = (length(subsets) > 1 ? SparseFock : AnyFock)) = new{T}(subsets, ordering)
    FockSpace(subset::Subset{Mode}; T::Type{<: AnyFock} = AnyFock) = FockSpace(
        Subset([subset]),
        Dict(mode => order for (order, mode) in enumerate(subset)),
        T=T)
    FockSpace(fockspace::FockSpace; T::Type{<: AnyFock} = AnyFock) = FockSpace(rep(fockspace), fockspace.ordering, T=T)
end

function Base.:union(focks::FockSpace...)::FockSpace
    partitions::Subset{Subset{Mode}} = union([rep(fock) for fock in focks]...)
    modes::Subset{Mode} = flatten(partitions)
    ordering::Dict{Mode, Integer} = Dict(mode => index for (index, mode) in enumerate(modes))
    return FockSpace(partitions, ordering)
end

Base.:iterate(fock_space::FockSpace, i...) = iterate(flatten(rep(fock_space)), i...)

"""
    Spaces.:dimension(fockspace::FockSpace)

Returns the number of unique member modes within the `fockspace`, each of those represents a vector from the Hilbert space.
"""
Spaces.:dimension(fockspace::FockSpace) = length(fockspace.ordering) # This is a short cut to retrieve the length, which is the dimension.

"""
    order(fockspace::FockSpace, mode::Mode)::Int64

Since the fockspace have implicit ordering, this function returns the order index of the `mode` in the `fockspace`.
"""
order(fockspace::FockSpace, mode::Mode)::Int64 = fockspace.ordering[mode]

"""
    modes(fockspace::FockSpace)::Set{Mode}

Returns an unordered set of modes of `fockspace`, this is a more efficient way to retrieve the underlying modes for a fockspace of type `SparseFock`.
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
The returned vector that at each mode position `n` of `tospace`, contains the order index a mode in `fromspace` to be moved/permuted into slot `n`, the mode ordering
of the permuted list of modes will matches the ordering of `tospace`.
"""
function orderingrule(fromspace::FockSpace, tospace::FockSpace)::Vector{Int64}
    @assert(hassamespan(fromspace, tospace))
    return [tospace.ordering[mode] for mode in orderedmodes(fromspace)]
end

"""
    hassamespan(a::FockSpace, b::FockSpace)::Bool

Check if fockspaces of `a` and `b` shares the same set of underlying modes regardless of partitions.
"""
hassamespan(a::FockSpace, b::FockSpace)::Bool = modes(a) == modes(b)

Base.:(==)(a::FockSpace, b::FockSpace)::Bool = rep(a) == rep(b)

Base.:convert(::Type{Subset}, source::FockSpace) = source.rep
Base.:convert(::Type{Subset{Subset}}, source::FockSpace) = convert(Subset, source)
Base.:convert(::Type{Subset{Subset{Mode}}}, source::FockSpace) = convert(Subset, source)
"""
    quantize(index::Integer, identifier::Symbol, point::Point, flavor::Integer; group::ModeGroup)::Mode

Quantizing a mode from a given `Point`.

### Input
- `index` The index of the `Mode` within it's own group.
- `identifier` The identifying atttibute key which the `point` object will be linked to.
- `point` The `Point` as the physical attribute or object to be quantized.
- `flavor` The flavor index of the `Mode`, don't mistaken this with the flavor count.
- `group` Optional parameter that stores the information about the origins of the mode group, or the actions they have been through,
  defaults to `ModeGroup(quantized, "physical")`.

### Output
The quantized `Mode` object.
"""
    @assert(identifier == :offset || identifier == :pos)
    home::Point = origin(spaceof(point))
    # Since there are only one of the attribute :offset or :pos will take the point, the left over shall take the origin.
    couple::Pair{Symbol, Point} = identifier == :offset ? :pos => home : :offset => home
    # The new mode will take a group of q:$(name).
    return Mode([:groups => [ModeGroup(quantized, name)], :index => index, identifier => point, :flavor => flavor, couple])
end

"""
    quantize(identifier::Symbol, subset::Subset{Point}, count::Integer; group::ModeGroup)::Subset{Mode}

Quantizing a set of mode from a given set of `Point`.

### Input
- `identifier` The identifying atttibute key which the `point` object will be linked to.
- `subset` The set of `Point` provided as the physical attributes or objects to be quantized.
- `count` The flavor count of the quantization, if it is greater than `1`, it means the given site defined by a `Point` in `subset` has more
  than one fermionic degree of freedom.
- `group` Optional parameter that stores the information about the origins of the mode group, or the actions they have been through,
  defaults to `ModeGroup(quantized, "physical")`.

### Output
The quantized set of
"""

struct FockMap <: Element{SparseMatrixCSC{ComplexF64, Int64}}
    outspace::FockSpace
    inspace::FockSpace
    rep::SparseMatrixCSC{ComplexF64, Int64}

    FockMap(outspace::FockSpace, inspace::FockSpace, rep::SparseMatrixCSC{ComplexF64, Int64}) = new(outspace, inspace, rep)
    FockMap(outspace::FockSpace, inspace::FockSpace, rep::AbstractArray{<:Number}) = new(outspace, inspace, SparseMatrixCSC{ComplexF64, Int64}(rep))

    function FockMap(outspace::FockSpace, inspace::FockSpace, mapping::Dict{Tuple{Mode, Mode}, ComplexF64})::FockMap
        rep::SparseMatrixCSC{ComplexF64, Int64} = spzeros(dimension(outspace), dimension(inspace))
        for ((out_mode::Mode, in_mode::Mode), value::ComplexF64) in mapping
            rep[outspace.ordering[out_mode], inspace.ordering[in_mode]] = value
        end
        return new(outspace, inspace, rep)
    end
end

function idmap(outspace::FockSpace, inspace::FockSpace)::FockMap
    @assert(dimension(outspace) == dimension(inspace))
    FockMap(outspace, inspace, SparseMatrixCSC(Matrix{Float64}(I(dimension(outspace)))))
end

Base.:convert(::Type{SparseMatrixCSC{ComplexF64, Int64}}, source::FockMap) = source.rep

function columns(fockmap::FockMap, restrictspace::FockSpace)::FockMap
    restrictindices::Vector{Integer} = [order(fockmap.inspace, mode) for mode in orderedmodes(restrictspace)]
    spmat::SparseMatrixCSC{ComplexF64, Int64} = sparse(view(rep(fockmap), :, restrictindices))
    return FockMap(fockmap.outspace, restrictspace, spmat)
end

function rows(fockmap::FockMap, restrictspace::FockSpace)::FockMap
    restrictindices::Vector{Integer} = [order(fockmap.outspace, mode) for mode in orderedmodes(restrictspace)]
    spmat::SparseMatrixCSC{ComplexF64, Int64} = sparse(view(rep(fockmap), restrictindices, :))
    return FockMap(restrictspace, fockmap.inspace, spmat)
end

function restrict(fockmap::FockMap, out_restrictspace::FockSpace, in_restrictspace::FockSpace)::FockMap
    out_restrictindices::Vector{Integer} = [order(fockmap.outspace, mode) for mode in orderedmodes(out_restrictspace)]
    in_restrictindices::Vector{Integer} = [order(fockmap.inspace, mode) for mode in orderedmodes(in_restrictspace)]
    spmat::SparseMatrixCSC{ComplexF64, Int64} = sparse(view(rep(fockmap), out_restrictindices, in_restrictindices))
    return FockMap(out_restrictspace, in_restrictspace, spmat)
end

function permute(source::FockMap, outspace::FockSpace=source.outspace, inspace::FockSpace=source.inspace)::FockMap
    row_rule::Vector{Int64} = orderingrule(source.outspace, outspace)
    col_rule::Vector{Int64} = orderingrule(source.inspace, inspace)
    return FockMap(outspace, inspace, SparseArrays.permute(rep(source), row_rule, col_rule))
end

Base.:-(target::FockMap)::FockMap = FockMap(target.outspace, target.inspace, -rep(target))

Base.:+(a::FockMap, b::FockMap)::FockMap = focksum([a, b])

Base.:-(a::FockMap, b::FockMap)::FockMap = a + (-b)

function Base.:*(a::FockMap, b::FockMap)::FockMap
    @assert(hassamespan(a.inspace, b.outspace)) # Even if the fockspaces are different, composition works as long as they have same span.
    return FockMap(a.outspace, b.inspace, rep(a) * rep(permute(b, a.inspace, b.inspace)))
end

Base.:*(fockmap::FockMap, number::Number)::FockMap = FockMap(fockmap.outspace, fockmap.inspace, rep(fockmap) * number)
Base.:*(number::Number, fockmap::FockMap)::FockMap = fockmap * number

Base.:/(fockmap::FockMap, number::Number)::FockMap = FockMap(fockmap.outspace, fockmap.inspace, rep(fockmap) / number)

Base.:transpose(source::FockMap)::FockMap = FockMap(source.inspace, source.outspace, rep(source)')
Base.:adjoint(source::FockMap)::FockMap = FockMap(source.inspace, source.outspace, rep(source)')

function eigmodes(fockmap::FockMap, attrs::Pair{Symbol}...)::Subset{Mode}
    @assert(fockmap.inspace == fockmap.outspace)
    return Subset([
        Mode([:groups => [ModeGroup(transformed, "eigh")], :index => index, :flavor => 1, attrs...]) for index in 1:dimension(fockmap.inspace)])
end

function eigvecsh(fockmap::FockMap, attrs::Pair{Symbol}...)::FockMap
    @assert(fockmap.inspace == fockmap.outspace)
    evecs::Matrix{ComplexF64} = eigvecs(Hermitian(Matrix(rep(fockmap))))
    return FockMap(fockmap.outspace, FockSpace(eigmodes(fockmap, attrs...)), SparseMatrixCSC{ComplexF64, Int64}(evecs))
end

function eigvalsh(fockmap::FockMap, attrs::Pair{Symbol}...)::Vector{Pair{Mode, Float64}}
    @assert(fockmap.inspace == fockmap.outspace)
    evs::Vector{Number} = eigvals(Hermitian(Matrix(rep(fockmap))))
    return [first(tup) => last(tup) for tup in Iterators.zip(eigmodes(fockmap, attrs...), evs)]
end

function eigh(fockmap::FockMap, attrs::Pair{Symbol}...)::Tuple{Vector{Pair{Mode, Float64}}, FockMap}
    decomposed::Eigen = eigen(Hermitian(Matrix(rep(fockmap))))
    modes::Subset{Mode} = eigmodes(fockmap, attrs...)
    evals::Vector{Pair{Mode, Float64}} = [m => v for (m, v) in Iterators.zip(modes, decomposed.values)]
    evecs::FockMap = FockMap(fockmap.outspace, FockSpace(modes), decomposed.vectors)
    return evals, evecs
end

"""
    fourier(momentums::Subset{Point}, inmodes::Subset{Mode})::FockMap

Create a `FockMap` corresponds to a Fourier transform of a `Subset{Mode}`. This map assumed orthogonality of modes with the attribute `:offset` dropped, which means
entries corresponds to different fermionic site within the same translational invariant unit cell will be default to `0 + 0im`.

### Input
- `momentums` The momentums which spans the output reciprocal subspace, for `spaceof(momentum) isa MomentumSpace`.
- `inspace` The input space, all contituent modes must have the attribute `:offset` defined or result in errors.

### Output
The `FockMap` represents this specific Fourier transform, with sizes `(N, M)` which `N = length(momentums) * count` for `count` is the number of fermionic site
within the translational invariant unit cell, supplied by `inmodes`; `M = length(inmodes)`. The `inspace` of the returned `FockMap` will equates to `FockSpace(inmodes)`;
the `outspace` is the product of `momentums` and the supplied fermionic sites.
"""
function fourier(momentums::Subset{Point}, inspace::FockSpace)::FockMap
    ∑𝑘::Matrix{Float64} = hcat([pos(euclidean(𝑘)) for 𝑘 in momentums]...)
    inmodes::Subset{Mode} = orderedmodes(inspace)
    basismodes::Subset{Mode} = removeattr(inmodes, :offset)
    outspace::FockSpace = sparsefock(basismodes, momentums)
    return fourier(outspace, inspace, ∑𝑘, basismodes)
end

function fourier(outspace::FockSpace, inspace::FockSpace)::FockMap
    ∑𝑘::Matrix{Float64} = hcat([pos(euclidean(getattr(first(partition), :offset))) for partition in rep(outspace)]...)
    basismodes::Subset{Mode} = removeattr(first(rep(outspace)), :offset) # Assumed the similarity in structure for each partitions.
    return fourier(outspace, inspace, ∑𝑘, basismodes)
end

function fourier(outspace::FockSpace, inspace::FockSpace, momentummatrix::Matrix{Float64}, basismodes::Subset{Mode})::FockMap
    values::Array{ComplexF64} = zeros(ComplexF64, length(basismodes), size(momentummatrix, 2), dimension(inspace))
    for ((n, basismode), (m, inmode)) in Iterators.product(enumerate(basismodes), enumerate(orderedmodes(inspace)))
        if removeattr(inmode, :offset) != basismode continue end
        𝑟ₑ::Point = euclidean(getattr(inmode, :offset))
        values[n, :, m] = exp.(-1im * momentummatrix' * pos(𝑟ₑ))
    end
    spmat::SparseMatrixCSC = SparseMatrixCSC(reshape(values, (length(basismodes) * size(momentummatrix, 2), dimension(inspace))))
    return FockMap(outspace, inspace, spmat)
end

function focksum(fockmaps::Vector{FockMap})::FockMap
    parts = [(fockmap, outpart, inpart) for fockmap in fockmaps for (outpart, inpart) in zip(rep(fockmap.outspace), rep(fockmap.inspace))]
    outfock = union([fockmap.outspace for fockmap in fockmaps]...)
    infock = union([fockmap.inspace for fockmap in fockmaps]...)
    spmat::SparseMatrixCSC{ComplexF64, Int64} = spzeros(dimension(outfock), dimension(infock))
    for (fockmap, outpart, inpart) in parts
        source::SparseMatrixCSC{ComplexF64, Int64} = rep(fockmap)
        for (outmode, inmode) in Iterators.product(outpart, inpart)
            spmat[order(outfock, outmode), order(infock, inmode)] += source[order(fockmap.outspace, outmode), order(fockmap.inspace, inmode)]
        end
    end
    return FockMap(outfock, infock, spmat)
end

function columnspec(fockmap::FockMap)::Vector{Pair{Mode, ComplexF64}}
    @assert(dimension(fockmap.inspace) == 1)
    mat::SparseMatrixCSC{ComplexF64, Int64} = rep(fockmap)
    return [outmode => mat[order(fockmap.outspace, outmode), 1] for outmode in orderedmodes(fockmap.outspace)]
end

end
