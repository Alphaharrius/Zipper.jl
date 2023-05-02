if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end

module Quantum

using LinearAlgebra, SparseArrays, OrderedCollections
using ..Spaces, ..Geometries

export quantized, transformed, symmetrized
export ModeGroupType, ModeGroup, Mode, FockSpace, FockMap
export groupname, hasattr, getattr, identify, unidentify, reidentify, setattr, removeattr, addgroup
export dimension, orderedmodes, orderingrule, modes, hassamespan, quantize, columns, rows, restrict, eigvecsh, eigvalsh, eigh, fourier, focksum
export hermitian_partitioning, columnspec, spanbasis

"""
    ModeGroupType

Classifiers of the type of a `ModeGroup`.
"""
@enum ModeGroupType begin
    quantized
    transformed
    symmetrized
end


"""
    ModeGroup(type::ModeGroupType, name::String)

A structure that holds the information of a specific group of mode.
"""
struct ModeGroup
    "The type of the group."
    type::ModeGroupType
    "The name of the group."
    name::String
end

const MODE_IDENTIFIERS_TYPE::Base.ImmutableDict{Symbol, DataType} = Base.ImmutableDict(
    :groups => Vector{ModeGroup},
    :index => Integer,
    :offset => Point,
    :pos => Point,
    :flavor => Integer)
const MODE_DEFAULT_IDENTIFIERS::Set{Symbol} = Set([:groups, :index, :pos, :flavor])

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

hasattr(mode::Mode, key::Symbol) = haskey(mode.attrs, key)

getattr(mode::Mode, key::Symbol) = mode.attrs[key]

removeattr(mode::Mode, keys::Symbol...) = Mode(Dict(filter(p -> !(p.first ∈ keys), mode.attrs)))

removeattr(modes::Subset{Mode}, keys::Symbol...) = Subset(OrderedSet{Mode}(removeattr(mode, keys...) for mode in modes))

function setattr(mode::Mode, attrs::Pair{Symbol}...)::Mode
    # TODO: The @assert is removed temporarily.
    return Mode(Dict(mode.attrs..., attrs...))
end

setattr(subset::Subset{Mode}, attrs::Pair{Symbol}...)::Subset{Mode} = Subset([setattr(mode, attrs...) for mode in subset])

spanbasis(basismodes::Subset{Mode}, points::Subset{Point})::Subset{Mode} = Subset([setattr(mode, :offset => point) for point in points for mode in basismodes])

Base.:convert(::Type{Point}, source::Mode)::Point = getattr(source, :offset) + getattr(source, :pos)

Base.:(==)(a::Mode, b::Mode)::Bool = a.attrs == b.attrs

Base.:hash(mode::Mode)::UInt = hash(mode.attrs)
Base.:isequal(a::Mode, b::Mode) = a == b

struct FockSpace <: AbstractSpace{Subset{Subset{Mode}}}
    rep::Subset{Subset{Mode}}
    """ Implicit ordering is required as matrices is used for mapping between `FockSpace`. """
    ordering::Dict{Mode, Integer}

    FockSpace(subsets::Subset{Subset{Mode}}, ordering::Dict{Mode, <: Integer}) = new(subsets, ordering)
    FockSpace(subset::Subset{Mode}) = FockSpace(
        Subset([subset]),
        Dict([mode => order for (order, mode) in enumerate(subset)]...))
end

Base.:iterate(fock_space::FockSpace, i...) = iterate(flatten(rep(fock_space)), i...)

Spaces.:dimension(fockspace::FockSpace) = length(fockspace.ordering)

order(fockspace::FockSpace, mode::Mode)::Int64 = fockspace.ordering[mode]

modes(fockspace::FockSpace)::Set{Mode} = Set(keys(fockspace.ordering)) # This is the most efficient way to get all distinct modes.

orderedmodes(fockspace::FockSpace)::Subset{Mode} = flatten(rep(fockspace))

function orderingrule(fromspace::FockSpace, tospace::FockSpace)::Vector{Int64}
    @assert(hassamespan(fromspace, tospace))
    return [tospace.ordering[mode] for mode in orderedmodes(fromspace)]
end

hassamespan(a::FockSpace, b::FockSpace)::Bool = modes(a) == modes(b)

Base.:(==)(a::FockSpace, b::FockSpace)::Bool = rep(a) == rep(b)

Base.:convert(::Type{Subset}, source::FockSpace) = source.rep
Base.:convert(::Type{Subset{Subset}}, source::FockSpace) = convert(Subset, source)
Base.:convert(::Type{Subset{Subset{Mode}}}, source::FockSpace) = convert(Subset, source)

Base.:convert(::Type{FockSpace}, source::Subset{Mode}) = FockSpace(source)

function quantize(name::String, index::Integer, identifier::Symbol, point::Point, flavor::Integer)::Mode
    @assert(identifier == :offset || identifier == :pos)
    home::Point = origin(spaceof(point))
    # Since there are only one of the attribute :offset or :pos will take the point, the left over shall take the origin.
    couple::Pair{Symbol, Point} = identifier == :offset ? :pos => home : :offset => home
    # The new mode will take a group of q:$(name).
    return Mode([:groups => [ModeGroup(quantized, name)], :index => index, identifier => point, :flavor => flavor, couple])
end

quantize(name::String, identifier::Symbol, subset::Subset{Point}, count::Integer)::Subset{Mode} = (
    Subset([quantize(name, index, identifier, point, flavor) for (index, point) in enumerate(subset) for flavor in 1:count]))

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
    𝑅ₖ::MomentumSpace = euclidean(MomentumSpace, dimension(first(momentums)))
    ∑𝑘::Matrix{Float64} = hcat([pos(lineartransform(𝑅ₖ, 𝑘)) for 𝑘 in momentums]...)

    inmodes::Subset{Mode} = orderedmodes(inspace)
    𝑅::RealSpace = euclidean(RealSpace, dimension(first(momentums)))
    basismodes::Subset{Mode} = removeattr(inmodes, :offset)
    values::Array{ComplexF64} = zeros(ComplexF64, length(basismodes), length(momentums), dimension(inspace))
    for ((n, basismode), (m, inmode)) in Iterators.product(enumerate(basismodes), enumerate(inmodes))
        if removeattr(inmode, :offset) != basismode continue end
        𝑟ₑ::Point = lineartransform(𝑅, getattr(inmode, :offset))
        values[n, :, m] = exp.(-1im * ∑𝑘' * pos(𝑟ₑ))
    end

    spmat::SparseMatrixCSC = SparseMatrixCSC(reshape(values, (length(basismodes) * length(momentums), dimension(inspace))))
    outspace::FockSpace = FockSpace(Subset(OrderedSet{Mode}(setattr(m, :offset => 𝑘) for 𝑘 in momentums for m in basismodes)))
    return FockMap(outspace, inspace, spmat)
end

function focksum(fockmaps::Vector{FockMap})::FockMap
    outparts::Vector{Subset{Mode}} = [part for fockmap in fockmaps for part in rep(fockmap.outspace)]
    inparts::Vector{Subset{Mode}} = [part for fockmap in fockmaps for part in rep(fockmap.inspace)]
    outmodes::OrderedSet{Mode} = OrderedSet([mode for subset in outparts for mode in subset])
    inmodes::OrderedSet{Mode} = OrderedSet([mode for subset in inparts for mode in subset])
    outordering::Dict{Mode, Int64} = Dict(map(tup -> last(tup) => first(tup), enumerate(outmodes)))
    inordering::Dict{Mode, Int64} = Dict(map(tup -> last(tup) => first(tup), enumerate(inmodes)))
    spmat::SparseMatrixCSC{ComplexF64, Int64} = spzeros(length(outmodes), length(inmodes))
    for (outpart, inpart, fockmap) in Iterators.zip(outparts, inparts, fockmaps)
        source::SparseMatrixCSC{ComplexF64, Int64} = rep(fockmap)
        for (outmode, inmode) in Iterators.product(outpart, inpart)
            spmat[outordering[outmode], inordering[inmode]] += source[order(fockmap.outspace, outmode), order(fockmap.inspace, inmode)]
        end
    end
    return FockMap(FockSpace(Subset([outparts...]), outordering), FockSpace(Subset([inparts...]), inordering), spmat)
end

"""
    hermitian_partitioning(hermitian::FockMap, partitions::Pair{Symbol}...)::Dict{Symbol, FockMap}

Partition the Hermitian `FockMap` into column partitions by a given set of rules based on the eigenvalues from eigenvalue decompositions.

### Input
- `hermitian` The source `FockMap` to be partitioned, this is assumed to be a Hermitian.
- `partition_rules` The rules for each partition in form of `Pair{Symbol, T}` for `T` is a predicate function with input of type `<: Float`,
                    each rule `Pair` have the first `Symbol` typed element as the partition name.

### Output
A `Dict{Symbol, FockMap}` keyed by the partition name specified from the `partition_rules` with value of the associated partitioned `FockMap`.
"""
function hermitian_partitioning(hermitian::FockMap, partition_rules::Pair{Symbol}...)::Dict{Symbol, FockMap}
    eigenvalues::Vector{Pair{Mode, Float64}} = eigvalsh(hermitian)
    𝔘::FockMap = eigvecsh(hermitian)
    ret::Dict{Symbol, FockMap} = Dict()
    for rule in partition_rules
        modes::Vector{Mode} = map(p -> p.first, filter(p -> rule.second(p.second), eigenvalues))
        ret[rule.first] = columns(𝔘, FockSpace(Subset(modes)))
    end
    return ret
end

function columnspec(fockmap::FockMap)::Vector{Pair{Mode, ComplexF64}}
    @assert(dimension(fockmap.inspace) == 1)
    mat::SparseMatrixCSC{ComplexF64, Int64} = rep(fockmap)
    return [outmode => mat[order(fockmap.outspace, outmode), 1] for outmode in orderedmodes(fockmap.outspace)]
end

end
