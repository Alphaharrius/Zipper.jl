if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end

module Quantum

using LinearAlgebra, SparseArrays
using ..Spaces, ..Geometries

const MODE_IDENTIFIERS_TYPE::Base.ImmutableDict{Symbol, DataType} = Base.ImmutableDict(
    :groups => Vector{String},
    :index => Integer,
    :offset => Point,
    :pos => Point,
    :flavor => Integer)
const MODE_DEFAULT_IDENTIFIERS::Set{Symbol} = Set([:groups, :index, :pos, :flavor])

groupname(type::Symbol, name::String)::String = "$(type):$(name)"

struct Mode <: AbstractSubset{Mode}
    attrs::Dict{Symbol}

    Mode(attrs::Dict{Symbol}) = new(attrs)
    Mode(datas::Vector{Pair{Symbol, T}}) where {T} = Mode(Dict(datas...))
end

"""
The space of a `Mode` comes from the physical quantities its defined on, such as `:offset` and `:pos`, if none of those are defined,
it will be `euclidean(RealSpace, 1)` as the mode and it's siblings can always be parameterized by a scalar.
"""
function Spaces.spaceof(mode::Mode)::AbstractSpace
    if hasattr(mode, :pos) return spaceof(getattr(mode, :pos)) end
    if hasattr(mode, :offset) return spaceof(getattr(mode, :offset)) end
    # If the mode does not based on any physical position or quantities for associating with any space, then it will simply lives
    # in a 1D euclidean space as the mode and it's siblings can always be parameterized by a scalar.
    return euclidean(RealSpace, 1)
end

hasattr(mode::Mode, key::Symbol) = haskey(mode.attrs, key)

getattr(mode::Mode, key::Symbol) = mode.attrs[key]

removeattr(mode::Mode, keys::Symbol...) = Mode(Dict(filter(p -> !(p.first ∈ keys), mode.attrs)))

function setattr(mode::Mode, attrs::Pair{Symbol, T}...)::Mode where {T}
    # TODO: The @assert is removed temporarily.
    return Mode(Dict(mode.attrs..., attrs...))
end

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

Base.:length(fock_space::FockSpace) = length(fock_space.ordering)

function ordered_modes(fock_space::FockSpace)::Array{Mode}
    modes = Array{Mode}(undef, dimension(fock_space))
    for mode in flatten(rep(fock_space))
        modes[fock_space.ordering[mode]] = mode
    end
    return modes
end

function ordering_rule(from_space::FockSpace, to_space::FockSpace)::Vector{Int64}
    @assert(from_space == to_space)
    return [to_space.ordering[mode] for mode in ordered_modes(from_space)]
end

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
    return Mode([:groups => [groupname(:q, name)], :index => index, identifier => point, :flavor => flavor, couple])
end

quantize(name::String, identifier::Symbol, subset::Subset{Point}, count::Integer)::Subset{Mode} = (
    Subset(Set([quantize(name, index, identifier, point, flavor) for (index, point) in enumerate(subset) for flavor in 1:count])))

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

function columns(fock_map::FockMap, subset::Subset{Mode})::FockMap
    inspace::FockSpace = FockSpace(subset)
    matrix::SparseMatrixCSC{ComplexF64, Int64} = spzeros(dimension(fock_map.outspace), dimension(inspace))
    for mode in rep(subset)
        matrix[:, inspace.ordering[mode]] = rep(fock_map)[:, fock_map.outspace.ordering[mode]]
    end
    return FockMap(inspace, fock_map.outspace, matrix)
end

function permute(source::FockMap, outspace::FockSpace=source.outspace, inspace::FockSpace=source.inspace)::FockMap
    row_rule::Vector{Int64} = ordering_rule(source.outspace, outspace)
    col_rule::Vector{Int64} = ordering_rule(source.inspace, inspace)
    return FockMap(outspace, inspace, SparseArrays.permute(rep(source), row_rule, col_rule))
end

function Base.:+(a::FockMap, b::FockMap)::FockMap
    @assert(a.inspace == b.inspace && a.outspace == b.outspace)
    return FockMap(a.outspace, a.inspace, rep(a) + rep(permute(b, a.outspace, a.inspace)))
end

function Base.:-(a::FockMap, b::FockMap)::FockMap
    @assert(a.inspace == b.inspace && a.outspace == b.outspace)
    return FockMap(a.outspace, a.inspace, rep(a) - rep(permute(b, a.outspace, a.inspace)))
end

function Base.:*(a::FockMap, b::FockMap)::FockMap
    @assert(a.inspace == b.outspace)
    return FockMap(a.outspace, b.inspace, rep(a) * rep(permute(b, a.inspace, b.inspace)))
end

transpose(source::FockMap)::FockMap = FockMap(source.inspace, source.outspace, rep(source)')

dagger(source::FockMap)::FockMap = FockMap(source.inspace, source.outspace, rep(source)') # `'` operator is dagger.

function eigvalsh(fock_map::FockMap, attrs::Pair{Symbol, T}...)::Vector{Pair{Mode, Float64}} where {T}
    @assert(fock_map.inspace == fock_map.outspace)
    evs::Vector{Number} = eigvals(Hermitian(Matrix(rep(fock_map))))
    return [Mode([:groups => [groupname(:t, "eigh")], :index => index, :flavor => 1, attrs...]) => ev for (index, ev) in enumerate(evs)]
end

"""
    fourier(momentums::Subset{Point}, inmodes::Subset{Mode})::FockMap

Create a `FockMap` corresponds to a Fourier transform of a `Subset{Mode}`. This map assumed orthogonality of modes with the attribute `:offset` dropped, which means
entries corresponds to different fermionic site within the same translational invariant unit cell will be default to `0 + 0im`.

### Input
- `momentums` The momentums which spans the output reciprocal subspace, for `spaceof(momentum) isa MomentumSpace`.
- `inmodes` The modes from the input subspace, all modes must have the attribute `:offset` defined or result in errors.

### Output
The `FockMap` represents this specific Fourier transform, with sizes `(N, M)` which `N = length(momentums) * count` for `count` is the number of fermionic site
within the translational invariant unit cell, supplied by `inmodes`; `M = length(inmodes)`. The `inspace` of the returned `FockMap` will equates to `FockSpace(inmodes)`;
the `outspace` is the product of `momentums` and the supplied fermionic sites.
"""
function fourier(momentums::Subset{Point}, inmodes::Subset{Mode})::FockMap
    # Disable the identification by :offset so that all inmodes collapes to the basis mode. 
    basismodes::Set{Mode} = Set(removeattr(inmode, :offset) for inmode in inmodes)
    # Enable the identification by :offset when adding momentum as :offset for each basis mode.
    outmodes = [Iterators.map(tup -> setattr(tup[1], :offset => tup[2]), Iterators.product(basismodes, momentums))...]
    outspace_orderings::Vector{Pair{Mode, Integer}} = []
    inspace_orderings::Vector{Pair{Mode, Integer}} = []
    mat::SparseMatrixCSC{ComplexF64, Int64} = spzeros(length(outmodes), length(inmodes))
    for ((oi, outmode), (ii, inmode)) in Iterators.product(enumerate(outmodes), enumerate(inmodes))
        # Different fermionic site within the same translational invariant unit cell will not be considered.
        if removeattr(outmode, :offset) != removeattr(inmode, :offset) continue end
        momentum::Point = getattr(outmode, :offset)
        euc_momentum::Point = linear_transform(euclidean(MomentumSpace, length(momentum)), momentum)
        inoffset::Point = getattr(inmode, :offset)
        euc_inoffset::Point = linear_transform(euclidean(RealSpace, length(inoffset)), inoffset)
        mat[oi, ii] = exp(-1im * dot(euc_momentum, euc_inoffset))
        push!(outspace_orderings, outmode => oi)
        push!(inspace_orderings, inmode => ii)
    end
    return FockMap(
        FockSpace(Subset([Subset(outmodes)]), Dict(outspace_orderings...)),
        FockSpace(Subset([inmodes]), Dict(inspace_orderings...)),
        mat)
end

export Mode, FockSpace, FockMap
export groupname, hasattr, getattr, identify, unidentify, reidentify, setattr, removeattr, addgroup
export dimension, ordered_modes, ordering_rule, quantize, columns, transpose, dagger, eigvalsh, fourier

end
