if !isdefined(Main, :Spaces) include("spaces.jl") end

module Quantum

using LinearAlgebra, SparseArrays, Arpack
using ..Spaces

MODE_IDENTIFIERS_TYPE::Base.ImmutableDict{Symbol, DataType} = Base.ImmutableDict(
    :groups => Vector{String},
    :index => Integer,
    :offset => Point,
    :pos => Point,
    :flavor => Integer)
MODE_DEFAULT_IDENTIFIERS::Set{Symbol} = Set([:groups, :index, :pos, :flavor])

groupname(type::Symbol, name::String)::String = "$(type):$(name)"

struct Mode <: AbstractSubset{Mode}
    attributes::Base.ImmutableDict{Symbol}
    identifiers::Set{Symbol}

    Mode(attributes::Base.ImmutableDict{Symbol, T}, identifiers::Set{Symbol}) where {T} = new(attributes, identifiers)

    function Mode(datas::Vector{Pair{Symbol, T}}, identifiers::Set{Symbol} = MODE_DEFAULT_IDENTIFIERS) where {T}
        @assert(any([data.second isa MODE_IDENTIFIERS_TYPE[data.first] for data in datas]))
        attributes = Base.ImmutableDict(datas...)
        @assert(any([haskey(attributes, identifier) for identifier in identifiers]))
        return new(attributes, identifiers)
    end
end

function Spaces.space_of(mode::Mode)::AbstractSpace
    if hasattr(mode, :pos) return space_of(getattr(mode, :pos)) end
    if hasattr(mode, :offset) return space_of(getattr(mode, :offset)) end
    # If the mode does not based on any physical position or quantities for associating with any space, then it will simply lives
    # in a 1D euclidean space as the mode and it's siblings can always be parameterized by a scalar.
    return euclidean(RealSpace, 1)
end

hasattr(mode::Mode, identifier::Symbol) = haskey(mode.attributes, identifier)

getattr(mode::Mode, identifier::Symbol) = mode.attributes[identifier]

identify(mode::Mode, identifier::Symbol...)::Mode = Mode(mode.attributes, Set([mode.identifiers..., identifier...]))

function unidentify(mode::Mode, identifier::Symbol...)::Mode
    identifiers::Set{Symbol} = Set(mode.identifiers)
    foreach(el -> delete!(identifiers, el), identifier)
    return Mode(mode.attributes, identifiers)
end

reidentify(mode::Mode, identifiers::Set{Symbol})::Mode = Mode(mode.attributes, identifiers)

function addattr(mode::Mode, attribute::Pair{Symbol, T})::Mode where {T}
    @assert(attribute.second isa MODE_IDENTIFIERS_TYPE[attribute.first])
    return Mode(Base.ImmutableDict(mode.attributes..., attribute), mode.identifiers)
end

function deleteattr(mode::Mode, identifier::Symbol)::Mode
    intermediate = Dict(mode.attributes)
    delete!(intermediate, identifier)
    return Mode(Base.ImmutableDict(intermediate...), mode.identifiers)
end

function addgroup(mode::Mode, group::String)::Mode
    groups::Vector{String} = [getattr(mode, :groups)..., group]
    return addattr(deleteattr(mode, :groups), :groups => groups)
end

Base.:hash(mode::Mode)::UInt = hash([getattr(mode, identifier) for identifier in mode.identifiers])
Base.:isequal(a::Mode, b::Mode) = a == b

function Base.:(==)(a::Mode, b::Mode)::Bool
    common_identifiers::Set{Symbol} = intersect(a.identifiers, b.identifiers)
    if isempty(common_identifiers) return false end # Cannot determine if they are equal if they aren't identified by the same set of identifiers.
    for identifier in common_identifiers
        if !isequal(getattr(a, identifier), getattr(b, identifier)) return false end
    end
    return true
end

struct FockSpace <: AbstractSpace{Subset{Subset{Mode}}}
    rep::Subset{Subset{Mode}}
    """ Implicit ordering is required as matrices is used for mapping between `FockSpace`. """
    ordering::Dict{Mode, Int64}

    FockSpace(subsets::Subset{Subset{Mode}}) = new(
        subsets,
        Dict(mode => order
             for (order, mode)
             in enumerate([element for subset in rep(subsets) for element in rep(subset)])))

    FockSpace(subset::Subset{Mode}) = FockSpace(convert(Subset{Subset{Mode}}, subset))
end

Base.:iterate(fock_space::FockSpace, i...) = iterate(flatten(rep(fock_space)), i...)

Base.:show(io::IO, fock_space::FockSpace) = print(io, string(rep(flatten(rep(fock_space)))))

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

quantize(name::String, index::Integer, point::Point, flavor::Integer)::Mode = (
    # The new mode will take a group of q:$(name).
    Mode([:groups => [groupname(:q, name)], :index => index, :pos => point, :flavor => flavor]))

quantize(name::String, subset::Subset{Point}, count::Integer)::Subset{Mode} = (
    Subset(Set([quantize(name, index, point, flavor) for (index, point) in enumerate(subset) for flavor in 1:count])))

struct FockMap <: Element{SparseMatrixCSC{ComplexF64, Int64}}
    out_space::FockSpace
    in_space::FockSpace
    rep::SparseMatrixCSC{ComplexF64, Int64}

    FockMap(out_space::FockSpace, in_space::FockSpace, rep::SparseMatrixCSC{ComplexF64, Int64}) = new(out_space, in_space, rep)
    FockMap(out_space::FockSpace, in_space::FockSpace, rep::AbstractArray{<:Number}) = new(out_space, in_space, SparseMatrixCSC{ComplexF64, Int64}(rep))

    function FockMap(out_space::FockSpace, in_space::FockSpace, mapping::Dict{Tuple{Mode, Mode}, ComplexF64})::FockMap
        rep::SparseMatrixCSC{ComplexF64, Int64} = spzeros(dimension(out_space), dimension(in_space))
        for ((out_mode::Mode, in_mode::Mode), value::ComplexF64) in mapping
            rep[out_space.ordering[out_mode], in_space.ordering[in_mode]] = value
        end
        return new(out_space, in_space, rep)
    end
end

Base.:convert(::Type{SparseMatrixCSC{ComplexF64, Int64}}, source::FockMap) = source.rep

function columns(fock_map::FockMap, subset::Subset{Mode})::FockMap
    in_space::FockSpace = FockSpace(subset)
    matrix::SparseMatrixCSC{ComplexF64, Int64} = spzeros(dimension(fock_map.out_space), dimension(in_space))
    for mode in rep(subset)
        matrix[:, in_space.ordering[mode]] = rep(fock_map)[:, fock_map.out_space.ordering[mode]]
    end
    return FockMap(in_space, fock_map.out_space, matrix)
end

function permute(source::FockMap, out_space::FockSpace=source.out_space, in_space::FockSpace=source.in_space)::FockMap
    row_rule::Vector{Int64} = ordering_rule(source.out_space, out_space)
    col_rule::Vector{Int64} = ordering_rule(source.in_space, in_space)
    return FockMap(out_space, in_space, SparseArrays.permute(rep(source), row_rule, col_rule))
end

function Base.:+(a::FockMap, b::FockMap)::FockMap
    @assert(a.in_space == b.in_space && a.out_space == b.out_space)
    return FockMap(a.out_space, a.in_space, rep(a) + rep(permute(b, a.out_space, a.in_space)))
end

function Base.:-(a::FockMap, b::FockMap)::FockMap
    @assert(a.in_space == b.in_space && a.out_space == b.out_space)
    return FockMap(a.out_space, a.in_space, rep(a) - rep(permute(b, a.out_space, a.in_space)))
end

transpose(source::FockMap)::FockMap = FockMap(source.in_space, source.out_space, rep(source)')

dagger(source::FockMap)::FockMap = FockMap(source.in_space, source.out_space, conj(rep(source)'))

function eigvalsh(fock_map::FockMap)::Vector{Pair{Mode, Float64}}
    @assert(fock_map.in_space == fock_map.out_space)
    evs::Vector{Number} = eigvals(Hermitian(Matrix(rep(fock_map))))
    return [Mode([:groups => [groupname(:t, "eigh")], :index => index, :flavor => 1]) => ev for (index, ev) in enumerate(evs)]
end

export Mode, FockSpace, FockMap
export groupname, hasattr, getattr, identify, unidentify, reidentify, addattr, deleteattr, addgroup
export dimension, ordered_modes, ordering_rule, quantize, columns, transpose, dagger, eigvalsh

end
