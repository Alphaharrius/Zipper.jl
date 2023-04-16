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

# This design of Mode have performance issue espectially on hashing.
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

Base.:show(io::IO, mode::Mode) = show(io, "$(string(typeof(mode)))($(mode.identifiers))")

function Spaces.space_of(mode::Mode)::AbstractSpace
    if hasattr(mode, :pos) return space_of(getattr(mode, :pos)) end
    if hasattr(mode, :offset) return space_of(getattr(mode, :offset)) end
    # If the mode does not based on any physical position or quantities for associating with any space, then it will simply lives
    # in a 1D euclidean space as the mode and it's siblings can always be parameterized by a scalar.
    return euclidean(RealSpace, 1)
end

hasattr(mode::Mode, identifier::Symbol) = haskey(mode.attributes, identifier)

getattr(mode::Mode, identifier::Symbol) = mode.attributes[identifier]

identifiers(mode::Mode)::Set{Symbol} = mode.identifiers

identify(mode::Mode, identifier::Symbol...)::Mode = Mode(mode.attributes, Set([mode.identifiers..., identifier...]))

function unidentify(mode::Mode, identifier::Symbol...)::Mode
    identifiers::Set{Symbol} = Set(mode.identifiers)
    foreach(el -> delete!(identifiers, el), identifier)
    return Mode(mode.attributes, identifiers)
end

reidentify(mode::Mode, identifiers::Set{Symbol})::Mode = Mode(mode.attributes, identifiers)

function setattr(mode::Mode, attribute::Pair{Symbol, T})::Mode where {T}
    @assert(attribute.second isa MODE_IDENTIFIERS_TYPE[attribute.first])
    intermediate::Mode = hasattr(mode, attribute.first) ? removeattr(mode, attribute.first) : mode
    return Mode(Base.ImmutableDict(intermediate.attributes..., attribute), intermediate.identifiers)
end

function removeattr(mode::Mode, identifier::Symbol)::Mode
    intermediate = Dict(mode.attributes)
    delete!(intermediate, identifier)
    return Mode(Base.ImmutableDict(intermediate...), mode.identifiers)
end

function addgroup(mode::Mode, group::String)::Mode
    groups::Vector{String} = [getattr(mode, :groups)..., group]
    return setattr(removeattr(mode, :groups), :groups => groups)
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
    ordering::Base.ImmutableDict{Mode, Int64}

    FockSpace(subsets::Subset{Subset{Mode}}, ordering::Base.ImmutableDict{Mode, Int64}) = new(subsets, ordering)

    FockSpace(subset::Subset{Mode}) = new(
        Subset(Set([subset])),
        Base.ImmutableDict([mode => order for (order, mode) in enumerate(subset)]...))
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
    home::Point = origin(space_of(point))
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

transpose(source::FockMap)::FockMap = FockMap(source.inspace, source.outspace, rep(source)')

dagger(source::FockMap)::FockMap = FockMap(source.inspace, source.outspace, rep(source)') # `'` operator is dagger.

function eigvalsh(fock_map::FockMap, eigenattr::Pair{Symbol, T})::Vector{Pair{Mode, Float64}} where {T}
    @assert(fock_map.inspace == fock_map.outspace)
    evs::Vector{Number} = eigvals(Hermitian(Matrix(rep(fock_map))))
    return [Mode([eigenattr, :groups => [groupname(:t, "eigh")], :index => index, :flavor => 1]) => ev for (index, ev) in enumerate(evs)]
end

function fourier(momentums::Subset{Point}, inmodes::Subset{Mode})::FockMap
    dict::Dict{Tuple{Mode, Mode}, ComplexF64} = Dict()
    # Disable the identification by :offset so that all inmodes collapes to the basis mode. 
    basismodes::Set{Mode} = Set(unidentify(inmode, :offset) for inmode in inmodes)
    # Enable the identification by :offset when adding momentum as :offset for each basis mode.
    outmodes = Iterators.map(tup -> identify(setattr(tup[1], :offset => tup[2]), :offset), Iterators.product(basismodes, momentums))
    for (outmode, inmode) in Iterators.product(outmodes, inmodes)
        momentum::Point = getattr(outmode, :offset)
        euc_momentum::Point = linear_transform(euclidean(MomentumSpace, length(momentum)), momentum)
        inoffset::Point = getattr(inmode, :offset)
        euc_inoffset::Point = linear_transform(euclidean(RealSpace, length(inoffset)), inoffset)
        dict[(outmode, inmode)] = exp(-1im * dot(euc_momentum, euc_inoffset))
        dict[(outmode, inmode)] = 1im
    end
    return FockMap(FockSpace(Subset(Set(outmodes))), FockSpace(inmodes), dict)
end

export Mode, FockSpace, FockMap
export groupname, hasattr, getattr, identify, unidentify, reidentify, setattr, removeattr, addgroup
export dimension, ordered_modes, ordering_rule, quantize, columns, transpose, dagger, eigvalsh, fourier

end
