if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :DataStructures) include("datastructures.jl") end

module Quantum

using LinearAlgebra, SparseArrays, Arpack
using ..Spaces, ..DataStructures

struct Mode <: AbstractSubset{Mode}
    name::String
    indexes::Tuple{Vararg{Int64}}
    base::Subset
    flavor::Int64
    space::AbstractSpace

    Mode(name::String, indexes::Tuple{Vararg{Int64}}, base::Subset, flavor::Int64) = new(name, indexes, base, flavor, space_of(base))
end

base_of(mode::Mode)::Subset = mode.base

Base.:hash(mode::Mode)::UInt = hash((mode.name, mode.indexes, mode.flavor))
Base.:isequal(a::Mode, b::Mode) = a == b

Base.:show(io::IO, mode::Mode) = print(io, "mode:$(mode.name)$(string(mode.indexes))")

Base.:(==)(a::Mode, b::Mode)::Bool = a.name == b.name && a.indexes == b.indexes && a.base == b.base && a.flavor == b.flavor

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

quantize(name::String, indexes::Tuple{Vararg{Int64}}, point::Point, flavor::Integer)::Mode = (
    Mode(name, indexes, convert(Subset, point), flavor))

quantize(name::String, subset::Subset, count::Int64)::Subset = (
    quantize(name, tuple(), subset, count))

function quantize(name::String, indexes::Tuple{Vararg{Int64}}, subset::Subset{T}, count::Integer)::Subset where {T <: AbstractSubset}
    quantized::Vector = [quantize(name, (indexes..., index), element, flavor) for (index, element) in enumerate(rep(subset)) for flavor in 1:count]
    return Subset(Set(quantized))
end

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

function eigvalsh(fock_map::FockMap)::Vector{Composite{Mode, Float64}}
    @assert(fock_map.in_space == fock_map.out_space)
    evals::Vector{Number} = eigvals(Hermitian(Matrix(rep(fock_map))))
    base::Subset = Subset(Set([data for mode in fock_map.in_space for data in base_of(mode)]))
    return [Composite(Mode("eigmode", (index,), base, 1), real(eval)) for (index, eval) in enumerate(evals)]
end

export Mode, FockSpace, FockMap
export dimension, ordered_modes, ordering_rule, quantize, columns, transpose, dagger, eigvalsh

end
