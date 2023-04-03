if !isdefined(Main, :Spaces) include("spaces.jl") end

module Quantum

using ..Spaces

struct Mode <: AbstractSubset{Mode}
    base::Subset
    fermion_count::Int64
    space::AbstractSpace

    Mode(base::Subset, fermion_count::Int64) = new(base, fermion_count, space_of(base))
end

base(mode::Mode)::Subset = mode.base

struct FockSpace <: AbstractSpace{Subset{Subset{Mode}}}
    rep::Subset{Subset{Mode}}
    """ Implicit ordering is required as matrices is used for mapping between `FockSpace`. """
    ordering::Dict{Mode, Int64}

    FockSpace(subsets::Subset{Subset{Mode}}) = new(
        subsets,
        Dict(mode => order
             for (order, mode)
             in enumerate([element for subset in representation(subsets) for element in representation(subset)])))

    FockSpace(subset::Subset{Mode}) = FockSpace(convert(Subset{Subset{Mode}}, subset))
end

ordering(fock_space::FockSpace)::Dict{Mode, Int64} = fock_space.ordering

Base.:(==)(a::FockSpace, b::FockSpace)::Bool = representation(a) == representation(b)

Base.:convert(::Type{Subset}, source::FockSpace) = source.rep
Base.:convert(::Type{Subset{Subset}}, source::FockSpace) = convert(Subset, source)
Base.:convert(::Type{Subset{Subset{Mode}}}, source::FockSpace) = convert(Subset, source)

Base.:convert(::Type{FockSpace}, source::Subset{Mode}) = FockSpace(source)

quantize(point::Point, fermion_count)::Mode = Mode(convert(Subset, point), fermion_count)

function quantize(subset::Subset{T}, fermion_count::Int64)::Subset where {T <: AbstractSubset}
    quantized::Vector = [quantize(element, fermion_count) for element in representation(subset)]
    return Subset(Set(quantized), space_of(quantized[1]))
end

export Mode, FockSpace
export ordering, quantize

end