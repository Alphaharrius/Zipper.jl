"""
Refers to an object with `dimension` spanned by some other objects like `Vector`.
"""
abstract type AbstractSpace{T} <: Element{T} end
export AbstractSpace

"""
This is a less general subtype of `AbstractSpace` that is spanned by `Vector` in a continuous manner.
"""
abstract type AffineSpace <: AbstractSpace{Matrix{Float64}} end
export AffineSpace

"""
    RealSpace(rep::Matrix{Float64})

A real space concrete type of `AffineSpace`, this is present in tandem with `MomentumSpace` to distinguish between real space and momentum space.

### Input
- `rep` The representation of the space which is the basis vectors.
"""
struct RealSpace <: AffineSpace
    rep::Matrix{Float64}
end
export RealSpace

"""
    MomentumSpace(rep::Matrix{Float64})

A momentum space concrete type of `AffineSpace`, this is present in tandem with `RealSpace` to distinguish between real space and momentum space.

### Input
- `rep` The representation of the space which is the basis vectors.
"""
struct MomentumSpace <: AffineSpace
    rep::Matrix{Float64}
end
export MomentumSpace

"""
Subset of element type `T <: AbstractSubset` defined within a full set.
"""
# TODO: Is this type still needed?
abstract type AbstractSubset{T} <: Element{OrderedSet{T}} end
export AbstractSubset

"""
    Point(pos::Vector{Float64}, space::AbstractSpace)

A point in an `AffineSpace`.

### Input
- `pos` The vector representation of the point in the spacial unit of the parent `AffineSpace`.
- `space` The parent `AffineSpace`.
"""
struct Point{T <: AffineSpace} <: AbstractSubset{Point}
    pos::Vector{Float64}
    space::T

    Point(pos, space::T) where {T <: AffineSpace} = new{space |> typeof}([pos...], space)
end
export Point

""" Alias of a real space point. """
Offset = Point{RealSpace}
""" Alias of a momentum space point. """
Momentum = Point{MomentumSpace}
export Offset, Momentum

"""
A concrete type of `AbstractSubset` act as a container of elements of type `<: AbstractSubset` in ordered manner, 
as a subset of a larger set, this set is assumed to reside within a given `AbstractSpace`, which all of the elements 
also belongs to the `AbstractSpace`.
"""
struct Subset{T} <: Element{Vector}
    elements::Vector{T}
    orderings::Dict{T, Integer}
end

""" Create an empty `Subset` with `eltype` of `T`. """
Subset{T}() where T = Subset{T}(tuple(), Dict())

""" Create a `Subset` from an iterable `iter` of `T`. """
function Subset(iter)
    T::Type = iter|>first|>typeof
    orderings::Dict{T, Integer} = Dict()
    for (n, v) in iter|>enumerate
        haskey(orderings, v) || (orderings[v] = n)
    end
    sorted = sort(orderings|>collect, by=last)
    # If there are repeative elements in iter, we need to reassign the orderings.
    orderings = Dict(v=>n for (n, (v, _)) in sorted|>enumerate)
    elements::Vector{T} = [key for (key, _) in sorted]
    return Subset{T}(elements, orderings)
end

""" Allow using shorthand of `Subset(v0, v1, ...)`. """
Subset(elements::Element...) = Subset(elements)

export Subset

""" Interface to get the dimension of an `Element`. """
dimension(::Element)::Integer = notimplemented()
export dimension

""" Interface to check if two `Element` objects have the same span. """
hassamespan(::Element, ::Element)::Bool = notimplemented()
export hassamespan

""" Interface to retrieve the space of an `Element`. """
getspace(::Element)::AffineSpace = notimplemented()
export getspace

""" Interface to retrieve the position of an `Element`. """
getpos(::Element) = notimplemented()
export getpos
