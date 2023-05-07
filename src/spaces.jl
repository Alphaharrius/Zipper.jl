module Spaces

using LinearAlgebra, OrderedCollections

export Element, AbstractSpace, AffineSpace, RealSpace, MomentumSpace, AbstractSubset, Point, Subset
export rep, euclidean, basis, dimension, spaceof, rpos, pos, lineartransform, fourier_coef, distance, flatten, members

"""
Simply means an distinct type of object that can be represented by a concrete type `R`.
"""
abstract type Element{R <: Any} end

""" Restrict all print out of elements to just printing the type. """
Base.:show(io::IO, element::Element) = print(io, string(typeof(element)))

""" Reflexive relation. """
Base.:convert(::Type{T}, source::T) where {T <: Element} = source
""" The set representation of an `Element` is a set containing itself. """
Base.:convert(::Type{OrderedSet{T}}, source::T) where {T <: Element} = OrderedSet{T}([source])

"""
    rep(source::Element{R}) where {R <: Any}

Retrieve the representation of the `source`, for each subtype of `Element{R}` are ensured to have a representation of type `R`.
This method will throw error if the subtype `T <: Element{R}` does not overload `Base.convert` to convert `T` to `R`.
"""
rep(source::Element{R}) where {R <: Any} = convert(R, source)

"""
Refers to an object with `dimension` spanned by some other objects like `Vector`.
"""
abstract type AbstractSpace{T} <: Element{T} end

"""
    dimension(space::T)::Integer where {T <: AbstractSpace}

The dimension of a `AbstractSpace` is the number of objects that spans this space.

### Input
- `space` The source space.

### Output
The dimension of the `space`, since this method is arbitary to all types `<: AbstractSpace`, the top root returns `0` as `AbstractSpace` does not refer to any
concrete space.
"""
dimension(space::T) where {T <: AbstractSpace} = 0

"""
    hassamespan(a::T, b::F)::Bool where {T <: AbstractSpace, F <: AbstractSpace}

Check whether the two space is spanned by the same set of elements, regardless of transformations.

### Input
- `a` & `b` The two spaces to be checked.

### Output
`true` if the two spaces `a` & `b` have same span; `false` otherwise. This method can have different implementation varies by the type `T` & `F`.
"""
hassamespan(a::T, b::F) where {T <: AbstractSpace, F <: AbstractSpace} = false

"""
This is a less general subtype of `AbstractSpace` that is spanned by `Vector` in a continuous manner.
"""
abstract type AffineSpace <: AbstractSpace{Matrix{Float64}} end

"""
    hassamespan(a::T, b::F) where {T <: AffineSpace, F <: AffineSpace}

Check between two `<: AffineSpace`, returns `false` as `T` and `F` are different subtype of `AffineSpace`.
"""
hassamespan(a::T, b::F) where {T <: AffineSpace, F <: AffineSpace} = false

"""
    hassamespan(a::T, b::T) where {T <: AffineSpace}

Check between two `<: AffineSpace`, returns `true` if they have same dimension.
"""
hassamespan(a::T, b::T) where {T <: AffineSpace} = dimension(a::T) == dimension(b::T)

Base.:(==)(a::AffineSpace, b::AffineSpace)::Bool = typeof(a) == typeof(b) && isapprox(basis(a), basis(b))

Base.:convert(::Type{Matrix{Float64}}, source::AffineSpace) = source.rep

"""
    euclidean(S::Type{<: AffineSpace}, n::Int64)

Create a Euclidean space (with basis vectors of a identity matrix) of type `S <: AffineSpace` with dimension `n`.
"""
euclidean(S::Type{<: AffineSpace}, n::Int64) = S(Matrix{Float64}(I(n)))

"""
    basis(space::AffineSpace)::Matrix{Float64}

Get the basis vectors in form of `Matrix{Float64}` from the `AffineSpace` `space`.
"""
basis(space::AffineSpace)::Matrix{Float64} = rep(space)

"""
    dimension(space::T)::Integer where {T <: AffineSpace}

Get the dimension of the `AffineSpace`, the dimension will be the row / column length of the representation `Matrix`.
"""
dimension(space::T) where {T <: AffineSpace} = size(basis(space), 1)

Base.:hash(space::AffineSpace) = hash(map(v -> Rational{Int64}(round(v * 1000000007)) // 1000000007, basis(space)))

"""
Subset of element type `T <: AbstractSubset` defined within a full set.
"""
abstract type AbstractSubset{T} <: Element{OrderedSet{T}} end

"""
    spaceof(subset::AbstractSubset)

Get the space of the parameter `subset`. If `subset isa Point`, then the output will be its parent space where its position is defined; if `subset isa Subset`, then
the output will be the common space of all elements within the `subset`.
"""
spaceof(subset::AbstractSubset) = subset.space

"""
    RealSpace(rep::Matrix{Float64})

A real space concrete type of `AffineSpace`, this is present in tandem with `MomentumSpace` to distinguish between real space and momentum space.

### Input
- `rep` The representation of the space which is the basis vectors.
"""
struct RealSpace <: AffineSpace
    rep::Matrix{Float64}
end

"""
    MomentumSpace(rep::Matrix{Float64})

A momentum space concrete type of `AffineSpace`, this is present in tandem with `RealSpace` to distinguish between real space and momentum space.

### Input
- `rep` The representation of the space which is the basis vectors.
"""
struct MomentumSpace <: AffineSpace
    rep::Matrix{Float64}
end

Base.:convert(::Type{RealSpace}, source::Matrix{Float64})::RealSpace = RealSpace(source)
Base.:convert(::Type{MomentumSpace}, source::Matrix{Float64})::MomentumSpace = MomentumSpace(source)

""" Performing conversion from `RealSpace` to `MomentumSpace` using formula 𝑅ₖ ≝ 2π⋅(𝑅ᵣ⁻¹)ᵀ."""
Base.:convert(::Type{MomentumSpace}, source::RealSpace)::MomentumSpace = MomentumSpace(2.0 * π * Matrix(transpose(inv(basis(source)))))
""" Performing conversion from `MomentumSpace` to `RealSpace`."""
Base.:convert(::Type{RealSpace}, source::MomentumSpace)::RealSpace = RealSpace(Matrix(transpose(inv(basis(source) / (2.0 * π)))))

"""
    Point(pos::Vector{Float64}, space::AbstractSpace)

A point in an `AffineSpace`.

### Input
- `pos` The vector representation of the point in the spacial unit of the parent `AffineSpace`.
- `space` The parent `AffineSpace`.
"""
struct Point <: AbstractSubset{Point}
    pos::Vector{Float64}
    space::AbstractSpace
end

Base.:show(io::IO, point::Point) = print(io, string("$([trunc(v, digits=5) for v in pos(point)])"))

"""
    rpos(point::Point, denominator::Int64)

Generates a `Vector{Rational{Int64}}` that round the position of the `Vector` in a fixed precision, this method can be used for hashing `Point`.

### Input
- `point`       The source point.
- `denominator` The denominator of the `Rational{Int64}` representation for each element of the position, defaults to `1000000007` which round all error with order
                below or equals `10e-11`.

### Output
A `Vector{Rational{Int64}}` which stores the rounded elements.
"""
rpos(point::Point, denominator::Int64 = 10000000)::Vector{Rational{Int64}} = [Rational{Int64}(round(el * denominator)) // denominator for el in pos(point)]

"""
    dimension(space::T)::Integer where {T <: AffineSpace}

Get the dimension of the `Point`, which is the dimension of its vector representation.
"""
dimension(point::Point)::Integer = length(pos(point))

# rpos is used to to equate the points to round off slight differences.
Base.:(==)(a::Point, b::Point)::Bool = spaceof(a) == spaceof(b) && rpos(a) == rpos(b)
LinearAlgebra.:norm(point::Point) = norm(pos(point))

LinearAlgebra.:dot(a::Point, b::Point)::Float64 = dot(collect(Float64, pos(a)), collect(Float64, pos(b)))

# Since it is not possible to hash the AffineSpace with Matrix{Float64}, we will hash only the Vector representation and leave to the Base.isequal for identification.
Base.:hash(point::Point)::UInt = hash(rpos(point))
Base.:isequal(a::Point, b::Point)::Bool = a == b

"""
    pos(point::Point)::Vector{Float64}

Get the vector representation of a `Point`.
"""
pos(point::Point)::Vector{Float64} = point.pos

"""
    lineartransform(new_space::AffineSpace, point::Point)::Point

Perform linear transformation on the point from the original space `spaceof(point)` to `newspace`, the vector representation of the `point` will be transformed.
"""
lineartransform(newspace::AffineSpace, point::Point)::Point = Point(inv(basis(newspace)) * basis(spaceof(point)) * pos(point), newspace)

"""
    euclidean(point::Point)::Point

Perform linear transformation on the point from the original space `spaceof(point)` to Euclidean space.
"""
euclidean(point::Point)::Point = Point(basis(spaceof(point)) * pos(point), euclidean(typeof(spaceof(point)), dimension(point)))

function Base.:+(a::Point, b::Point)::Point
    @assert(typeof(spaceof(a)) == typeof(spaceof(b)))
    return Point(pos(a) + pos(b), spaceof(a))
end

function Base.:-(a::Point, b::Point)::Point
    @assert(typeof(spaceof(a)) == typeof(spaceof(b)))
    return Point(pos(a) - pos(b), spaceof(a))
end

Base.:*(point::Point, val::T) where {T <: Real} = Point(collect(Float64, pos(point)) * val, spaceof(point))
Base.:/(point::Point, val::T) where {T <: Real} = Point(collect(Float64, pos(point)) / val, spaceof(point))

"""
    Subset(elements::OrderedSet{T})
    Subset(elements::Vector{T})

A concrete type of `AbstractSubset` act as a container of elements of type `<: AbstractSubset` in ordered manner, as a subset of a larger set, this set is assumed
to reside within a given `AbstractSpace`, which all of the elements also belongs to the `AbstractSpace`. To 

### Input
- `elements` The set of elements to be stored, can be type of `OrderedSet{T}` or `Vector{T}` depends on use case, normally `OrderedSet{T}` is suggested if possible
             if it is created before the constructor of this type.
"""
struct Subset{T <: AbstractSubset} <: AbstractSubset{T}
    rep::OrderedSet{T}

    Subset(elements::OrderedSet{T}) where {T <: AbstractSubset} = new{T}(elements)
    Subset(elements::Vector{T}) where {T <: AbstractSubset} = Subset(OrderedSet(elements))
end

"""
    spaceof(subset::Subset)::AbstractSpace

The space of a `Subset` have to be determined by all of it's underlying elements.

### Error
If the underlying elements of `subset` does not belongs to the same space from `spaceof(element)`.
"""
function spaceof(subset::Subset)
    @assert(length(subset) > 0, "Cannot be determine the space of an empty set.")
    space = spaceof(first(subset))
    @assert(all(Iterators.map(el -> spaceof(el) == space, subset)))
    return space
end

# Since it is not possible to hash the AffineSpace with Matrix{Float64}, we will hash only the OrderedSet representation and leave to the Base.isequal for identification.
Base.:hash(subset::Subset)::UInt = hash(rep(subset))
Base.:isequal(a::Subset, b::Subset)::Bool = a == b

# Overloads for iterators.
Base.:iterate(subset::T, i...) where {T <: Subset} = iterate(rep(subset), i...)
Base.:length(subset::T) where {T <: Subset} = length(rep(subset))
Base.:lastindex(subset::T) where {T <: Subset} = length(subset)
Base.:getindex(subset::T, index) where {T <: Subset} = rep(subset)[index]
Base.:getindex(subset::T, inds...) where {T <: Subset} = Subset([rep(subset)...][inds...])

Base.:filter(f, subset::T) where {T <: Subset} = Subset(filter(f, rep(subset)))

"""
    flatten(subset::Subset{T})::Subset

Flatten a higher order subset to order `1`, and preserves the element orderings of each order.
"""
function flatten(subset::Subset{T})::Subset where {T <: AbstractSubset}
    if T <: Subset
        return union([flatten(element) for element in rep(subset)]...)
    end
    return subset
end

"""
    members(subset::Subset)::Tuple

Convert the `Subset` into a `Tuple` with all its elements for easier accessing.

### Examples
- `e₀, e₁ = members`
"""
members(subset::Subset)::Tuple = (rep(subset)...,)

Base.:(==)(a::Subset, b::Subset)::Bool = spaceof(a) == spaceof(b) && Set(rep(a)) == Set(rep(b))

Base.:convert(::Type{OrderedSet{T}}, source::Subset{T}) where {T <: AbstractSubset} = source.rep
""" Reflexive relation. """
Base.:convert(::Type{Subset}, source::Subset{T}) where {T <: AbstractSubset} = source
Base.:convert(::Type{Subset{T}}, source::Subset{T}) where {T <: AbstractSubset} = source
""" Element-wise conversions. """
Base.:convert(::Type{Subset{A}}, source::Subset{B}) where {A <: AbstractSubset, B <: AbstractSubset} = (
    Subset(OrderedSet{A}([convert(A, el) for el in rep(source)])))
""" Convert to the generalized type of `AbstractSubset`. """
Base.:convert(::Type{Subset}, source::T) where {T <: AbstractSubset} = Subset(rep(source))
Base.:convert(::Type{Subset{T}}, source::T) where {T <: AbstractSubset} = convert(Subset, source)

function Base.:union(input::Subset...)::Subset
    @assert(length(OrderedSet{AbstractSpace}([spaceof(subset) for subset in input])) == 1)
    return Subset(union([rep(subset) for subset in input]...))
end

Base.:union(input::T...) where {T <: AbstractSubset} = union((convert(Subset{T}, element) for element in input)...)

function Base.:intersect(input::Subset...)::Subset
    @assert(length(OrderedSet{AbstractSpace}([spaceof(subset) for subset in input])) == 1)
    return Subset(intersect([rep(subset) for subset in input]...))
end

Base.:intersect(input::T...) where {T <: AbstractSubset} = intersect((convert(Subset{T}, element) for element in input)...)

end
