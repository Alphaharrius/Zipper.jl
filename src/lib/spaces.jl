# Overload the hashing function for vector of points to yield the same value with the same data content.
Base.:hash(elements::Vector{T}) where {T <: Element} = hash([element |> hash for element in elements])

"""
    dimension(space::T)::Integer where {T <: AbstractSpace}

The dimension of a `AbstractSpace` is the number of objects that spans this space.

### Input
- `space` The source space.

### Output
The dimension of the `space`, since this method is arbitary to all types `<: AbstractSpace`, the top root returns `0` as `AbstractSpace` does not refer to any
concrete space.
"""
Zipper.:dimension(space::T) where {T <: AbstractSpace} = error("Dimension resolving is not defined for `$(space |> typeof)`!")

"""
    hassamespan(a::T, b::F)::Bool where {T <: AbstractSpace, F <: AbstractSpace}

Check whether the two space is spanned by the same set of elements, regardless of transformations.

### Input
- `a` & `b` The two spaces to be checked.

### Output
`true` if the two spaces `a` & `b` have same span; `false` otherwise. This method can have different implementation varies by the type `T` & `F`.
"""
Zipper.:hassamespan(a::T, b::F) where {T <: AbstractSpace, F <: AbstractSpace} = false

""" Shorthand for creating a `Point` within the specified space. """
Base.:∈(data, space::AffineSpace)::Point = Point(data |> collect, space)

"""
    hassamespan(a::T, b::F) where {T <: AffineSpace, F <: AffineSpace}

Check between two `<: AffineSpace`, returns `false` as `T` and `F` are different subtype of `AffineSpace`.
"""
Zipper.:hassamespan(a::T, b::F) where {T <: AffineSpace, F <: AffineSpace} = false

"""
    hassamespan(a::T, b::T) where {T <: AffineSpace}

Check between two `<: AffineSpace`, returns `true` if they have same dimension.
"""
Zipper.:hassamespan(a::T, b::T) where {T <: AffineSpace} = dimension(a::T) == dimension(b::T)

Base.:(==)(a::AffineSpace, b::AffineSpace)::Bool = typeof(a) == typeof(b) && isapprox(getbasis(a), getbasis(b))

Base.:convert(::Type{Matrix{Float64}}, source::AffineSpace) = source.rep

"""
    euclidean(S::Type{<: AffineSpace}, n::Int64)

Create a Euclidean space (with basis vectors of a identity matrix) of type `S <: AffineSpace` with dimension `n`.
"""
euclidean(S::Type{<: AffineSpace}, n::Int64) = S(Matrix{Float64}(I(n)))
export euclidean

"""
    getbasis(space::AffineSpace)::Matrix{Float64}

Get the basis vectors in form of `Matrix{Float64}` from the `AffineSpace` `space`.
"""
getbasis(space::AffineSpace)::Matrix{Float64} = rep(space)
export getbasis

"""
    dimension(space::T)::Integer where {T <: AffineSpace}

Get the dimension of the `AffineSpace`, the dimension will be the row / column length of the representation `Matrix`.
"""
Zipper.:dimension(space::T) where {T <: AffineSpace} = size(getbasis(space), 1)

Base.:hash(space::AffineSpace) = hash(map(v -> Rational{Int64}(round(v * 10000000)) // 10000000, getbasis(space)))

""" Shorthand to get the euclidean space of the same span. """
euclidean(space::AffineSpace) = euclidean(space |> typeof, space |> dimension)

getbasisvectors(space::AffineSpace) = ((space |> rep)[:, d] ∈ euclidean(space) for d in axes(space |> rep, 2))
export getbasisvectors

Base.:convert(::Type{RealSpace}, source::Matrix{Float64})::RealSpace = RealSpace(source)
Base.:convert(::Type{MomentumSpace}, source::Matrix{Float64})::MomentumSpace = MomentumSpace(source)

""" Performing conversion from `RealSpace` to `MomentumSpace` using formula 𝑅ₖ ≝ 2π⋅(𝑅ᵣ⁻¹)ᵀ."""
Base.:convert(::Type{MomentumSpace}, source::RealSpace)::MomentumSpace = MomentumSpace(2.0 * π * Matrix(transpose(inv(getbasis(source)))))
""" Performing conversion from `MomentumSpace` to `RealSpace`."""
Base.:convert(::Type{RealSpace}, source::MomentumSpace)::RealSpace = RealSpace(Matrix(transpose(inv(getbasis(source) / (2.0 * π)))))

Base.:show(io::IO, point::Point) = print(io, string("$([trunc(v, digits=5) for v in vec(point)]) ∈ $(point |> getspace |> typeof)"))

"""
    getspace(point::Point)

Retrieve the `AffineSpace` of the `Point`.
"""
Zipper.:getspace(point::Point) = point.space

hashablereal(v::Real, denominator::Integer = 10000000)::Rational = ((v * denominator) |> round |> Integer) // denominator
export hashablereal

hashablecomplex(z::Complex, denominator::Integer = 10000000)::Tuple = (hashablereal(z |> real, denominator), hashablereal(z |> imag, denominator))
export hashablecomplex

hashablenumber(v::Real, denominator::Integer = 10000000) = hashablereal(v, denominator)
hashablenumber(z::Complex, denominator::Integer = 10000000) = hashablecomplex(z, denominator)
export hashablenumber

global spatialhashdenominators::Vector{Integer} = [128, 128, 128] # Default to 128 for each dimension.
global reciprocalhashdenominators::Vector{Integer} = [128, 128, 128]

"""
    spatialsnappingcalibration(positions)

Analyse the given `positions` to determine the spatial denominators for each dimension to be used for hashing.

### Input
- `positions` The set of positions to be analysed that is iterable.
"""
function spatialsnappingcalibration(positions)
    values::Matrix = hcat(map(p -> p |> euclidean |> vec, positions)...)
    foreach(n -> spatialhashdenominators[n] = snappingdenominator(values[n, :]).denominator, axes(values, 1))
    @warn "Updated position hash denominators to $spatialhashdenominators."
end
export spatialsnappingcalibration

function spatialsnappingcalibration(; x::Integer = 128, y::Integer = 128, z::Integer = 128)
    spatialhashdenominators[1] = x
    spatialhashdenominators[2] = y
    spatialhashdenominators[3] = z
    @warn "Updated position hash denominators to $spatialhashdenominators."
end

"""
    reciprocalhashcalibration(crystalsizes::Vector{<:Integer})

Analyse the given `crystalsizes` to determine the reciprocal denominators for each dimension to be used for hashing.

### Input
- `crystalsizes` The maximum crystal sizes.
"""
function reciprocalhashcalibration(crystalsizes::Vector{<:Integer})
    foreach(e -> reciprocalhashdenominators[e |> first] = e |> last, crystalsizes |> enumerate)
    @warn "Updated momentum hash denominators to $reciprocalhashdenominators."
end
export reciprocalhashcalibration

function reciprocalhashcalibration(; x::Integer = 128, y::Integer = 128, z::Integer = 128)
    reciprocalhashdenominators[1] = x
    reciprocalhashdenominators[2] = y
    reciprocalhashdenominators[3] = z
    @warn "Updated momentum hash denominators to $reciprocalhashdenominators."
end

"""
    rpos(position::Offset)
    rpos(momentum::Momentum)

Generates a `Vector{Rational{Int64}}` that round the position of the `Vector` in a fixed precision, this method can be used for hashing `Point`.
Noted that `rpos` cannot be used for computation since the denominators defined globally are not intented to provide a precise representation of the `Point`,
this method is preserved for hashing purposes only.

### Input
- `position`    The source position.
- `momentum`    The source momentum.
- `denominator` The denominator of the `Rational{Int64}` representation for each element of the position, defaults to `10000000` which round all error with order
                below or equals `10e-7`.

### Output
A `Vector{Rational{Int64}}` which stores the rounded elements.
"""
rpos(position::Offset)::Vector{Rational{Integer}} = [hashablereal(r, d) for (r, d) in zip(position |> euclidean |> vec, spatialhashdenominators[1:dimension(position)])]
rpos(momentum::Momentum)::Vector{Rational{Integer}} = [hashablereal(r, d) for (r, d) in zip(momentum |> vec, reciprocalhashdenominators[1:dimension(momentum)])]
export rpos

"""
    dimension(space::T)::Integer where {T <: AffineSpace}

Get the dimension of the `Point`, which is the dimension of its vector representation.
"""
Zipper.:dimension(point::Point)::Integer = point |> vec |> length

# rpos is used to to equate the points to round off slight differences.
function Base.:(==)(a::Point, b::Point)::Bool
    @assert(typeof(a |> getspace) == typeof(b |> getspace))
    return rpos(a) == rpos(b)
end

LinearAlgebra.:norm(point::Point) = point |> vec |> norm
LinearAlgebra.:normalize(point::Point) = Point(point |> vec |> normalize, point |> getspace)
LinearAlgebra.:dot(a::Point, b::Point)::Float64 = dot(collect(Float64, a |> vec), collect(Float64, b |> vec))

# Since it is not possible to hash the AffineSpace with Matrix{Float64}, we will hash only the Vector representation and leave to the Base.isequal for identification.
Base.:hash(point::Point)::UInt = hash(rpos(point))
Base.:isequal(a::Point, b::Point)::Bool = a == b

""" Get the vector representation of a `Point`. """
Base.:vec(point::Point)::Vector{Float64} = point.pos

"""
    lineartransform(new_space::AffineSpace, point::Point)::Point

Perform linear transformation on the point from the original space `getspace(point)` to `newspace`, the vector representation of the `point` will be transformed.
"""
lineartransform(newspace::AffineSpace, point::Point)::Point = Point(inv(getbasis(newspace)) * getbasis(getspace(point)) * vec(point), newspace)
export lineartransform

"""
    orthospace(space::AffineSpace)::AffineSpace

Given a affine space, get it't corresponding orthogonal space of the same scale.
"""
function orthospace(space::AffineSpace)::AffineSpace
    basisvectors::Base.Generator = ((space |> rep)[:, d] for d in axes(space |> rep, 2))
    scalings::Vector = [norm(v) for v in basisvectors]
    return (space |> typeof)(scalings |> diagm)
end
export orthospace

"""
    affinespace(points::Point...)::AffineSpace

Create a `AffineSpace` with the `basisvectors` provided.
"""
function affinespace(basisvectors::Point...)::AffineSpace
    basis::Matrix = hcat((p |> vec for p in basisvectors)...)
    isapprox(basis |> det, 0, atol=1e-10) && error("Cannot create affine space with linearly dependent basis vectors!")
    return (basisvectors[1] |> getspace |> typeof)(basis)
end
export affinespace

"""
    euclidean(point::Point)::Point

Perform linear transformation on the point from the original space `getspace(point)` to Euclidean space.
"""
Zipper.:euclidean(point::Point)::Point = Point(getbasis(getspace(point)) * vec(point), euclidean(typeof(getspace(point)), dimension(point)))

Base.:-(point::Point)::Point = Point(-vec(point), getspace(point))

function Base.:+(a::Point, b::Point)::Point
    @assert(typeof(getspace(a)) == typeof(getspace(b)))
    return Point(vec(a) + vec(b), getspace(a))
end

function Base.:-(a::Point, b::Point)::Point
    @assert(typeof(getspace(a)) == typeof(getspace(b)))
    return Point(vec(a) - vec(b), getspace(a))
end

""" Performing linear transform of a `point` to `space`. """
Base.:*(space::T, point::Point{T}) where {T <: AffineSpace} = lineartransform(space, point)

Base.:*(point::Point, val::T) where {T <: Real} = Point(collect(Float64, vec(point)) * val, getspace(point))
Base.:/(point::Point, val::T) where {T <: Real} = Point(collect(Float64, vec(point)) / val, getspace(point))
Base.:*(val::T, point::Point) where {T <: Real} = point * val

Zipper.:Subset(points::Point...) = Subset(p for p in points)

Base.:show(io::IO, subset::Subset) = print(io, string("$(typeof(subset))(len=$(length(subset)))"))

Base.:in(item, subset::Subset)::Bool = item in (subset |> rep)

""" Allows additions for a subset with another type. """
# TODO: Requires revision.
Base.:+(subset::Subset, val::Element)::Subset = Subset(p + val for p in subset)
Base.:-(subset::Subset, val::Element)::Subset = subset + (-val)

"""
    getspace(subset::Subset)::AbstractSpace

The space of a `Subset` have to be determined by all of it's underlying elements.

### Error
If the underlying elements of `subset` does not belongs to the same space from `getspace(element)`.
"""
function Zipper.:getspace(subset::Subset)
    @assert(length(subset) > 0, "Cannot be determine the space of an empty set.")
    space = getspace(first(subset))
    @assert(all(Iterators.map(el -> getspace(el) == space, subset)))
    return space
end

# Since it is not possible to hash the AffineSpace with Matrix{Float64}, we will hash only the OrderedSet representation and leave to the Base.isequal for identification.
Base.:hash(subset::Subset)::UInt = hash(rep(subset))
Base.:isequal(a::Subset, b::Subset)::Bool = a == b

# Overloads for iterators.
Base.:iterate(subset::T, i...) where {T <: Subset} = iterate(rep(subset), i...)
Base.:length(subset::T) where {T <: Subset} = subset |> rep |> length
Base.:lastindex(subset::T) where {T <: Subset} = length(subset)
Base.:getindex(subset::T, index) where {T <: Subset} = rep(subset)[index]
Base.:getindex(subset::T, range::UnitRange) where {T <: Subset} = Subset([rep(subset)...][range])

Base.:filter(f, subset::T) where {T <: Subset} = Subset(Iterators.filter(f, rep(subset)))
Iterators.:filter(f, subset::T) where {T <: Subset} = Base.filter(f, subset)

"""
    flatten(subset::Subset{T})::Subset

Flatten a higher order subset to order `1`, and preserves the element orderings of each order.
"""
function Iterators.flatten(subset::Subset{T})::Subset where {T}
    if T <: Subset
        return union((flatten(element) for element in rep(subset))...)
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
export members

Base.:(==)(a::Subset, b::Subset)::Bool = Set(rep(a)) == Set(rep(b))

Base.:convert(::Type{OrderedSet{T}}, source::Subset{T}) where {T} = source.rep
""" Reflexive relation. """
Base.:convert(::Type{Subset}, source::Subset{T}) where {T} = source
Base.:convert(::Type{Subset{T}}, source::Subset{T}) where {T} = source
""" Element-wise conversions. """
Base.:convert(::Type{Subset{A}}, source::Subset{B}) where {A, B} = (
    Subset(OrderedSet{A}([convert(A, el) for el in rep(source)])))
""" Convert to the generalized type of `AbstractSubset`. """
Base.:convert(::Type{Subset}, source::T) where {T} = Subset(rep(source))
Base.:convert(::Type{Subset{T}}, source::T) where {T} = convert(Subset, source)

Base.:eltype(subset::Subset) = subset |> rep |> eltype

subsetunion(subsets)::Subset = union(OrderedSet(), (v for subset in subsets for v in subset |> rep)) |> Subset
export subsetunion

function subsetintersect(subsets)::Subset
    intersections::OrderedSet = intersect((subset |> rep for subset in subsets)...)
    if intersections |> isempty
        return Subset{subsets |> first |> first |> typeof}()
    end
    return Subset(intersections)
end
export subsetintersect

Base.:union(subsets::Subset...)::Subset = subsetunion(subsets)
Base.:intersect(subsets::Subset...)::Subset = subsetintersect(subsets)

Base.:+(a::Subset, b::Subset)::Subset = union(a, b)

function Base.:setdiff(a::Subset, b::Subset)::Subset
    diffresult::OrderedSet = setdiff(a |> rep, b |> rep)
    if diffresult |> isempty
        return Subset{a |> first |> typeof}()
    end
    return Subset(diffresult)
end

Base.:-(a::Subset, b::Subset)::Subset = setdiff(a, b)

struct SnappingResult
    forvalues::Vector{Number}
    denominator::Integer
end
export SnappingResult

"""
    snappingdenominator(values; denominatorrange::UnitRange = 2:128, tolerantscalepercent::Number = 1/4)::SnappingResult

Find the smallest denominator that yield a measurement scale `1/denominator` that every value in `values` can be snapped to.

### Input
- `values`               The values to be snapped.
- `denominatorrange`     The range of denominator to be searched.
- `tolerantscalepercent` The tolerance of the snapping, the denominator will be accepted iff all values are at max `measurementscale * tolerantscalepercent` away from
                         the nearest marker defined by `n * measurementscale` where `n` ∈ ℝ.

### Output
A `SnappingResult` that contains the `forvalues` which contains the decimal part of `values`, and the `denominator` which is the denominator used to snap the values.
"""
function snappingdenominator(values; denominatorrange::UnitRange = 2:128, tolerantscalepercent::Number = 1/4)::SnappingResult
    unitvalues::Base.Generator = (abs(v % 1) for v in values)

    function trydenominator(denominator::Integer)::Bool
        measurementscale::Number = 1 / denominator
        markers::LinRange = LinRange(0, 1, denominator + 1)
        tolerance::Number = measurementscale * abs(tolerantscalepercent)
        checker::Function = v -> findall(v -> v, abs(marker - v) < tolerance for marker in markers) |> length == 1 # Make sure that all values can be snapped to a marker.
        return all(checker, unitvalues)
    end

    for denominator in denominatorrange
        if trydenominator(denominator) return SnappingResult([unitvalues...], denominator) end
    end

    error("No denominator found for values $(values)!")
end
export snappingdenominator

"""
    findcomplexdenominator(values; denominatorrange::UnitRange = 2:128, tolerantscalepercent::Number = 1/4)::SnappingResult

Find the smallest denominator that yield a measurement scale `1/denominator` that every complex value in `values` can be snapped to.

### Input
- `values`               The values to be snapped.
- `denominatorrange`     The range of denominator to be searched.
- `tolerantscalepercent` The tolerance of the snapping, the denominator will be accepted iff all values are at max `measurementscale * tolerantscalepercent` away from
                         the nearest marker defined by `n * measurementscale` where `n` ∈ ℝ.

### Output
A `SnappingResult` that contains the `forvalues` which contains the decimal part of each complex real & imaginary part in `values`, and the `denominator`
which is the denominator used to snap the values.
"""
function findcomplexdenominator(values; denominatorrange::UnitRange = 2:128, tolerantscalepercent::Number = 1/4)::SnappingResult
    flattenvalues::Base.Iterators.Flatten = (r for v in values for r in (v |> real, v |> imag)) # Flatten the complex values into real & imaginary parts.
    return snappingdenominator(flattenvalues; denominatorrange=denominatorrange, tolerantscalepercent=tolerantscalepercent)
end
export findcomplexdenominator

function dosnf(matrix::Matrix)::Tuple{Matrix, Matrix, Matrix}
    snf = smith(matrix)
    return snf.S, diagm(snf), snf.T
end
export dosnf
