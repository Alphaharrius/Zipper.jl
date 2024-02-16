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
Base.:show(io::IO, ::Type{Offset}) = print(io, "Offset")
Base.:show(io::IO, ::Type{Momentum}) = print(io, "Momentum")

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
Base.:(==)(a::Offset, b::Offset)::Bool = rpos(a) == rpos(b)
# Reciprocal space points are equivalent only if they reside in the same space.
Base.:(==)(a::Momentum, b::Momentum)::Bool = getspace(a) == getspace(b) && rpos(a) == rpos(b)

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

LinearAlgebra.:det(space::AffineSpace)::Real = getbasis(space)|>det

"""
    euclidean(point::Point)::Point

Perform linear transformation on the point from the original space `getspace(point)` to Euclidean space.
"""
Zipper.:euclidean(point::Point)::Point = Point(getbasis(getspace(point)) * vec(point), euclidean(typeof(getspace(point)), dimension(point)))

Base.:-(point::Point)::Point = Point(-vec(point), getspace(point))

function Base.:+(a::Point, b::Point)::Point
    @assert(typeof(getspace(a)) == typeof(getspace(b)))
    localb::Point = (a|>getspace) * b
    return Point(vec(a) + vec(localb), getspace(a))
end

Base.:-(a::Point, b::Point)::Point = a + (-b)

""" Performing linear transform of a `point` to `space`. """
Base.:*(space::T, point::Point{T}) where {T <: AffineSpace} = lineartransform(space, point)

Base.:*(point::Point, val::T) where {T <: Real} = Point(collect(Float64, vec(point)) * val, getspace(point))
Base.:/(point::Point, val::T) where {T <: Real} = Point(collect(Float64, vec(point)) / val, getspace(point))
Base.:*(val::T, point::Point) where {T <: Real} = point * val

# ========================================================
# Subset essentials
Base.:convert(::Type{Vector}, set::Subset) = set.elements
Base.:eltype(::Subset{T}) where T = T
Base.:hash(set::Subset) = hash(set.orderings|>keys)
# ========================================================

# ===============================================================================
# Subset logical methods
""" Equality like a mathematical set ignores orderings. """
Base.:(==)(a::Subset, b::Subset)::Bool = keys(a.orderings) == keys(b.orderings)
Base.:in(v, set::Subset)::Bool = haskey(set.orderings, v)

function Base.:union(set::Subset, iters...)
    orderings = copy(set.orderings)
    headelements = copy(set|>rep)
    headlength = headelements|>length
    others = iters|>Iterators.flatten|>enumerate
    for (n, v) in others
        if haskey(orderings, v)
            continue
        end
        orderings[v] = headlength + n
        push!(headelements, v)
    end
    sorted = sort(orderings|>collect, by=last)
    headorderings = Dict(v=>n for (n, (v, _)) in enumerate(sorted))
    return Subset{set|>eltype}(headelements, headorderings)
end

"""
    subsetunion(iters)

Treating every iterator in `iters` as a mathematical set, 
perform a union operation on all of the sets and return the 
result as a `Subset`.
"""
function subsetunion(iters)
    orderings = Dict()
    elements = []
    items = iters|>Iterators.flatten|>enumerate
    for (n, v) in items
        if haskey(orderings, v)
            continue
        end
        orderings[v] = n
        push!(elements, v)
    end
    sorted = sort(orderings|>collect, by=last)
    orderings = Dict(v=>n for (n, (v, _)) in enumerate(sorted))
    return Subset{elements|>eltype}(elements, orderings)
end
export subsetunion

function Base.:intersect(set::Subset, iters...)
    orderings = copy(set.orderings)
    others = Set(iters|>Iterators.flatten)
    for (v, _) in set.orderings
        v ∈ others || delete!(orderings, v)
    end
    sorted = sort(orderings|>collect, by=last)
    headelements = [v for (v, _) in sorted]
    headorderings = Dict(v=>n for (n, v) in enumerate(headelements))
    return Subset{set|>eltype}(headelements, headorderings)
end

function Base.:setdiff(set::Subset, iters...)
    orderings = copy(set.orderings)
    others = iters|>Iterators.flatten
    for v in others
        if !haskey(orderings, v)
            continue
        end
        delete!(orderings, v)
    end
    sorted = sort(orderings|>collect, by=last)
    headelements = [v for (v, _) in sorted]
    headorderings = Dict(v=>n for (n, v) in enumerate(headelements))
    return Subset{set|>eltype}(headelements, headorderings)
end

Base.:issubset(a::Subset, b::Subset)::Bool = issubset(a.orderings|>keys, b.orderings|>keys)
# ===============================================================================

# ======================================================================================
# Subset iterators methods
Base.:length(set::Subset) = length(set.elements)
Base.:iterate(set::Subset, i...) = iterate(set|>rep, i...)
Base.:lastindex(set::Subset) = length(set)
Base.:getindex(set::Subset, index) = rep(set)[index]
Base.:getindex(set::Subset, range::UnitRange) = Subset(rep(set)[range])

Base.:filter(f, set::Subset) = Subset(Iterators.filter(f, set|>rep))

Iterators.:flatten(set::Subset{<:Subset}) = Subset(v for inner in set for v in inner)
# ======================================================================================

# =====================================================
# Subset indexing
Base.:getindex(set::Subset{T}, v::T) where T = set.orderings[v]
# =====================================================

# ========================================================================================================
# Subset arithmetics
Base.:+(a::Subset, b::Subset) = union(a, b)
Base.:-(a::Subset, b::Subset) = setdiff(a, b)

""" Performs element-wise additions. """
Base.broadcasted(::typeof(+), left::Subset, right::Subset) = Subset(a+b for (a, b) in zip(left, right))
""" Performs element-wise subtractions. """
Base.broadcasted(::typeof(-), left::Subset, right::Subset) = Subset(a-b for (a, b) in zip(left, right))

""" Add `v` to each element in the `Subset`. """
Base.broadcasted(::typeof(+), set::Subset, v) = Subset(el+v for el in set)
""" Subtract `v` from each element in the `Subset`. """
Base.broadcasted(::typeof(-), set::Subset, v) = set .+ (-v)
# ========================================================================================================

# ========================================================================================
# Subset display
Base.:show(io::IO, set::Subset) = print(io, string("$(typeof(set))(len=$(set|>length))"))
# ========================================================================================

"""
    getspace(set::SortSet)

With the assumption that all elements in the `SortSet` resides in the same space, 
this method returns the space of the first element in the `SortSet`.
"""
Zipper.:getspace(set::Subset) = set|>first|>getspace

"""
    members(subset::Subset)::Tuple

Convert the `Subset` into a `Tuple` with all its elements for easier accessing.

### Examples
- `e₀, e₁ = members`
"""
members(set::Subset)::Tuple = (rep(set)...,)
export members

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
