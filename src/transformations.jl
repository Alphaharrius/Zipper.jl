module Transformations

using LinearAlgebra, SmithNormalForm, OrderedCollections, Combinatorics, Statistics
using ..Spaces, ..Geometries

abstract type Transformation{T} <: Element{T} end
export Transformation

struct BasisFunction <: Element{Vector{Complex}}
    rep::Vector{Complex}
    dimension::Integer
    rank::Integer
end
export BasisFunction

Base.:(==)(a::BasisFunction, b::BasisFunction) = a.rank == b.rank && (a |> dimension) == (b |> dimension) && isapprox(a |> rep, b |> rep)

Base.:hash(basisfunction::BasisFunction)::UInt = (map(hashablecomplex, basisfunction |> rep), basisfunction |> dimension, basisfunction.rank) |> hash
Base.:isequal(a::BasisFunction, b::BasisFunction)::Bool = a == b

swave::BasisFunction = BasisFunction([1], 0, 0)
export swave

LinearAlgebra.normalize(basis::BasisFunction)::BasisFunction = BasisFunction(basis.rep |> normalize, basis.dimension, basis.rank)
end

function Base.:show(io::IO, basisfunction::BasisFunction)
    if basisfunction.rank == 0 # Handles the case of rank 0 basis function.
        print(io, string("$(typeof(basisfunction))(rank=$(basisfunction.rank), $(basisfunction |> rep |> first))"))
        return
    end

    invindexmap::Dict{Integer, Tuple} = Dict(i => t for (t, i) in entryindexmap(basisfunction |> dimension, basisfunction.rank))
    invindextable::Dict{Integer, String} = Dict(i => s for (s, i) in INDEXTABLE)
    function generatesymbol(info::Tuple)::String
        if isapprox(info |> last |> abs, 0, atol=1e-10) return "" end
        coords::Tuple = invindexmap[info |> first]
        return "($(info |> last))*" * reduce(*, Iterators.map(c -> invindextable[c], coords))
    end
    expression::String = join(filter(v -> v != "", map(generatesymbol, basisfunction |> rep |> enumerate)), " + ")
    print(io, string("$(typeof(basisfunction))(rank=$(basisfunction.rank), $(expression))"))
end

Base.:convert(::Type{Vector{Complex}}, source::BasisFunction) = source.rep

Spaces.dimension(basisfunction::BasisFunction)::Integer = basisfunction.dimension

INDEXTABLE::Dict{String, Integer} = Dict("x" => 1, "y" => 2, "z" => 3) # The mappings for axis (x y z) in basis functions to tensor indicies.

resolveentrycoords(expression::Symbol)::Vector{Integer} = map(v -> INDEXTABLE[v], split(expression |> string, ""))

entryindexmap(dimension::Integer, rank::Integer)::Dict{Tuple, Integer} = Dict(v => n for (n, v) in Iterators.map(reverse, Iterators.product(repeat([1:dimension], rank)...)) |> enumerate)

function BasisFunction(expressions::Pair{Symbol, <:Number}...; dimension::Integer)::BasisFunction
    ranks::Tuple = map(expression -> expression |> first |> string |> length, expressions)
    maxrank::Integer = max(ranks...)
    if min(ranks...) != maxrank
        error("Function with mixed order elements is not a valid basis function!")
    end

    indexmap::Dict{Tuple, Integer} = entryindexmap(dimension, maxrank)
    data::Vector{Complex} = zeros(ComplexF64, indexmap |> length)

    for expression in expressions
        coords::Tuple = Tuple(expression |> first |> resolveentrycoords)
        data[indexmap[coords]] = expression |> last
    end

    return BasisFunction(data, dimension, maxrank)
end

eigenfunctionsignature(eigenfunction::BasisFunction, eigenvalue::Number)::Tuple = eigenfunctionsignature(
    eigenfunction.rank, eigenfunction.dimension, eigenvalue)

eigenfunctionsignature(rank::Integer, dimension::Integer, eigenvalue::Number)::Tuple = (
    rank, dimension, eigenvalue |> ComplexF64 |> hashablecomplex)

function Base.:+(a::BasisFunction, b::BasisFunction)::BasisFunction
    @assert(a.dimension == b.dimension)
    @assert(a.rank == b.rank)
    return BasisFunction((a |> rep) + (b |> rep), a.dimension, a.rank)
end

Base.:iszero(basis::BasisFunction)::Bool = basis |> normalize |> rep |> iszero

struct AffineTransform <: Transformation{Matrix{Float64}}
    localspace::AffineSpace
    shiftvector::Vector{Float64}
    transformmatrix::Matrix{Float64}
    eigenfunctions::Dict{Tuple, BasisFunction}
    eigenvaluehashdenominator::Integer
    antiunitary::Bool
end
export AffineTransform
pointgrouptransform(
    pointgroupmatrix::Matrix;
    dimension::Integer = pointgroupmatrix |> size |> first,
    localspace::RealSpace = euclidean(RealSpace, dimension),
    referencepoint::Position = localspace |> origin,
    antiunitary::Bool = false)::AffineTransform = AffineTransform(pointgroupmatrix, transformationshift(pointgroupmatrix, localspace, referencepoint);
    localspace=localspace, antiunitary=antiunitary)
export pointgrouptransform

recenter(transformation::AffineTransform, center::Position)::AffineTransform = AffineTransform(
    transformation.transformmatrix, transformationshift(transformation.transformmatrix, transformation.localspace, center);
    localspace=transformation |> getspace, antiunitary=transformation.antiunitary)
export recenter

recenter(center::Position) = transformation::PointGroupTransformation -> recenter(transformation, center)

translation(
    shiftvector::Vector;
    dimension::Integer = shiftvector |> length,
    localspace::RealSpace = euclidean(RealSpace, dimension),
    antiunitary::Bool = false)::AffineTransform = AffineTransform(
    Matrix{Float64}(I, dimension, dimension), shiftvector; localspace=localspace, antiunitary=antiunitary)
export translation

function affinematrix(transformation::AffineTransform)::Matrix{Float64}
    shiftrow::Vector = vcat(transformation.shiftvector, [1])
    leftcolumns::Matrix = vcat(
        transformation.transformmatrix, zeros(Float64, transformation.transformmatrix |> size |> first, 1) |> transpose)
    return hcat(leftcolumns, shiftrow)
end
export affinematrix

function transformationshift(transformmatrix::Matrix, localspace::AffineSpace, reference::Point)::Vector
    localreference::Point = lineartransform(localspace, reference)
    return (localreference |> pos) - (transformmatrix * (localreference |> pos))
end

Spaces.:getspace(transformation::AffineTransform)::AffineSpace = transformation.localspace

Base.:convert(::Type{Matrix{Float64}}, source::AffineTransform) = source |> affinematrix

function Base.:(==)(a::AffineTransform, b::AffineTransform)::Bool
    localb::AffineTransform = (a |> getspace) * b
    return isapprox(a |> rep, localb |> rep) && a.antiunitary == localb.antiunitary
end

function Base.:*(space::RealSpace, transformation::AffineTransform)::AffineTransform
    if space |> dimension != transformation |> dimension
        error("Dimension mismatch!")
    end
    relativebasis::Matrix = (space |> basis |> inv) * (transformation |> getspace |> basis)
    transformmatrix::Matrix = relativebasis * (transformation.transformmatrix) * (relativebasis |> inv)
    shiftvector::Vector = lineartransform(space, transformation.localspace & transformation.shiftvector) |> pos
    return AffineTransform(
        transformmatrix, shiftvector;
        localspace=space, antiunitary=transformation.antiunitary)
end

function Base.:*(a::AffineTransform, b::AffineTransform)::AffineTransform
    dim::Integer = a |> dimension
    if dim != b |> dimension
        error("Dimension mismatch!")
    end
    localb::AffineTransform = (a |> getspace) * b
    affinematrix::Matrix = (a |> rep) * (localb |> rep)
    transformmatrix::Matrix = affinematrix[1:dim, 1:dim]
    shiftvector::Vector = affinematrix[1:dim, end]
    antiunitary::Bool = a.antiunitary ⊻ localb.antiunitary
    return AffineTransform(
        transformmatrix, shiftvector;
        localspace=a |> getspace, antiunitary=antiunitary)
end
Base.:*(transformation::AffineTransform, space::RealSpace)::RealSpace = RealSpace((transformation.transformmatrix) * (space |> rep))

function Base.:^(source::AffineTransform, exponent::Number)::AffineTransform
    if exponent == 0
        return AffineTransform(Matrix{Float64}(I, source |> dimension, source |> dimension))
    end
    return reduce(*, repeat([source], exponent))
end

function Base.:*(transformation::AffineTransform, region::Subset{Position})::Subset{Position}
    nativetransformation::AffineTransform = (region |> getspace) * transformation
    transformed::Base.Generator = (
        Point(((nativetransformation |> rep) * vcat(point |> pos, [1]))[1:end - 1], point |> getspace)
        for point in region)
    return Subset(transformed)
end

function Base.:*(transformation::AffineTransform, zone::Subset{Momentum})::Subset{Momentum}
    kspace::MomentumSpace = zone |> getspace
    realspace::RealSpace = convert(RealSpace, kspace)
    nativetransformation::AffineTransform = realspace * transformation
    kspacerep::Matrix = Matrix(nativetransformation.transformmatrix |> transpose |> inv)
    return Subset(Point(kspacerep * (k |> pos), kspace) for k in zone)
end

Base.:*(transformation::AffineTransform, point::Point) = (transformation * Subset(point)) |> first

Spaces.dimension(transformation::AffineTransform)::Integer = transformation.transformmatrix |> size |> first

pointgrouprepresentation(transformation::AffineTransform; rank::Integer = 1)::Matrix = pointgrouprepresentation(transformation.transformmatrix; rank=rank)

function findeigenfunction(transformation::AffineTransform;
    rankrange::UnitRange = 0:3, dimensionrange::UnitRange = 0:3, eigenvalue::Number = 1)::BasisFunction

    lookupsignatures::Base.Generator = (
        functionsignature(rank, dimension, eigenvalue |> Complex; denominator=transformation.eigenvaluehashdenominator)
        for (rank, dimension) in Iterators.product(rankrange, dimensionrange))

    for signature in lookupsignatures
        if haskey(transformation.eigenfunctions, signature)
            return transformation.eigenfunctions[signature]
        end
    end

    error("No eigenfunction is found for the provided constraints!")
end
export findeigenfunction

function Base.:*(transformation::AffineTransform, basisfunction::BasisFunction)::BasisFunction
    if basisfunction.rank == 0
        return basisfunction
    end

    if dimension(transformation) != dimension(basisfunction)
        error("Dimension mismatch!")
    end

    matrix::Matrix = pointgrouprepresentation(transformation; rank=basisfunction.rank)
    # TODO: Is this really the case?
    functionrep::Vector = transformation.antiunitary ? basisfunction |> rep |> conj : basisfunction |> rep  # Handling anti-unitary.
    transformed::Vector = matrix * functionrep

    return BasisFunction(transformed, basisfunction |> dimension, basisfunction.rank)
end

function pointgroupelements(pointgroup::AffineTransform; maxelements=128)::Vector{AffineTransform}
    identity::AffineTransform = pointgroup ^ 0
    elements::Vector{AffineTransform} = [identity]
    for n in 1:(maxelements - 1)
        current = pointgroup ^ n
        if current == identity
            break
        end
        push!(elements, current)
    end
    return elements
end
export pointgroupelements

pointgrouporder(pointgroup::AffineTransform; maxorder=128)::Integer = pointgroupelements(pointgroup; maxelements=maxorder) |> length
export pointgrouporder

function relativephase(target::BasisFunction, ref::BasisFunction)::Complex
    normalizedref::Vector = ref |> rep |> normalize
    normalizedtarget::Vector = target |> rep |> normalize
    trialphase::Number = dot(normalizedref, normalizedtarget)
    if !isapprox(dot(normalizedtarget ./ trialphase, normalizedref), 1)
        error("Target function is not a gauged version of the reference function.")
    end
    return trialphase
end
export relativephase

relativephase(ref::BasisFunction) = target::BasisFunction -> relativephase(target, ref)

struct Scale <: Transformation{Matrix{Float64}}
    rep::Matrix{Float64}
end
export Scale

Base.:convert(::Type{Matrix{Float64}}, source::Scale) = source.rep

Base.:inv(scale::Scale)::Scale = scale |> rep |> inv |> Scale
Base.:*(a::Scale, b::Scale)::Scale = Scale((a |> rep) * (b |> rep))
Base.:*(scale::Scale, space::RealSpace)::RealSpace = RealSpace((scale |> rep) * (space |> rep))
Base.:*(scale::Scale, space::MomentumSpace)::MomentumSpace = MomentumSpace((scale |> inv |> rep) * (space |> rep))
Base.:*(scale::Scale, point::Point)::Point = lineartransform(scale * (point |> getspace), point)
Base.:*(scale::Scale, subset::Subset)::Subset = Subset(scale * element for element in subset)

function Base.:*(scale::Scale, crystal::Crystal)::Crystal
    oldspace::RealSpace = crystal |> getspace
    snf = vcat(scale |> rep, crystal.sizes |> diagm) |> smith
    boundarysnf = snf.S[end - dimension(oldspace) + 1:end, 1:dimension(oldspace)] |> smith
    Δ::Matrix{Float64} = snf |> diagm
    newbasiscoords::Matrix{Float64} = boundarysnf.T * (Δ |> diag |> diagm) * snf.T
    blockingpoints::Base.Generator = (Point(collect(coord), oldspace) for coord in Iterators.product((0:size - 1 for size in diag(Δ))...))
    relativescale::Scale = Scale(newbasiscoords)
    scaledunitcell::Subset{Position} = Subset(relativescale * (a + b) for (a, b) in Iterators.product(blockingpoints, crystal.unitcell))
    return Crystal(scaledunitcell, diag(diagm(boundarysnf)))
end

end
