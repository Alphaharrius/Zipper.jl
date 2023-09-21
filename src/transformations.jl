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

struct PointGroupTransformation <: Transformation{Matrix{Float64}}
    localspace::AffineSpace
    shiftvector::Vector{Float64}
    pointgroupmatrix::Matrix{Float64}
    eigenfunctions::Dict{Tuple, BasisFunction}
    antiunitary::Bool
end
export PointGroupTransformation

function PointGroupTransformation(
    pointgroupmatrix::Matrix, shiftvector::Vector = zeros(Float64, pointgroupmatrix |> size |> first);
    dimension::Integer = pointgroupmatrix |> size |> first, antiunitary::Bool = false,
    localspace::AffineSpace = euclidean(RealSpace, dimension))::PointGroupTransformation

    eigenfunctions::Dict{Tuple, BasisFunction} = Dict(eigenfunctionsignature(swave, 1) => swave)
    invariantfound::Bool = false
    for rank in 1:2 # The maximum rank will be 2 since it is sufficient for the considered point groups.
        rankmatrix::Matrix = pointgrouprepresentation(pointgroupmatrix; rank=rank)
        evals, evecs = rankmatrix |> eigen
        invariantindexs = findall(v -> isapprox(v, 1), evals)
        if invariantindexs isa Nothing continue end
        invariantfound = true
        for (n, eval) in enumerate(evals)
            basis::BasisFunction = BasisFunction(antiunitary ? evecs[:, n] |> conj : evecs[:, n], dimension, rank) |> normalize
            if basis |> iszero continue end
            eigenfunctions[eigenfunctionsignature(basis, eval)] = basis
        end
    end

    if invariantfound
        return PointGroupTransformation(localspace, shiftvector, pointgroupmatrix, eigenfunctions, antiunitary)
    end
    error("Invalid point group transformation, invariant basis is not found!")
end

isometrictransformation(
    pointgroupmatrix::Matrix;
    dimension::Integer = pointgroupmatrix |> size |> first,
    localspace::RealSpace = euclidean(RealSpace, dimension),
    referencepoint::Position = localspace |> origin,
    antiunitary::Bool = false)::PointGroupTransformation = PointGroupTransformation(pointgroupmatrix, transformationshift(pointgroupmatrix, localspace, referencepoint);
    dimension=dimension, localspace=localspace, antiunitary=antiunitary)
export isometrictransformation

recenter(transformation::PointGroupTransformation, center::Position)::PointGroupTransformation = PointGroupTransformation(
    transformation.pointgroupmatrix, transformationshift(transformation.pointgroupmatrix, transformation.localspace, center);
    dimension=transformation |> dimension, localspace=transformation |> spaceof, antiunitary=transformation.antiunitary)
export recenter

recenter(center::Position) = transformation::PointGroupTransformation -> recenter(transformation, center)

translation(
    shiftvector::Vector;
    dimension::Integer = shiftvector |> length,
    localspace::RealSpace = euclidean(RealSpace, dimension),
    antiunitary::Bool = false)::PointGroupTransformation = PointGroupTransformation(
    Matrix{Float64}(I, dimension, dimension), shiftvector; dimension=dimension, localspace=localspace, antiunitary=antiunitary)
export translation

function affinematrix(transformation::PointGroupTransformation)::Matrix{Float64}
    shiftrow::Vector = vcat(transformation.shiftvector, [1])
    leftcolumns::Matrix = vcat(
        transformation.pointgroupmatrix, zeros(Float64, transformation.pointgroupmatrix |> size |> first, 1) |> transpose)
    return hcat(leftcolumns, shiftrow)
end
export affinematrix

function transformationshift(pointgroupmatrix::Matrix, localspace::AffineSpace, reference::Point)::Vector
    return (reference |> pos) - (pointgroupmatrix * (lineartransform(localspace, reference) |> pos))
end

Spaces.:spaceof(transformation::PointGroupTransformation)::AffineSpace = transformation.localspace

Base.:convert(::Type{Matrix{Float64}}, source::PointGroupTransformation) = source |> affinematrix

Base.:(==)(a::PointGroupTransformation, b::PointGroupTransformation) = isapprox(a |> rep, b |> rep)

function Base.:*(space::RealSpace, transformation::PointGroupTransformation)::PointGroupTransformation
    if space |> dimension != transformation |> dimension
        error("Dimension mismatch!")
    end
    relativebasis::Matrix = (space |> basis |> inv) * (transformation |> getspace |> basis)
    pointgroupmatrix::Matrix = relativebasis * (transformation.pointgroupmatrix) * (relativebasis |> inv)
    shiftvector::Vector = lineartransform(space, transformation.localspace & transformation.shiftvector) |> pos
    return PointGroupTransformation(
        pointgroupmatrix, shiftvector;
        localspace=space, dimension=(transformation |> dimension), antiunitary=transformation.antiunitary)
end

function Base.:*(a::PointGroupTransformation, b::PointGroupTransformation)::PointGroupTransformation
    dim::Integer = a |> dimension
    if dim != b |> dimension
        error("Dimension mismatch!")
    end
    localb::PointGroupTransformation = (a |> spaceof) * b
    affinematrix::Matrix = (a |> rep) * (localb |> rep)
    pointgroupmatrix::Matrix = affinematrix[1:dim, 1:dim]
    shiftvector::Vector = affinematrix[1:dim, end]
    antiunitary::Bool = a.antiunitary ⊻ localb.antiunitary
    return PointGroupTransformation(
        pointgroupmatrix, shiftvector;
        localspace=(a |> spaceof), dimension=dim, antiunitary=antiunitary)
end

Base.:*(transformation::PointGroupTransformation, space::RealSpace)::RealSpace = RealSpace((transformation.pointgroupmatrix) * (space |> rep))

Base.:^(source::PointGroupTransformation, exponent::Number)::PointGroupTransformation = reduce(*, repeat([source], exponent))

function Base.:*(transformation::PointGroupTransformation, region::Subset{Position})::Subset{Position}
    nativetransformation::PointGroupTransformation = (region |> spaceof) * transformation
    transformed::Base.Generator = (
        Point(((nativetransformation |> rep) * vcat(point |> pos, [1]))[1:end - 1], point |> spaceof)
        for point in region)
    return Subset(transformed)
end

function Base.:*(transformation::PointGroupTransformation, zone::Subset{Momentum})::Subset{Momentum}
    kspace::MomentumSpace = zone |> getspace
    realspace::RealSpace = convert(RealSpace, kspace)
    nativetransformation::PointGroupTransformation = realspace * transformation
    kspacerep::Matrix = Matrix(nativetransformation.pointgroupmatrix |> transpose |> inv)
    return Subset(Point(kspacerep * (k |> pos), kspace) for k in zone)
end

Base.:*(transformation::PointGroupTransformation, point::Point) = (transformation * Subset(point)) |> first

Spaces.dimension(transformation::PointGroupTransformation)::Integer = transformation.pointgroupmatrix |> size |> first

pointgrouprepresentation(matrix::Matrix; rank::Integer = 1)::Matrix = reduce(kron, repeat([matrix], rank))
pointgrouprepresentation(transformation::PointGroupTransformation; rank::Integer = 1)::Matrix = pointgrouprepresentation(transformation.pointgroupmatrix; rank=rank)
export pointgrouprepresentation

function findeigenfunction(transformation::PointGroupTransformation;
    rankrange::UnitRange = 0:2, dimensionrange::UnitRange = 0:3, eigenvalue::Number = 1)::BasisFunction

    lookupsignatures::Base.Generator = (
        eigenfunctionsignature(rank, dimension, eigenvalue) for (rank, dimension) in Iterators.product(rankrange, dimensionrange))

    for signature in lookupsignatures
        if haskey(transformation.eigenfunctions, signature)
            return transformation.eigenfunctions[signature]
        end
    end

    error("No eigenfunction is found for the provided constraints!")
end
export findeigenfunction

function Base.:*(transformation::PointGroupTransformation, basisfunction::BasisFunction)::BasisFunction
    if basisfunction.rank == 0
        return basisfunction
    end

    if dimension(transformation) != dimension(basisfunction)
        error("Dimension mismatch!")
    end

    matrix::Matrix = pointgrouprepresentation(transformation; rank=basisfunction.rank)
    functionrep::Vector = transformation.antiunitary ? basisfunction |> rep |> conj : basisfunction |> rep  # Handling anti-unitary.
    transformed::Vector = matrix * functionrep

    return BasisFunction(transformed, basisfunction |> dimension, basisfunction.rank)
end

function pointgroupelements(transformation::PointGroupTransformation; maxelements=128)::Vector{PointGroupTransformation}
    elements::Vector{PointGroupTransformation} = [transformation]
    for n in 2:maxelements
        current = transformation ^ n
        if current == transformation
            break
        end
        push!(elements, current)
    end
    return elements
end
export pointgroupelements

pointgrouporder(transformation::PointGroupTransformation; maxorder=128)::Integer = pointgroupelements(transformation; maxelements=maxorder) |> length
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
