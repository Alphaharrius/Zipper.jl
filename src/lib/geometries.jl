"""
    distance(a::Point, b::Point)::Float64

Compute the distance between two points `a` & `b` within the same parent `AffineSpace`.
"""
distance(a::Point, b::Point)::Float64 = norm(a - b)
export distance

function interpolate(from::Point, to::Point, count::T)::Array{Point} where {T <: Integer}
    @assert(getspace(from) == getspace(to))
    march::Point = (to - from) / (count + 1)
    points::Array{Point} = [from + march * n for n in 0:(count + 1)]
    points[end] = to
    return points
end
export interpolate

getorigin(space::AffineSpace)::Point = Point(zeros(Float64, dimension(space)), space)
export getorigin

getcenter(region::Subset{<:Point})::Point = sum(p for p in region) / length(region)
export getcenter

function getradius(region::Subset{<:Point}; metricspace::AffineSpace = region |> getspace |> euclidean)::Real
    center::Point = region |> getcenter
    return maximum(lineartransform(metricspace, center - p) |> norm for p in region)
end
export getradius

BoundaryCondition(space::RealSpace, sizes::Integer...) = BoundaryCondition(
    space, Matrix(1.0I, sizes|>length, sizes|>length), sizes|>collect)

getspace(bc::BoundaryCondition)::RealSpace = bc.localspace

function Base.hash(bc::BoundaryCondition)::UInt
    spacehash = hash(bc.localspace)
    ruleshash = hash(map(v -> Rational{Int64}(round(v * 10000000)) // 10000000, bc.rules))
    boundshash = hash(bc.bounds)
    return hash((spacehash, ruleshash, boundshash))
end

Base.:(==)(a::BoundaryCondition, b::BoundaryCondition)::Bool = (
    a.localspace == b.localspace && isapprox(a.rules, b.rules) && a.bounds == b.bounds)

getbounds(bc::BoundaryCondition) = bc.bounds
export getbounds

Base.:show(io::IO, crystal::Crystal) = print(io, string("$(crystal|>typeof)"))

getunitcell(crystal::Crystal) = crystal.unitcell
export getunitcell

getbc(crystal::Crystal) = crystal.bc
export getbc

"""
    getbarecrystal(crystal::Crystal)::Crystal

A bare crystal defines a crystal that have a singlar zero offset position as the unit-cell, and the provided 
argument `crystal` is used to supply the `RealSpace` and the size information to the bare crystal.
"""
getbarecrystal(crystal::Crystal)::Crystal = Crystal(crystal|>getspace|>getorigin|>Subset, crystal|>getbc)
return getbarecrystal

Base.:(==)(a::Crystal, b::Crystal)::Bool = a.bc == b.bc && a.unitcell == b.unitcell

latticeoff(point::Point)::Point = Point([trunc(v) for v in point |> vec], getspace(point))
export latticeoff

getorigin(crystal::Crystal)::Offset = crystal|>getspace|>getorigin

"""
    basispoint(point::Point)::Point

Convert a point back to the corresponding basis point within the unitcell.
"""
function basispoint(point::Point)::Point
    rationalized::Vector = [hashablereal(v) for v in point |> vec]
    return Point([mod(v |> numerator, v |> denominator) / denominator(v) for v in rationalized], point |> getspace)
end
export basispoint

Zipper.:getspace(crystal::Crystal) = crystal|>getunitcell|>getspace
Zipper.:dimension(crystal::Crystal) = crystal|>getspace|>dimension

mesh(sizes::Vector{Int64})::Matrix{Int64} = hcat([collect(tup) for tup in collect(Iterators.product([0:(d - 1) for d in sizes]...))[:]]...)

vol(crystal::Crystal)::Integer = crystal|>brillouinzone|>length
export vol

@memoize function getboundsize(crystal::Crystal)
    _, S, _ = crystal.bc.rules|>dosnf
    return crystal.bc.bounds./(S|>diag)
end
export getboundsize

@memoize function brillouinzone(crystal::Crystal)::Subset{Momentum}
    bc = crystal|>getbc
    imesh = hcat((
        transpose(bc.rules)*collect(v)./bc.bounds 
        for v in Base.product((0:d-1 for d in bc.bounds)...)|>collect)...)
    kspace = convert(MomentumSpace, crystal|>getspace)
    return Subset(kspace*collect(p)|>basispoint for p in eachcol(imesh))
end
export brillouinzone

@memoize function computemomentummatrix(crystal::Crystal)
    function computekmatrix(batch)
        matrix::SparseMatrixCSC = spzeros(crystal|>dimension, crystal|>vol)
        for (n, k) in batch
            matrix[:, n] = k|>euclidean|>vec
        end
        return matrix
    end

    batchsize::Integer = ((crystal|>vol) / getmaxthreads())|>ceil
    batches = Iterators.partition((el for el in crystal|>brillouinzone|>enumerate), batchsize)
    matrixparts = paralleltasks(
        name="computemomentummatrix",
        tasks=(()->computekmatrix(batch) for batch in batches),
        count=getmaxthreads())|>parallel

    momentummatrix::SparseMatrixCSC = spzeros(crystal|>dimension, crystal|>vol)
    watchprogress(desc="computemomentummatrix")
    for part in matrixparts
        momentummatrix[:, :] += part
        updateprogress()
    end
    unwatchprogress()

    return momentummatrix
end

buildregion(space::AffineSpace, sizes::Integer...) = Subset(
    space*collect(v) for v in Base.product((0:d-1 for d in sizes)...))

function buildregion(crystal::Crystal, sizes::Integer...)
    @assert dimension(crystal) == length(sizes)

    space = crystal|>getspace
    offsets = buildregion(space, sizes...)
    return Subset(r+b for (r, b) in Base.product(offsets, crystal|>getunitcell))
end
export buildregion

function getsphericalregion(; crystal::Crystal, radius::Real, metricspace::RealSpace)
    genradius = ceil(Int, radius * 1.5)
    genlengths = (genradius*2/norm(metricspace*bv)|>ceil|>Integer for bv in crystal|>getspace|>getbasisvectors)
    seedoffsets = buildregion(crystal|>getspace, genlengths...)
    seedoffsets = seedoffsets .- getcenter(seedoffsets)
    seedregion = Subset(r+b for (r, b) in Base.product(seedoffsets, crystal|>getunitcell))
    return Subset(r for r in seedregion if norm(metricspace*r) <= radius)
end
export getsphericalregion

geometricfilter(f, metricspace::AffineSpace) = p -> (metricspace * p) |> f
export geometricfilter

linearscale(space::AffineSpace) = log.([v|>norm for v in space|>getbasisvectors])|>mean|>exp
export linearscale

function orthodirections(vector::Point)
    Q, _ = vector|>vec|>qr
    space::RealSpace = vector|>getspace
    return Iterators.drop((Q[:, n] ∈ space for n in axes(Q, 2)), 1)
end
export orthodirections

function getcrosssection(; crystal::Crystal, normalvector::Offset, radius::Real, metricspace::RealSpace = crystal|>getspace|>orthospace, minbottomheight::Real = 0)
    height::Real = (*(metricspace, normalvector)|>norm)
    sphericalregion::Region = getsphericalregion(crystal=crystal, radius=sqrt(height^2 + radius^2)*2, metricspace=metricspace)

    normaldirection::Offset = *(metricspace, normalvector)|>normalize
    function crosssectionfilter(point::Point)::Bool
        metricpoint::Point = metricspace * point
        iswithinheight::Bool = minbottomheight <= dot(metricpoint, normaldirection) <= height
        orthoreminder::Point = metricpoint - dot(metricpoint, normaldirection) * normaldirection
        iswithinradius::Bool = norm(orthoreminder) < radius
        return iswithinheight && iswithinradius
    end

    rawregion::Region = sphericalregion|>filter(crosssectionfilter)
    proximityregion::Region = rawregion .- normalvector

    return rawregion - intersect(rawregion, proximityregion)
end
export getcrosssection
