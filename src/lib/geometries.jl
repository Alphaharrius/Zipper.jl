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

Base.:show(io::IO, crystal::Crystal) = print(io, string("$(crystal |> typeof)(sizes=$(crystal.sizes))"))

Zipper.:getregion(crystal::Crystal) = sitepoints(crystal)

getunitcell(crystal::Crystal) = crystal.unitcell
export getunitcell

"""
    getbarecrystal(crystal::Crystal)::Crystal

A bare crystal defines a crystal that have a singlar zero offset position as the unit-cell, and the provided 
argument `crystal` is used to supply the `RealSpace` and the size information to the bare crystal.
"""
getbarecrystal(crystal::Crystal)::Crystal = Crystal(crystal|>getspace|>getorigin|>Subset, crystal|>size)
return getbarecrystal

Base.:size(crystal::Crystal) = crystal.sizes

Base.:(==)(a::Crystal, b::Crystal)::Bool = a.sizes == b.sizes && a.unitcell == b.unitcell

"""
    pbc(crystal::Crystal, point::Point)::Point

Apply periodic boundary conditions to a point `point` within a `Crystal`.
To do: fix mod in pbc
"""
pbc(crystal::Crystal, point::Point)::Point = Point([mod(p, s) for (p, s) in zip(point |> vec, crystal.sizes)], getspace(point))
export pbc

"""
    pbc(crystal::Crystal)::Function

Return a function that applies periodic boundary conditions to a point within a `Crystal`.
"""
pbc(crystal::Crystal)::Function = p -> pbc(crystal, p)

latticeoff(point::Point)::Point = Point([trunc(v) for v in point |> vec], getspace(point))
export latticeoff

getorigin(crystal::Crystal)::Offset = crystal|>getspace|>getorigin

getkbasis(crystal::Crystal, k::Momentum) = vec(k) .* size(crystal)
export getkbasis

kindexoperator(crystal::Crystal) = [1, (size(crystal)[1:end-1]|>cumprod)...]
export kindexoperator

Base.getindex(crystal::Crystal, k::Momentum)::Integer = kindexoperator(crystal)' * getkbasis(crystal, k) + 1

Base.getindex(crystal::Crystal, i::Integer) = brillouinzone(crystal)[i]

Base.:in(k::Momentum, crystal::Crystal) = (
    getspace(k) == convert(MomentumSpace, getspace(crystal)) && getkbasis(crystal, k)|>sum|>isinteger)

"""
    basispoint(point::Point)::Point

Convert a point back to the corresponding basis point within the unitcell.
"""
function basispoint(point::Point)::Point
    rationalized::Vector = [hashablereal(v) for v in point |> vec]
    return Point([mod(v |> numerator, v |> denominator) / denominator(v) for v in rationalized], point |> getspace)
end
export basispoint

Zipper.:getspace(crystal::Crystal) = getspace(crystal.unitcell)
Zipper.:dimension(crystal::Crystal) = crystal.sizes |> length

resize(crystal::Crystal, sizes::Vector{Int64})::Crystal = Crystal(crystal.unitcell, sizes)
export resize

mesh(sizes::Vector{Int64})::Matrix{Int64} = hcat([collect(tup) for tup in collect(Iterators.product([0:(d - 1) for d in sizes]...))[:]]...)

vol(crystal::Crystal)::Integer = prod(crystal.sizes)
export vol

@memoize function latticepoints(crystal::Crystal)::Subset{Offset}
    realspace::RealSpace = getspace(crystal.unitcell)
    crystalmesh::Matrix{Int64} = mesh(crystal.sizes)
    return Subset(Point(pos, realspace) for pos in eachcol(crystalmesh))
end
export latticepoints

@memoize sitepoints(crystal::Crystal)::Subset{Offset} = Subset(
    latticepoint + basispoint for latticepoint in latticepoints(crystal) for basispoint in crystal.unitcell)
export sitepoints

@memoize function brillouinzone(crystal::Crystal)::Subset{Momentum}
    momentum_space::MomentumSpace = convert(MomentumSpace, getspace(crystal.unitcell))
    crystal_mesh::Matrix{Int64} = mesh(crystal.sizes)
    tiled_sizes::Matrix{Int64} = hcat([crystal.sizes for i in 1:size(crystal_mesh, 2)]...)
    recentered_mesh::Matrix{Float64} = crystal_mesh ./ tiled_sizes
    return Subset(Point(pos, momentum_space) for pos in eachcol(recentered_mesh))
end
export brillouinzone

function brillouinmesh(crystal::Crystal)::Array{Point}
    kspace::MomentumSpace = getspace(crystal)
    return [Point(collect(p) ./ crystal.sizes, kspace) for p in Iterators.product([0:d - 1 for d in crystal.sizes]...)]
end
export brillouinmesh

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

function getsphericalregion(; crystal::Crystal, radius::Real, metricspace::RealSpace)
    generatingradius::Integer = ceil(Int, radius * 1.5) # Multiply by 1.5 to ensure all unitcell points fits.
    generatinglength::Integer = generatingradius * 2
    generatingcrystal::Crystal = Crystal(crystal|>getunitcell, [generatinglength, generatinglength])
    crystalregion::Region = generatingcrystal|>sitepoints
    centeredregion::Region = crystalregion .- (crystalregion|>getcenter)
    return Subset(point for point in centeredregion if norm(metricspace * point) <= radius)
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
