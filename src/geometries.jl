if !isdefined(Main, :Spaces) include("spaces.jl") end

module Geometries

using LinearAlgebra, OrderedCollections
using ..Spaces

export Crystal
export distance, interpolate, origin, radius, pbc, basispoint, latticeoff, resize, mesh, vol, latticepoints, sitepoints, brillouinzone, brillouinmesh, geometricalfilter, circularfilter

"""
    distance(a::Point, b::Point)::Float64

Compute the distance between two points `a` & `b` within the same parent `AffineSpace`.
"""
distance(a::Point, b::Point)::Float64 = norm(a - b)

function interpolate(from::Point, to::Point, count::T)::Array{Point} where {T <: Integer}
    @assert(spaceof(from) == spaceof(to))
    march::Point = (to - from) / (count + 1)
    points::Array{Point} = [from + march * n for n in 0:(count + 1)]
    points[end] = to
    return points
end

origin(space::AffineSpace)::Point = Point(zeros(Float64, dimension(space)), space)

center(region::Subset{<:Point})::Point = sum(p for p in region) / length(region)
export center

radius(region::Subset, center::Point)::Float64 = maximum(distance(center, point) for point in rep(region))

struct Crystal <: AbstractSubset{Crystal}
    unitcell::Subset{Position}
    sizes::Vector{Int64}
end

Base.:show(io::IO, crystal::Crystal) = print(io, string("$(crystal |> typeof)(sizes=$(crystal.sizes))"))

pbc(crystal::Crystal, point::Point)::Point = Point([mod(p, s) for (p, s) in zip(pos(point), crystal.sizes)], spaceof(point))

latticeoff(point::Point)::Point = Point([trunc(v) for v in pos(point)], spaceof(point))

basispoint(point::Point)::Point = Point([mod(numerator(v), denominator(v)) / denominator(v) for v in rpos(point)], spaceof(point))

Spaces.spaceof(crystal::Crystal) = spaceof(crystal.unitcell)
Spaces.dimension(crystal::Crystal) = crystal.sizes |> length

resize(crystal::Crystal, sizes::Vector{Int64})::Crystal = Crystal(crystal.unitcell, sizes)

mesh(sizes::Vector{Int64})::Matrix{Int64} = hcat([collect(tup) for tup in collect(Iterators.product([0:(d - 1) for d in sizes]...))[:]]...)

vol(crystal::Crystal)::Integer = prod(crystal.sizes)

function latticepoints(crystal::Crystal)::Subset{Position}
    real_space::RealSpace = spaceof(crystal.unitcell)
    crystal_mesh::Matrix{Int64} = mesh(crystal.sizes)
    tiled_sizes::Matrix{Int64} = hcat([crystal.sizes for i in 1:size(crystal_mesh, 2)]...)
    recentered_mesh::Matrix{Float64} = (crystal_mesh - tiled_sizes / 2)
    return Subset(Point(pos, real_space) for pos in eachcol(recentered_mesh))
end

sitepoints(crystal::Crystal)::Subset{Position} = Subset(
    latticepoint + basispoint for latticepoint in latticepoints(crystal) for basispoint in crystal.unitcell)

function brillouinzone(crystal::Crystal)::Subset{Momentum}
    momentum_space::MomentumSpace = convert(MomentumSpace, spaceof(crystal.unitcell))
    crystal_mesh::Matrix{Int64} = mesh(crystal.sizes)
    tiled_sizes::Matrix{Int64} = hcat([crystal.sizes for i in 1:size(crystal_mesh, 2)]...)
    recentered_mesh::Matrix{Float64} = crystal_mesh ./ tiled_sizes
    return Subset(Point(pos, momentum_space) for pos in eachcol(recentered_mesh))
end

function brillouinmesh(crystal::Crystal)::Array{Point}
    kspace::MomentumSpace = spaceof(crystal)
    return [Point(collect(p) ./ crystal.sizes, kspace) for p in Iterators.product([0:d - 1 for d in crystal.sizes]...)]
end

geometricalfilter(f, center::Point) = input -> f(lineartransform(spaceof(center), convert(Point, input)), center)

circularfilter(center::Point, radius::Float64) = geometricalfilter((point, center) -> norm(point - center) <= radius, center)

end
