if !isdefined(Main, :Spaces) include("spaces.jl") end

module Geometries

using ..Spaces

origin(space::AffineSpace)::Point = Point(zeros(Float64, dimension(space)), space)

radius(region::Subset, center::Point)::Float64 = maximum(distance(center, point) for point in rep(region))

struct Crystal
    unit_cell::Subset{Point}
    sizes::Vector{Int64}
end

resize(crystal::Crystal, sizes::Vector{Int64})::Crystal = Crystal(crystal.unit_cell, sizes)

mesh(sizes::Vector{Int64})::Matrix{Int64} = hcat([collect(tup) for tup in collect(Iterators.product([0:(d - 1) for d in sizes]...))[:]]...)

vol(crystal::Crystal)::Integer = prod(crystal.sizes)

function points(crystal::Crystal)::Subset{Point}
    real_space::RealSpace = spaceof(crystal.unit_cell)
    crystal_mesh::Matrix{Int64} = mesh(crystal.sizes)
    mesh_points::Array{Point} = [Point(pos, real_space) for pos in eachcol(crystal_mesh)]
    points::Set{Point} = Set{Point}([mesh_point + point for mesh_point in mesh_points for point in rep(crystal.unit_cell)])
    return Subset(points)
end

function brillouin_zone(crystal::Crystal)::Subset{Point}
    momentum_space::MomentumSpace = convert(MomentumSpace, spaceof(crystal.unit_cell))
    crystal_mesh::Matrix{Int64} = mesh(crystal.sizes)
    tiled_sizes::Matrix{Int64} = hcat([crystal.sizes for i in 1:size(crystal_mesh, 2)]...)
    recentered_mesh::Matrix{Float64} = (crystal_mesh - tiled_sizes / 2) ./ tiled_sizes
    points::Set{Point} = Set{Point}([Point(pos, momentum_space) for pos in eachcol(recentered_mesh)])
    return Subset(points)
end

export Crystal
export origin, radius, resize, mesh, vol, points, brillouin_zone

end
