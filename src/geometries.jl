if !isdefined(Main, :Spaces) include("spaces.jl") end

module Geometries

using ..Spaces

radius(region::Subset, center::Point)::Float64 = maximum(distance(center, point) for point in rep(region))

struct Crystal
    unit_cell::Subset{Point}
    shape::Vector{Int64}
end

resize(crystal::Crystal, shape::Vector{Int64})::Crystal = Crystal(crystal.unit_cell, shape)

mesh(shape::Vector{Int64})::Matrix{Int64} = hcat([collect(tup) for tup in collect(Iterators.product([0:(d - 1) for d in shape]...))[:]]...)

vol(crystal::Crystal)::Integer = prod(crystal.shape)

function points(crystal::Crystal)::Subset{Point}
    real_space::RealSpace = space_of(crystal.unit_cell)
    crystal_mesh::Matrix{Int64} = mesh(crystal.shape)
    mesh_points::Array{Point} = [Point(pos, real_space) for pos in eachcol(crystal_mesh)]
    points::Set{Point} = Set{Point}([mesh_point + point for mesh_point in mesh_points for point in rep(crystal.unit_cell)])
    return Subset(points)
end

function brillouin_zone(crystal::Crystal)::Subset{Point}
    momentum_space::MomentumSpace = convert(MomentumSpace, space_of(crystal.unit_cell))
    crystal_mesh::Matrix{Int64} = mesh(crystal.shape)
    tiled_shape::Matrix{Int64} = hcat([crystal.shape for i in 1:size(crystal_mesh, 2)]...)
    recentered_mesh::Matrix{Float64} = (crystal_mesh - tiled_shape / 2) ./ tiled_shape
    points::Set{Point} = Set{Point}([Point(pos, momentum_space) for pos in eachcol(recentered_mesh)])
    return Subset(points)
end

export Crystal
export radius, resize, mesh, vol, points, brillouin_zone

end
