if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end

module Plotting

using PlotlyJS
using ..Spaces, ..Geometries

function visualize_region(region::Subset{Point}, visualization_space::AffineSpace)
    ordered_points::Array{Point} = [linear_transform(visualization_space, point) for point in rep(region)]
    positions::Array{Vector} = [pos(point) for point in ordered_points]
    visualize_vector_positions(positions)
end

function visualize_vector_positions(positions::Array{Vector})
    values::Matrix = hcat(positions...)
    @assert(0 < size(values, 1) <= 3)
    # Pad all positions to 3-dimensional to generalize all cases.
    padded_values::Matrix = zeros(3, size(values, 2))
    copyto!(view(padded_values, 1:size(values, 1), :), values)
    plot([scatter3d(x=padded_values[1,:], y=padded_values[2,:], z=padded_values[3,:], mode="markers")])
end

export visualize_region, visualize_vector_positions

end
