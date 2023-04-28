if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end

module Plotting

using PlotlyJS, ColorTypes, LinearAlgebra
using ..Spaces, ..Geometries, ..Quantum

function visualize_region(title::String, region::Subset{Point}, visualization_space::AffineSpace)
    ordered_points::Array{Point} = [linear_transform(visualization_space, point) for point in rep(region)]
    positions::Array{Vector} = [pos(point) for point in ordered_points]
    visualize_vector_positions(title, positions)
end

function visualize_vector_positions(title::String, positions::Array{Vector})
    values::Matrix = hcat(positions...)
    @assert(0 < size(values, 1) <= 3)
    # Pad all positions to 3-dimensional to generalize all cases.
    padded_values::Matrix = zeros(3, size(values, 2))
    copyto!(view(padded_values, 1:size(values, 1), :), values)
    trace = scatter3d(x=padded_values[1,:], y=padded_values[2,:], z=padded_values[3,:], mode="markers")
    # How to make all axis having the same scale:
    # https://stackoverflow.com/questions/52863305/plotly-scatter3d-how-can-i-force-3d-axes-to-have-the-same-scale-aspect-ratio
    layout = Layout(title=title, scene=attr(aspectmode="data"))
    plot([trace], layout)
end

function visualize_spectrum(title::String, spectrum::Vector{Pair{Mode, T}}) where {T <: Number}
    𝑁::Int64 = length(spectrum)
    ∑𝑝::Vector{Point} = [getattr(pair.first, :offset) + getattr(pair.first, :pos) for pair in spectrum]
    𝑀ₚ::Matrix{Float64} = hcat(map(𝑝 -> pos(linear_transform(euclidean(RealSpace, dimension(𝑝)), 𝑝)), ∑𝑝)...)
    sizes::Vector{Float64} = [abs(pair.second) for pair in spectrum]
    markersizes::Vector{Float64} = sizes / norm(sizes) * 50
    colors::Vector{RGB{Float64}} = [convert(RGB{Float64}, HSV(angle(pair.second) / 2π, 1, 1)) for pair in spectrum]
    markercolors::Vector{Tuple{Int16, Int64, Int16}} = map(c -> Tuple([c.r, c.g, c.b] * 255), colors)
    trace = scatter3d(x=𝑀ₚ[1, :], y=𝑀ₚ[2, :], z=𝑀ₚ[3, :], mode="marker", marker=attr(size=markersizes, color=markercolors))
    layout::Layout = Layout(title=title, scene=attr(aspectmode="data"))
    plot([trace], layout)
end

export visualize_region, visualize_vector_positions, visualize_spectrum

end
