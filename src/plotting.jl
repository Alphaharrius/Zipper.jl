if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end

module Plotting

using PlotlyJS, ColorTypes, LinearAlgebra
using ..Spaces, ..Geometries, ..Quantum

function visualize(fockmap::FockMap; title::String="", rowrange=:, colrange=:)
    subplots = [plot(heatmap(z=real(rep(fockmap)[rowrange, colrange]))) plot(heatmap(z=imag(rep(fockmap)[rowrange, colrange])))]
    relayout!(subplots, title_text=title)
    subplots
end

function visualize_region(title::String, region::Subset{Point}, visualization_space::AffineSpace)
    ordered_points::Array{Point} = [lineartransform(visualization_space, point) for point in rep(region)]
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
    𝑀ₚ::Matrix{Float64} = hcat(map(𝑝 -> pos(lineartransform(euclidean(RealSpace, dimension(𝑝)), 𝑝)), ∑𝑝)...)
    markerpositions::Matrix{Float64} = zeros(3, 𝑁)
    copyto!(view(markerpositions, 1:size(𝑀ₚ, 1), :), 𝑀ₚ)
    sizes::Vector{Float64} = [abs(pair.second) for pair in spectrum]
    markersizes::Vector{Float64} = sizes / norm(sizes) * 80
    colors::Vector{RGB{Float64}} = [convert(RGB{Float64}, HSV(angle(pair.second) / 2π * 360, 1, 1)) for pair in spectrum]
    markercolors::Vector{Tuple{Float32, Float32, Float32}} = map(c -> Tuple([c.r, c.g, c.b] * 255), colors)
    trace = scatter3d(x=markerpositions[1, :], y=markerpositions[2, :], z=markerpositions[3, :], mode="markers", marker=attr(
        size=markersizes,
        color=markercolors))
    layout::Layout = Layout(title=title, scene=attr(aspectmode="data"))
    plot([trace], layout)
end

export visualize_region, visualize_vector_positions, visualize_spectrum, visualize

end
