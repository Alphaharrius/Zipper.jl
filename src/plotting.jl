module Plotting

using Plotly, ColorTypes, LinearAlgebra, Compat
using ..Spaces, ..Geometries, ..Quantum

export visualize

function visualize(fockmap::FockMap; title::String = "", rowrange = :, colrange = :)
    source = rep(fockmap)[rowrange, colrange]
    layout = Layout(yaxis=attr(autorange="reversed"))
    subplots = [plot(heatmap(z=real(source)), layout) plot(heatmap(z=imag(source)), layout)]
    relayout!(subplots, title_text=title)
    subplots
end

function visualize(regions::Subset{<: Point}...; title::String = "", visualspace::AffineSpace)
    function generateregionplot(region::Subset{<: Point})
        positions::Array{Vector} = [lineartransform(visualspace, point) |> pos for point in region]
        return visualizevectorpositions(positions)
    end
    # How to make all axis having the same scale:
    # https://stackoverflow.com/questions/52863305/plotly-scatter3d-how-can-i-force-3d-axes-to-have-the-same-scale-aspect-ratio
    layout = Layout(title=title, scene=attr(aspectmode="data"))
    fig = plot([region |> generateregionplot for region in regions], layout)
    relayout!(fig, template="simple_white")
    fig
end

function visualizevectorpositions(positions::Array{Vector})
    values::Matrix = hcat(positions...)
    @assert(0 < size(values, 1) <= 3)
    # Pad all positions to 3-dimensional to generalize all cases.
    padded_values::Matrix = zeros(3, size(values, 2))
    copyto!(view(padded_values, 1:size(values, 1), :), values)
    trace = scatter3d(x=padded_values[1,:], y=padded_values[2,:], z=padded_values[3,:], mode="markers")
    return trace
end

function visualize(spectrum::Vector{Pair{Mode, T}}; title::String = "") where {T <: Number}
    𝑁::Int64 = length(spectrum)
    ∑𝑝::Vector{Point} = [getattr(pair.first, :offset) + getattr(pair.first, :pos) for pair in spectrum]
    𝑀ₚ::Matrix{Float64} = hcat(map(𝑝 -> pos(lineartransform(euclidean(RealSpace, dimension(𝑝)), 𝑝)), ∑𝑝)...)
    markerpositions::Matrix{Float64} = zeros(3, 𝑁)
    copyto!(view(markerpositions, 1:size(𝑀ₚ, 1), :), 𝑀ₚ)
    sizes::Vector{Float64} = [abs(pair.second) for pair in spectrum]
    markersizes::Vector{Float64} = sizes / norm(sizes) * 120
    colors::Vector{RGB{Float64}} = [convert(RGB{Float64}, HSV(angle(pair.second) / 2π * 360, 1, 1)) for pair in spectrum]
    markercolors::Vector{Tuple{Float32, Float32, Float32, Float32}} = map(c -> Tuple([c.r, c.g, c.b, 1.0] * 255), colors)
    trace = scatter3d(
        x=markerpositions[1, :], y=markerpositions[2, :], z=markerpositions[3, :], mode="markers",
        marker=attr(
            symbol="circle",
            size=markersizes,
            color=markercolors))
    layout::Layout = Layout(title=title, scene=attr(aspectmode="data"))
    fig = plot([trace], layout)
    relayout!(fig, template="simple_white")
    fig
end

function visualize_crystalspectrum_2d(spectrum::CrystalSpectrum)
    kspectrum::Dict{Momentum} = Dict(k => ([spectrum.eigenvalues[m] for m in modes] |> sort) for (k, modes) in spectrum.eigenmodes)
    mesh = spectrum.crystal |> brillouinmesh
    plottingspectrum = stack(map(k -> kspectrum[k], mesh))
    return [surface(z=plottingspectrum[n, :, :]) for n in axes(plottingspectrum, 1)]
end

CRYSTALSPECPLOTTINGFUNCTIONS::Dict = Dict(2 => visualize_crystalspectrum_2d)

function visualize(spectrum::CrystalSpectrum; title="")
    plottingdimension::Integer = spectrum.crystal |> dimension
    if !haskey(CRYSTALSPECPLOTTINGFUNCTIONS, plottingdimension)
        @error("Plotting crystal spectrum of dimension $(plottingdimension) is not supported!")
    end
    surfaces = spectrum |> CRYSTALSPECPLOTTINGFUNCTIONS[plottingdimension]
    layout::Layout = Layout(title=title)
    plot(surfaces, layout)
end

end
