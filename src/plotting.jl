if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end

module Plotting

using PlotlyJS, ColorTypes
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
    function generatetrace(element::Pair{Mode, T}) where {T <: Number}
        offset::Point = getattr(element.first, :offset)
        pos::Point = linear_transform(spaceof(offset), getattr(element.first, :pos))
        point::Point = offset + pos
        position::Vector{Float64} = Spaces.pos(linear_transform(euclidean(RealSpace, dimension(point)), point))
        padded_position::Vector{Float64} = vcat(position, zeros(Float64, 3 - length(position)))
        value::ComplexF64 = ComplexF64(element.second)
    
        # Generate the RGB color based on the complex value.
        hue = angle(value) / 2π
        rgb::RGB{Float64} = convert(RGB{Float64}, HSV(hue, 1, 1))
        colortext::String = "rgba($(rgb.r * 255), $(rgb.g * 255), $(rgb.b * 255), 0.3)"

        return scatter3d(x=[padded_position[1]], y=[padded_position[2]], z=[padded_position[3]], mode="marker", marker=attr(
            size=abs(value) * 50,
            color=colortext))
    end

    layout::Layout = Layout(title=title, scene=attr(aspectmode="data"))
    plot([generatetrace(element) for element in spectrum], layout)
end

export visualize_region, visualize_vector_positions, visualize_spectrum

end
