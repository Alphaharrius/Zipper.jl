module Plotting
export show

include("spaces.jl")
include("points.jl")

using PlotlyJS, ..Spaces, ..Points

function show(source::SpaceSubset, space::ClassicalSpace)
    ordered_representation::Array{Point} = collect(Point, source.representation)
    target_representation::Array{Point} = [map_to(point, space) for point in ordered_representation]
    x_values::Array{Rational{Int64}} = [point.position[1] for point in target_representation]
    y_values::Array{Rational{Int64}} = [point.position[2] for point in target_representation]

    points_scatter = scatter(
        mode="markers",
        x=x_values,
        y=y_values,
        marker=attr(
            color="LightSkyBlue",
            size=20,
            line=attr(
                color="MediumPurple",
                width=2
            )
        ),
        showlegend=false
    )

    plot([points_scatter])
end

end