include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/transformations.jl")
include("../src/quantum.jl")
include("../src/plotting.jl")

using ..Spaces, ..Geometries, ..Transformations, ..Plotting

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
c6 = triangular * isometrictransformation([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)]) |> recenter(triangular & [1., 0.])
c6 |> rep

p = triangular & [0, 0]

points = Subset(t * p for t in repeat([c6], 6) |> cumprod) + Subset(triangular & [1., 0.])

visualize(points, visualspace=euclidean(RealSpace, 2))

c6 * p

triangular |> rep