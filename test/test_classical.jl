include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/plotting.jl")

using ..Spaces
using ..Geometries
using ..Quantum
using ..Plotting

space = RealSpace([1. 1.5; 0.5 0.5])
point0 = Point([1//1, 3//1], space)
point1 = Point([-8//1, 3//1], space)
point2 = Point([0//1, 0//1], space)

region::Subset{Point} = union(convert(Subset, point0), convert(Subset, point1), convert(Subset, point2))

modes::Subset{Mode} = quantize(region, 2)
fock_space::FockSpace = FockSpace(modes)
