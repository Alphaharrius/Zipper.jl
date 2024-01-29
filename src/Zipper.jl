module Zipper

using OrderedCollections, Base.Iterators, SparseArrays
using LinearAlgebra, Combinatorics, Statistics, SmithNormalForm

include("core.jl")

include("interfaces/core.jl")
include("interfaces/spaces.jl")
include("interfaces/transformations.jl")
include("interfaces/geometries.jl")

include("lib/spaces.jl")
include("lib/geometries.jl")
include("lib/transformations.jl")
include("lib/quantum.jl")
include("lib/quantumtransformations.jl")
include("lib/physics.jl")
include("lib/renormalization.jl")

using Plotly, ColorTypes, Compat
include("lib/plotting.jl")

end
