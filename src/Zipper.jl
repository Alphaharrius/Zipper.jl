module Zipper

using OrderedCollections, Base.Iterators, SparseArrays
using LinearAlgebra, Combinatorics, Statistics, SmithNormalForm

using ExprTools
include("lib/memorization.jl")

using ProgressMeter
include("lib/parallel.jl")

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
include("lib/extendedrestrictions.jl")
include("lib/gmerarenormalization.jl")

using Plotly, ColorTypes, Compat
include("lib/plotting.jl")

using JSON, DataFrames, CSV, Dates
include("lib/fio.jl")
include("lib/fioinjection.jl")

end
