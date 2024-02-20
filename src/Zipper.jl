module Zipper

using OrderedCollections, Base.Iterators, SparseArrays
using LinearAlgebra, Combinatorics, Statistics, SmithNormalForm

using ExprTools, ConcurrentCollections
include("lib/memoization.jl")

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

include("lib/quantum/mode.jl")
include("lib/quantum/fockspace.jl")
include("lib/quantum/normalfock.jl")
include("lib/quantum/regionfock.jl")
include("lib/quantum/crystalfock.jl")
include("lib/quantum.jl")
include("lib/quantum/transformations.jl")

include("lib/physics.jl")
include("lib/renormalization.jl")
include("lib/extendedrestrictions.jl")
include("lib/gmerarenormalization.jl")

using Plotly, ColorTypes, Compat
include("lib/plotting.jl")

using JSON, DataFrames, CSV, Dates
include("lib/lzw.jl")
include("lib/fio.jl")
include("lib/fioinjection.jl")

end
