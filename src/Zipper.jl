module Zipper

using DataStructures, OrderedCollections, Base.Iterators, SparseArrays, PrettyTables
using AbstractAlgebra, LinearAlgebra, Combinatorics, Statistics

include("tools/collectiontools.jl")

using ExprTools, ConcurrentCollections
include("tools/memoization.jl")

using ProgressMeter
include("tools/parallel.jl")

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
include("lib/quantum/fockmap.jl")
include("lib/quantum/sparsefockmap.jl")
include("lib/quantum/crystalfockmap.jl")
include("lib/quantum/decompositions.jl")
include("lib/quantum/fouriertransform.jl")
include("lib/quantum/regionstate.jl")
include("lib/quantum/extensions.jl")
include("lib/quantum/transformations.jl")

include("lib/physics.jl")
include("lib/renormalization.jl")
include("lib/extendedrestrictions.jl")
include("lib/gmerarenormalization.jl")

using Plots
include("lib/plotting.jl")

using JSON, DataFrames, CSV, Dates
include("tools/lzw.jl")
include("tools/fio.jl")
include("lib/fioadapt.jl")

end
