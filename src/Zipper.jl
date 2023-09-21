module Zipper

using LinearAlgebra, SparseArrays, OrderedCollections, Base.Iterators, SmithNormalForm, Combinatorics, Statistics

include("spaces.jl")
using .Spaces
export Element, AbstractSpace, AffineSpace, RealSpace, MomentumSpace, AbstractSubset, Point, Position, Momentum, Subset
export rep, euclidean, basis, dimension, getspace, rpos, pos, lineartransform, fourier_coef, distance, flatten, members, subsetunion

include("geometries.jl")
include("transformations.jl")
include("quantum.jl")
include("physical.jl")

include("quantumtransformations.jl")
include("renormalization.jl")

using Plotly, ColorTypes, Compat

include("plotting.jl")

end
