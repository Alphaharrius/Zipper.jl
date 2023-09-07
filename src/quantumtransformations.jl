if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end
if !isdefined(Main, :Transformations) include("transformations.jl") end

module QuantumTransformations

using OrderedCollections
using ..Spaces, ..Geometries, ..Quantum, ..Transformations

end
