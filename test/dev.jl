include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/transformations.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/quantumtransformations.jl")
include("../src/renormalization.jl")
include("../src/plotting.jl")

using PlotlyJS, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using ..Spaces, ..Geometries, ..Quantum, ..Transformations, ..Plotting, ..QuantumTransformations, ..Physical, ..Renormalization

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = triangular & [1/3, 2/3]
pb = triangular & [2/3, 1/3]
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

c6 * findeigenfunction(c6, eigenvalue=-1)

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [24, 24])
reciprocalhashcalibration(crystal.sizes)

function fillregion(; crystal::Crystal, seedsizes::Vector{<:Integer}, seedcenter::Point, regionpredicate::Function)
    seedcrystal::Crystal = Crystal(crystal.unitcell, seedsizes)
    seedregion::Region = (seedcrystal |> sitepoints) + seedcenter
    return filter(regionpredicate, seedregion)
end

function circularpredicate(center::Point, radius::Real)::Function
    scaledeuclidean::AffineSpace = center |> getspace |> orthogonalspace
    return p -> norm(scaledeuclidean * p - scaledeuclidean * center) < radius
end

visualize(Crystal(unitcell, [6, 6]) |> sitepoints, visualspace=euclidean(RealSpace, 2))

visualize(fillregion(crystal=crystal, seedsizes=[4,4], seedcenter=triangular |> origin, regionpredicate=circularpredicate(triangular |> origin, 1.5)), visualspace=euclidean(RealSpace, 2))
