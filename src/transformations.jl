if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end

module Transformations

using LinearAlgebra, SmithNormalForm
using ..Spaces, ..Geometries

export Scale

abstract type Transformation{T} <: Element{T} end

Base.:*(transformation::T, element::F) where {T <: Transformation, F <: Element} = error(
    "Composition of `$(typeof(transformation))` with `$(typeof(element))` is not defined.")

struct Scale <: Transformation{Matrix{Float64}}
    rep::Matrix{Float64}
end

Base.:convert(::Type{Matrix{Float64}}, source::Scale) = source.rep

Base.:*(a::Scale, b::Scale)::Scale = Scale(rep(a) * rep(b))

Base.:*(scale::Scale, space::T) where {T <: AffineSpace} = convert(typeof(space), rep(scale) * rep(space))

function Base.:*(scale::Scale, point::Point)::Point
    oldspace::AffineSpace = spaceof(point)
    scaledspace::AffineSpace = scale * oldspace
    scalevector::Vector{Float64} = diag(rep(scale))
    return Point(pos(point) ./ scalevector, scaledspace)
end

Base.:*(scale::Scale, subset::Subset{T}) where {T} = Subset([scale * element for element in subset])

function Base.:*(scale::Scale, crystal::Crystal)::Crystal
    oldspace::RealSpace = spaceof(crystal)
    snf = smith(vcat(rep(scale), diagm(crystal.sizes)))
    boundary_snf = smith(snf.S[end - dimension(oldspace) + 1:end, 1:dimension(oldspace)])
    Δ::Matrix{Float64} = diagm(snf)
    newbasiscoords::Matrix{Float64} = boundary_snf.T * diagm(diag(Δ)) * snf.T
    blockingpoints::Vector{Point} = [Point(collect(coord), oldspace) for coord in [Iterators.product([0:size - 1 for size in diag(Δ)]...)...]]
    relativescale::Scale = Scale(newbasiscoords)
    scaledunitcell::Subset{Point} = Subset([relativescale * (a + b) for (a, b) in [Iterators.product(blockingpoints, crystal.unitcell)...]])
    return Crystal(scaledunitcell, diag(diagm(boundary_snf)))
end

end