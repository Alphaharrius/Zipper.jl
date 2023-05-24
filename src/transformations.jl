if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end

module Transformations

using LinearAlgebra, SmithNormalForm, OrderedCollections
using ..Spaces, ..Geometries, ..Quantum

export Scale, Symmetry, Irrep
export irrepof

abstract type Transformation{T} <: Element{T} end

Base.:*(transformation::T, element::F) where {T <: Transformation, F <: Element} = error(
    "Composition of `$(typeof(transformation))` with `$(typeof(element))` is not defined.")

struct Scale <: Transformation{Matrix{Float64}}
    rep::Matrix{Float64}
end

Base.:inv(scale::Scale)::Scale = Scale(inv(rep(scale)))

Base.:convert(::Type{Matrix{Float64}}, source::Scale) = source.rep

Base.:*(a::Scale, b::Scale)::Scale = Scale(rep(a) * rep(b))

Base.:*(scale::Scale, space::RealSpace) = convert(typeof(space), rep(scale) * rep(space))
Base.:*(scale::Scale, space::MomentumSpace) = convert(typeof(space), rep(inv(scale)) * rep(space))

Base.:*(scale::Scale, point::Point)::Point = lineartransform(scale * spaceof(point), point)

Base.:*(scale::Scale, subset::Subset{T}) where {T} = Subset(scale * element for element in subset)

function Base.:*(scale::Scale, crystal::Crystal)::Crystal
    oldspace::RealSpace = spaceof(crystal)
    snf = smith(vcat(rep(scale), diagm(crystal.sizes)))
    boundarysnf = smith(snf.S[end - dimension(oldspace) + 1:end, 1:dimension(oldspace)])
    Δ::Matrix{Float64} = diagm(snf)
    newbasiscoords::Matrix{Float64} = boundarysnf.T * (Δ |> diag |> diagm) * snf.T
    blockingpoints::Vector{Point} = [Point(collect(coord), oldspace) for coord in [Iterators.product([0:size - 1 for size in diag(Δ)]...)...]]
    relativescale::Scale = Scale(newbasiscoords)
    scaledunitcell::Subset{Point} = Subset(relativescale * (a + b) for (a, b) in [Iterators.product(blockingpoints, crystal.unitcell)...])
    return Crystal(scaledunitcell, diag(diagm(boundarysnf)))
end

function Base.:*(scale::Scale, crystalfock::FockSpace{Crystal})::FockMap
    crystal::Crystal = crystalof(crystalfock)
    newcrystal::Crystal = scale * crystal
    blockedregion::Subset{Point} = inv(scale) * newcrystal.unitcell
    BZ::Subset{Point} = brillouinzone(crystal)
    basismodes::Subset{Mode} = crystalfock |> rep |> first
    newBZ::Subset{Point} = brillouinzone(newcrystal)
    # Generate the Dict which keys each crystal fockspace partition by its momentum.
    momentumtopartition::Dict{Point, Subset{Mode}} = Dict(getattr(first(part), :offset) => part for part in rep(fockspace))
    momentummappings::Vector{Pair{Point, Point}} = [basispoint(scale * p) => p for p in BZ]
    mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point,Vector{Point}}()) do d,(k,v)
        mergewith!(append!, d, LittleDict(k=>[v]))
    end
    blockedcrystalpartitions::Subset{Subset{Mode}} = Subset(
        union([momentumtopartition[k] for k in mappingpartitions[scaled_k]]...) for scaled_k in newBZ)
    blockedcrystalordering::Dict{Mode, Integer} = Dict(mode => index for (index, mode) in enumerate(flatten(blockedcrystalpartitions)))
    blockedcrystalfock::FockSpace = FockSpace(blockedcrystalpartitions, blockedcrystalordering)

    blockedoffsets::Subset{Point} = Subset(pbc(crystal, p) |> latticeoff for p in blockedregion)
    blockedfock::FockSpace = FockSpace(spanoffset(basismodes, blockedoffsets))

    restrictedfourier::FockMap = fourier(BZ, blockedfock)
    volumeratio::Float64 = sqrt(vol(crystal) / vol(newcrystal))
    permutedmap::FockMap = Quantum.permute(restrictedfourier, outspace=blockedcrystalfock, inspace=restrictedfourier.inspace) / volumeratio
    function repack_fourierblocks(sourcemap::FockMap, scaled_k::Point, partition::Subset{Mode})::FockMap
        mappart::FockMap = rows(sourcemap, FockSpace(partition))
        inmodes::Subset{Mode} = Subset(
            setattr(m, :offset => scaled_k, :pos => scale * convert(Point, m))
            for m in orderedmodes(mappart.inspace))
        return FockMap(mappart.outspace, FockSpace(inmodes), rep(mappart))
    end
    mapblocks::Vector{FockMap} = [repack_fourierblocks(permutedmap, scaled_k, partition)
                                  for (scaled_k, partition) in Iterators.zip(newBZ, rep(blockedcrystalfock))]
    blockmap::FockMap = focksum(mapblocks)
    return FockMap(blockmap.outspace, FockSpace(blockmap.inspace, reflected=newcrystal), rep(blockmap))'
end

struct Irrep <: Element{ComplexF64}
    rep::ComplexF64
end

Base.:show(io::IO, irrep::Irrep) = print(io, string("$(typeof(irrep))($(rep(irrep)))"))

Base.:^(irrep::Irrep, n)::Irrep = Irrep(rep(irrep) ^ n)

Base.:(==)(a::Irrep, b::Irrep)::Bool = isapprox(rep(a), rep(b))

Base.:*(a::Irrep, b::Irrep)::Irrep = Irrep(rep(a) * rep(b))

Base.:convert(::Type{ComplexF64}, source::Irrep) = source.rep

function Base.:hash(source::Irrep)::UInt
    denom::Int64 = 10000000
    rationalizedrep::Complex{Rational{Int64}} = (
        Rational{Int64}(round((source |> rep |> real) * denom)) // denom + Rational{Int64}(round((source |> rep |> imag) * denom)) // denom * 1im)
    return hash(rationalizedrep)
end

Base.:isequal(a::Irrep, b::Irrep)::Bool = a == b

struct Symmetry <: Transformation{Vector{Matrix{Float64}}}
    group::Symbol
    firstorderrep::Matrix{Float64}
    order::Integer
    center::Point
    irreps::Dict{Symbol, Irrep}

    Symmetry(; group::Symbol, rep::Matrix{Float64}, order::Integer, center::Point, irreps::Vector{Pair{Symbol, Irrep}}) = new(
        group, rep, order, center, Dict(irreps...))
end

Base.:convert(::Type{Matrix{Float64}}, source::Symmetry) = [symmetry.firstorderrep ^ n for n in 0:symmetry.order - 1]

Spaces.:dimension(symmetry::Symmetry)::Integer = size(symmetry.firstorderrep, 1)

Base.:*(symmetry::Symmetry, region::Subset{Point})::Subset{Point} = Subset(
    Point((p |> spaceof |> rep |> inv) * sym * ((p - symmetry.center) |> euclidean |> pos), spaceof(p)) + symmetry.center
    for sym in groupreps(symmetry) for p in region)

Base.:*(symmetry::Symmetry, point::Point)::Subset{Point} = symmetry * Subset(point)

end
