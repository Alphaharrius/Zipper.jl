if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end

module Transformations

using LinearAlgebra, SmithNormalForm, OrderedCollections
using ..Spaces, ..Geometries, ..Quantum

export Recipient, Scale

abstract type Transformation{T} <: Element{T} end

struct Recipient{T}
    object::T
    parameters::Dict{Symbol}

    Recipient(object::T, parameters::Vararg{Pair{Symbol}}) where {T} = new{T}(object, Dict(parameters...))
end

function recipientparam(recipient::Recipient, key::Symbol)
    @assert(haskey(recipient.parameters, key), "Missing parameter of `$(key)`!")
    return recipient.parameters[key]
end

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

Base.:*(scale::Scale, subset::Subset{T}) where {T} = Subset([scale * element for element in subset])

function Base.:*(scale::Scale, crystal::Crystal)::Crystal
    oldspace::RealSpace = spaceof(crystal)
    snf = smith(vcat(rep(scale), diagm(crystal.sizes)))
    boundarysnf = smith(snf.S[end - dimension(oldspace) + 1:end, 1:dimension(oldspace)])
    Δ::Matrix{Float64} = diagm(snf)
    newbasiscoords::Matrix{Float64} = boundarysnf.T * diagm(diag(Δ)) * snf.T
    blockingpoints::Vector{Point} = [Point(collect(coord), oldspace) for coord in [Iterators.product([0:size - 1 for size in diag(Δ)]...)...]]
    relativescale::Scale = Scale(newbasiscoords)
    scaledunitcell::Subset{Point} = Subset([relativescale * (a + b) for (a, b) in [Iterators.product(blockingpoints, crystal.unitcell)...]])
    return Crystal(scaledunitcell, diag(diagm(boundarysnf)))
end

function Base.:*(scale::Scale, recipient::Recipient{FockSpace{CrystalFock}})::FockMap # For fockspace that is sparse under momentum partitions.
    fockspace::FockSpace{CrystalFock} = recipient.object
    crystal::Crystal = recipientparam(recipient, :crystal)
    newcrystal::Crystal = scale * crystal
    blockedregion::Subset{Point} = inv(scale) * newcrystal.unitcell
    BZ::Subset{Point} = brillouinzone(crystal)
    basismodes::Subset{Mode} = fockspace |> rep |> first
    newBZ::Subset{Point} = brillouinzone(newcrystal)
    # Generate the Dict which keys each fockspace partition by its momentum.
    momentumtopartition::Dict{Point, Subset{Mode}} = Dict(getattr(first(part), :offset) => part for part in rep(fockspace))
    momentummappings::Vector{Pair{Point, Point}} = [basispoint(scale * p) => p for p in BZ]
    mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point,Vector{Point}}()) do d,(k,v)
        mergewith!(append!, d, LittleDict(k=>[v]))
    end
    blockedcrystalpartitions::Subset{Subset{Mode}} = Subset([
        union([momentumtopartition[k] for k in mappingpartitions[scaled_k]]...) for scaled_k in newBZ])
    blockedcrystalordering::Dict{Mode, Integer} = Dict(mode => index for (index, mode) in enumerate(flatten(blockedcrystalpartitions)))
    blockedcrystalfock::FockSpace = FockSpace(blockedcrystalpartitions, blockedcrystalordering)

    blockedoffsets::Subset{Point} = Subset([pbc(crystal, p) |> latticeoff for p in blockedregion])
    blockedfock::FockSpace = FockSpace(spanoffset(basismodes, blockedoffsets))

    restrictedfourier::FockMap = fourier(BZ, blockedfock)
    volumeratio::Float64 = sqrt(vol(crystal) / vol(newcrystal))
    permutedmap::FockMap = Quantum.permute(restrictedfourier, blockedcrystalfock, restrictedfourier.inspace) / volumeratio
    function repack_fourierblocks(sourcemap::FockMap, scaled_k::Point, partition::Subset{Mode})::FockMap
        mappart::FockMap = rows(sourcemap, FockSpace(partition))
        inmodes::Subset{Mode} = Subset([
            setattr(m, :groups => ModeGroup(transformed, "scaled"), :index => i, :offset => scaled_k, :pos => scale * convert(Point, m))
            for (i, m) in enumerate(orderedmodes(mappart.inspace))])
        return FockMap(mappart.outspace, FockSpace(inmodes), rep(mappart))
    end
    mapblocks::Vector{FockMap} = [repack_fourierblocks(permutedmap, scaled_k, partition)
                                  for (scaled_k, partition) in Iterators.zip(newBZ, rep(blockedcrystalfock))]
    blockmap::FockMap = focksum(mapblocks)
    return FockMap(blockmap.outspace, FockSpace(blockmap.inspace, T=CrystalFock), rep(blockmap))'
end

end