if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end

module Transformations

using LinearAlgebra, SmithNormalForm, OrderedCollections
using ..Spaces, ..Geometries, ..Quantum

export Scale, Symmetry, Irrep
export getirrep, groupelement, groupelements

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
    scaledunitcell::Subset{Position} = Subset(relativescale * (a + b) for (a, b) in [Iterators.product(blockingpoints, crystal.unitcell)...])
    return Crystal(scaledunitcell, diag(diagm(boundarysnf)))
end

function Base.:*(scale::Scale, crystalfock::FockSpace{Crystal})::FockMap
    crystal::Crystal = crystalof(crystalfock)
    newcrystal::Crystal = scale * crystal
    blockedregion::Subset{Position} = inv(scale) * newcrystal.unitcell
    BZ::Subset{Momentum} = brillouinzone(crystal)
    basismodes::Subset{Mode} = crystalfock |> rep |> first
    newBZ::Subset{Momentum} = brillouinzone(newcrystal)
    # Generate the Dict which keys each crystal fockspace partition by its momentum.
    momentumtopartition::Dict{Point, Subset{Mode}} = Dict(commonattr(part, :offset) => part for part in rep(crystalfock))
    momentummappings::Vector{Pair{Point, Point}} = [basispoint(scale * p) => p for p in BZ]
    mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point,Vector{Point}}()) do d,(k,v)
        mergewith!(append!, d, LittleDict(k=>[v]))
    end
    blockedcrystalpartitions::Subset{Subset{Mode}} = Subset(
        union([momentumtopartition[k] for k in mappingpartitions[scaled_k]]...) for scaled_k in newBZ)
    blockedcrystalordering::Dict{Mode, Integer} = Dict(mode => index for (index, mode) in enumerate(flatten(blockedcrystalpartitions)))
    blockedcrystalfock::FockSpace = FockSpace(blockedcrystalpartitions, blockedcrystalordering)

    blockedoffsets::Subset{Position} = Subset(pbc(crystal, p) |> latticeoff for p in blockedregion)
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

struct Symmetry <: AbstractSubset{Symmetry}
    group::Symbol
    pointgrouprep::Matrix{Float64}
    order::Integer
    index::Integer
    shift::Point
    irreps::Dict{Symbol, Irrep}
end

Symmetry(; group::Symbol, pointgrouprep::Matrix{Float64}, order::Integer, shift::Point, irreps::Dict{Symbol, Irrep}) = Symmetry(
    group, pointgrouprep, order,
    1, # This symmetry element is assumed to be the generator of the group, thus making it the second element (The first is the identity, with index = 0).
    shift, irreps)

Spaces.spaceof(symmetry::Symmetry) = symmetry.shift |> spaceof

function groupelement(generator::Symmetry, index::Integer = 0)::Symmetry
    @assert(generator.index == 1, "The argument is not a generator!")
    pointgrouprep::Matrix{Float64} = generator.pointgrouprep ^ index
    irreps::Dict{Symbol, Irrep} = Dict(orbital => irrep ^ index for (orbital, irrep) in generator.irreps)
    return Symmetry(generator.group, pointgrouprep, generator.order, index, generator.shift, irreps)
end

groupelements(generator::Symmetry)::Vector{Symmetry} = [groupelement(generator, index) for index in 0:generator.order - 1]

function getirrep(symmetry::Symmetry, orbital::Symbol)::Irrep
    @assert(haskey(symmetry.irreps, orbital), "Symmetry `$(symmetry.group)` does not have orbital `$(orbital)`!")
    return symmetry.irreps[orbital]
end

function Base.:*(space::AffineSpace, symmetry::Symmetry)::Symmetry
    realspace::RealSpace = convert(RealSpace, space)
    relativebasis::Matrix{Float64} = (realspace |> basis |> inv) * (symmetry |> spaceof |> basis)
    pointgrouprep::Matrix{Float64} = relativebasis * symmetry.pointgrouprep * inv(relativebasis)
    shift::Point = Point(((symmetry.shift |> pos)' * (relativebasis |> inv))[1, :], space)
    return Symmetry(symmetry.group, pointgrouprep, symmetry.order, symmetry.index, shift, symmetry.irreps)
end

function Base.:*(symmetry::Symmetry, region::Subset{Position})::Subset{Position}
    nativesymmetry::Symmetry = spaceof(region) * symmetry
    shift::Point = nativesymmetry.shift
    return Subset(Point(nativesymmetry.pointgrouprep * (point |> pos), point |> spaceof) + shift for point in region)
end

function Base.:*(symmetry::Symmetry, zone::Subset{Momentum})::Subset{Momentum}
    kspace::MomentumSpace = zone |> spaceof
    realspace::RealSpace = convert(RealSpace, kspace)
    pointgrouprep::Matrix{Float64} = (realspace * symmetry).pointgrouprep
    kspacerep::Matrix{Float64} = Matrix(pointgrouprep |> transpose |> inv)
    return Subset(Point(kspacerep * (k |> pos), kspace) for k in zone)
    # TODO: Anti-unitary is not handled.
end

Base.:*(symmetry::Symmetry, point::P) where {P <: Point} = (symmetry * Subset(point)) |> first

function Base.:*(symmetry::Symmetry, mode::Mode)::FockMap
    # This is used to correct the :pos attribute, since the :pos as a Point will be symmetrized,
    # which the basis point set might not include the symmetrized :pos. Thus we would like to set
    # the :pos to its corresponding basis point, and offload the difference to :offset.
    function correctsymmetrizedmode(mode::Mode)::Mode
        currentoffset::Point = getattr(mode, :offset)
        currentpos::Point = getattr(mode, :pos)
        actualpos::Point = basispoint(currentpos)
        adjustoffset::Point = currentpos - actualpos
        return setattr(mode, :offset => currentoffset + adjustoffset, :pos => basispoint(currentpos))
    end
    newattrs::Dict{Symbol, Any} = Dict(mode.attrs)
    foreach(p -> newattrs[p.first] = symmetry * p.second, Iterators.filter(p -> hasmethod(*, Tuple{Symmetry, typeof(p.second)}), mode.attrs))
    newmode::Mode = Mode(newattrs) |> correctsymmetrizedmode
    irrep::Irrep = getirrep(symmetry, orbital(symmetry.group, mode))
    return FockMap(newmode |> FockSpace, mode |> FockSpace, Dict((newmode, mode) => irrep |> rep))
end

function Base.:*(symmetry::Symmetry, subset::Subset{Mode})::FockMap
    fockmap::FockMap = focksum([symmetry * mode for mode in subset])
    return FockMap(fockmap.outspace |> flattensubspaces, fockmap.inspace |> flattensubspaces, fockmap |> rep)
end

Base.:*(symmetry::Symmetry, fockspace::FockSpace{<: Any})::FockMap = focksum([symmetry * subspace for subspace in fockspace |> rep])

function Base.:*(symmetry::Symmetry, crystalfock::FockSpace{Crystal})::FockMap
    homefock::FockSpace = unitcellfock(crystalfock)
    homefocktransform::FockMap = symmetry * homefock
    momentumsubspaces::Dict{Momentum, FockSpace} = crystalsubspaces(crystalfock)
    fouriermap::FockMap = fourier(crystalfock, homefock)
    transformedfouriermap::FockMap = fourier(crystalfock, homefocktransform.outspace)
    transformer::FockMap = focksum([
        rows(transformedfouriermap, momentumsubspaces[(symmetry * k) |> basispoint]) * homefocktransform * rows(fouriermap, fockspace)'
        for (k, fockspace) in momentumsubspaces])
    crystal::Crystal = crystalof(crystalfock)
    return FockMap(
        transformer,
        outspace=FockSpace(transformer.outspace, reflected=crystal),
        inspace=FockSpace(transformer.inspace, reflected=crystal))
end

end
