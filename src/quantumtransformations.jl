if !isdefined(Main, :Spaces) include("spaces.jl") end
if !isdefined(Main, :Geometries) include("geometries.jl") end
if !isdefined(Main, :Quantum) include("quantum.jl") end
if !isdefined(Main, :Transformations) include("transformations.jl") end

module QuantumTransformations

using OrderedCollections
using ..Spaces, ..Geometries, ..Quantum, ..Transformations

function Base.:*(scale::Scale, crystalfock::FockSpace{Crystal})::FockMap
    crystal::Crystal = crystalfock |> crystalof
    scaledcrystal::Crystal = scale * crystal
    unscaledblockedregion::Subset{Position} = (scale |> inv) * scaledcrystal.unitcell
    bz::Subset{Momentum} = crystal |> brillouinzone
    basismodes::Subset{Mode} = crystalfock |> rep |> first
    scaledbz::Subset{Momentum} = scaledcrystal |> brillouinzone

    momentummappings::Base.Generator = (basispoint(scale * p) => p for p in bz)
    mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point, Vector{Point}}()) do d,(k,v)
        mergewith!(append!, d, LittleDict(k=>[v]))
    end

    ksubsets::Dict{Momentum, Subset{Mode}} = crystalfock |> crystalsubsets
    scaledfock::FockSpace = ((ksubsets[k] for k in mappingpartitions[kscaled]) |> subsetunion |> FockSpace
        for kscaled in scaledbz) |> fockspaceunion

    unscaledblockedlatticeoffsets::Subset{Position} = Subset(pbc(crystal, p) |> latticeoff for p in unscaledblockedregion)
    unscaledblockedunitcellfock::FockSpace = spanoffset(basismodes, unscaledblockedlatticeoffsets)

    restrictedfourier::FockMap = fourier(bz, unscaledblockedunitcellfock)
    volumeratio::Number = (crystal |> vol) / (scaledcrystal |> vol)
    permutedfourier::FockMap = Quantum.permute(restrictedfourier, outspace=scaledfock) / sqrt(volumeratio)

    function repackfourierblocks(source::FockMap, kscaled::Momentum, partition::Subset{Mode})::FockMap
        partitionrows::FockMap = rows(source, partition |> FockSpace)
        inspace::FockSpace = (Subset(setattr(mode, :offset => kscaled, :pos => scale * convert(Point, mode))
            for mode in partitionrows.inspace |> orderedmodes)
            |> FockSpace)
        return FockMap(partitionrows.outspace, inspace, partitionrows |> rep)
    end

    repackedblocks::Base.Generator = (
        repackfourierblocks(permutedfourier, kscaled, partition)
        for (kscaled, partition) in Iterators.zip(scaledbz, scaledfock |> rep))
    blocking::FockMap = directsum(repackedblocks)
    return FockMap(blocking.outspace, FockSpace(blocking.inspace, reflected=scaledcrystal), blocking |> rep)'
end
end
