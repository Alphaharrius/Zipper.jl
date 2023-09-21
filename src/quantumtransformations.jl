module QuantumTransformations

using OrderedCollections
using ..Spaces, ..Geometries, ..Quantum, ..Transformations

function Base.:*(scale::Scale, crystalfock::FockSpace{Crystal})::FockMap
    crystal::Crystal = crystalfock |> getcrystal
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

function Base.:*(transformation::PointGroupTransformation, subset::Subset{Mode})::FockMap
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

    modemapping::Dict{Mode, Mode} = Dict()

    function modesymmetrize(mode::Mode)::Mode
        newattrs::Dict{Symbol, Any} = Dict(mode.attrs)
        # TODO: There are some attributes that are not meant to be transformed with the mode.
        filterpredicate = p -> hasmethod(*, Tuple{PointGroupTransformation, p.second |> typeof})
        foreach(p -> newattrs[p.first] = transformation * p.second, Iterators.filter(filterpredicate, mode.attrs))
        newmode::Mode = Mode(newattrs) |> correctsymmetrizedmode
        modemapping[newmode] = mode
        return newmode
    end

    connections::Dict{Tuple{Mode, Mode}, ComplexF64} = Dict()

    function rebaseorbital(mode::Mode)::Mode
        frommode::Mode = modemapping[mode]
        tomode::Mode = mode |> setattr(:orbital => (frommode |> getorbital()))
        connections[(tomode, frommode)] = (mode |> getorbital()) |> relativephase(mode |> getorbital())
        return tomode
    end

    outmodes::Subset{Mode} = Subset(mode |> modesymmetrize |> rebaseorbital for mode in subset)
    
    return FockMap(outmodes |> FockSpace, subset |> FockSpace, connections)
end

Base.:*(transformation::PointGroupTransformation, fockspace::FockSpace)::FockMap = transformation * (fockspace |> orderedmodes)

function Base.:*(transformation::PointGroupTransformation, crystalfock::FockSpace{Crystal})::FockMap
    homefock::FockSpace = crystalfock |> unitcellfock
    homefocktransform::FockMap = transformation * homefock
    ksubspaces::Dict{Momentum, FockSpace} = crystalfock |> crystalsubspaces
    fouriertransform::FockMap = fourier(crystalfock, homefock)
    transformedfourier::FockMap = fourier(crystalfock, homefocktransform.outspace)
    transform::FockMap = directsum(
        rows(transformedfourier, ksubspaces[(transformation * k) |> basispoint]) * homefocktransform * rows(fouriertransform, fockspace)'
        for (k, fockspace) in ksubspaces)
    crystal::Crystal = crystalfock |> getcrystal
    return FockMap(transform, outspace=FockSpace(transform.outspace, reflected=crystal), inspace=FockSpace(transform.inspace, reflected=crystal))
end

end
