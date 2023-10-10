function Base.:*(scale::Scale, crystalfock::FockSpace{Crystal})::FockMap
    crystal::Crystal = crystalfock |> getcrystal
    scaledcrystal::Crystal = scale * crystal
    unscaledblockedregion::Subset{Offset} = (scale |> inv) * scaledcrystal.unitcell
    bz::Subset{Momentum} = crystal |> brillouinzone
    basismodes::Subset{Mode} = crystalfock |> rep |> first
    scaledbz::Subset{Momentum} = scaledcrystal |> brillouinzone

    momentummappings::Base.Generator = (basispoint(scale * p) => p for p in bz)
    mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point, Vector{Point}}()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end

    ksubsets::Dict{Momentum, Subset{Mode}} = crystalfock |> crystalsubsets
    scaledfock::FockSpace = ((ksubsets[k] for k in mappingpartitions[kscaled]) |> subsetunion |> FockSpace
        for kscaled in scaledbz) |> fockspaceunion

    # TODO: Use of latticeoff requires revisit, since latticeoff only truncates the decimals.
    unscaledblockedlatticeoffsets::Subset{Offset} = Subset(pbc(crystal, p) |> latticeoff for p in unscaledblockedregion)
    unscaledblockedunitcellfock::FockSpace = spanoffset(basismodes, unscaledblockedlatticeoffsets)

    restrictedfourier::FockMap = fourier(bz, unscaledblockedunitcellfock)
    volumeratio::Number = (crystal |> vol) / (scaledcrystal |> vol)
    permutedfourier::FockMap = Zipper.permute(restrictedfourier, outspace=scaledfock) / sqrt(volumeratio)

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

function Base.:*(transformation::AffineTransform, regionfock::FockSpace{Region})::FockMap
    # This is used to correct the :pos attribute, since the :pos as a Point will be symmetrized,
    # which the basis point set might not include the symmetrized :pos. Thus we would like to set
    # the :pos to its corresponding basis point, and offload the difference to :offset.
    function correctsymmetrizedmode(mode::Mode)::Mode
        actualposition::Offset = mode |> getattr(:pos)
        basisposition::Offset = actualposition |> basispoint
        offset::Offset = actualposition - basisposition
        return mode |> setattr(:pos => basisposition) |> setattr(:offset => offset)
    end

    modemapping::Dict{Mode, Mode} = Dict()

    function mergepositions(mode::Mode)::Mode
        latticeoffset::Point = mode |> getattr(:offset)
        latticeoffset isa Offset || error("Transforming a mode based on momentum must be done with crystal fockspace!")
        actualposition::Offset = getattr(mode, :pos) + latticeoffset
        return mode |> setattr(:pos => actualposition) |> removeattr(:offset)
    end

    function modesymmetrize(mode::Mode)::Mode
        fixedmode = mode |> mergepositions
        newattrs::Dict{Symbol, Any} = Dict(fixedmode.attrs)
        # TODO: There are some attributes that are not meant to be transformed with the mode.
        filterpredicate = p -> hasmethod(*, Tuple{AffineTransform, p.second |> typeof})
        foreach(p -> newattrs[p.first] = transformation * p.second, Iterators.filter(filterpredicate, fixedmode.attrs))
        newmode::Mode = Mode(newattrs) |> correctsymmetrizedmode
        modemapping[newmode] = mode
        return newmode
    end

    connections::Dict{Tuple{Mode, Mode}, ComplexF64} = Dict()

    function rebaseorbital(mode::Mode)::Mode
        frommode::Mode = modemapping[mode]
        tomode::Mode = mode |> setattr(:orbital => (frommode |> getorbital()))
        connections[(tomode, frommode)] = (mode |> getorbital()) |> relativephase(frommode |> getorbital())
        return tomode
    end

    outmodes::Subset{Mode} = Subset(mode |> modesymmetrize |> rebaseorbital for mode in regionfock)
    
    return FockMap(outmodes |> FockSpace, regionfock, connections)
end

function Base.:*(transformation::AffineTransform, crystalfock::FockSpace{Crystal})::FockMap
    homefock::FockSpace = crystalfock |> unitcellfock
    homefocktransform::FockMap = transformation * homefock
    ksubspaces::Dict{Momentum, FockSpace} = crystalfock |> crystalsubspaces
    fouriertransform::FockMap = fourier(crystalfock, homefock)
    transformedfourier::FockMap = fourier(crystalfock, homefocktransform.outspace)
    transform::FockMap = directsum(
        rows(transformedfourier, ksubspaces[transformation * k |> basispoint]) * homefocktransform * rows(fouriertransform, fockspace)'
        for (k, fockspace) in ksubspaces)
    crystal::Crystal = crystalfock |> getcrystal
    return FockMap(transform, outspace=FockSpace(transform.outspace, reflected=crystal), inspace=FockSpace(transform.inspace, reflected=crystal))
end

""" Defines the computation of the `transform` representation in the `inspace` with the full data provided by the data and `outspace` of the fockmap. """
function Base.:*(fockmap::FockMap, transform::AffineTransform)::FockMap
    outspacerep::FockMap = *(transform, fockmap |> getoutspace)
    hassamespan(outspacerep |> getoutspace, fockmap |> getoutspace) || error("The symmetry action on the outspace in not closed!")
    return fockmap' * outspacerep * fockmap
end

""" Defines the computation of the `transform` representation in the `outspace` with the full data provided by the data and `inspace` of the fockmap. """
function Base.:*(transform::AffineSpace, fockmap::FockMap)::FockMap
    inspacerep::FockMap = *(transform, fockmap |> getinspace)
    hassamespan(inspacerep |> getoutspace, fockmap |> getinspace) || error("The symmetry action on the inspace in not closed!")
    return fockmap' * inspacerep * fockmap
end

"""
    spatialmap(fockmap::FockMap)::FockMap

Given `fockmap` with a `outspace` of `FockMap{Region}`, determine the center position of the column function and generate a identity map that transforms
the `inspace` of `fockmap` to include the actual physical attribute of `:offset` and `:pos`.

### Output
The transformer `FockMap` with `inspace` of `fockmap` and the spatially decorated version as the `outspace`.
"""
function spatialmap(fockmap::FockMap)::FockMap
    function spatialinmode(colmap::FockMap)
        inmode::Mode = colmap |> getinspace |> first
        absmap::FockMap = colmap |> abs
        weights::FockMap = absmap / (absmap |> rep |> collect |> real |> sum)
        modecenter::Offset = reduce(+, (outmode |> getpos) * (weights[outmode, inmode] |> real) for outmode in weights |> getoutspace)
        basis::Offset = modecenter |> basispoint
        offset::Offset = modecenter - basis
        return inmode |> setattr(:offset => offset) |> setattr(:pos => basis)
    end

    spatialinspace::FockSpace{Region} = FockSpace{Region}(fockmap[:, m] |> spatialinmode for m in fockmap |> getinspace)
    return idmap(spatialinspace, fockmap |> getinspace)
end
export spatialmap

function getsymmetrizer(symmetry::AffineTransform, fockmap::FockMap)::FockMap
    inspacerep::FockMap = getinspacerep(symmetry, fockmap)
    phasespectrum::EigenSpectrum = inspacerep |> eigspec
    inspace::FockSpace = FockSpace(
        m |> setattr(:orbital => findeigenfunction(symmetry, eigenvalue=(phasespectrum |> geteigenvalues)[m]))
          |> removeattr(:eigenindex) # The :orbital can subsitute the :eigenindex.
        for m in phasespectrum |> geteigenvectors |> getinspace)
    return FockMap(phasespectrum |> geteigenvectors, inspace=inspace, performpermute=false)
end
export getsymmetrizer
