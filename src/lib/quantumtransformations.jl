function Base.:*(scale::Scale, crystalfock::FockSpace{Crystal})::FockMap
    crystal::Crystal = crystalfock |> getcrystal
    scaledcrystal::Crystal = scale * crystal
    unscaledblockedregion::Subset{Offset} = (scale |> inv) * scaledcrystal.unitcell
    bz::Subset{Momentum} = crystal |> brillouinzone
    basismodes::Subset{Mode} = crystalfock|>unitcellfock|>orderedmodes
    scaledbz::Subset{Momentum} = scaledcrystal |> brillouinzone

    momentummappings::Base.Generator = (basispoint(scale * p) => p for p in bz)
    mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point, Vector{Point}}()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end

    ksubsets::Dict{Momentum, Subset{Mode}} = crystalfock |> crystalsubsets
    scaledfock::FockSpace = ((ksubsets[k] for k in mappingpartitions[kscaled]) |> subsetunion |> FockSpace
        for kscaled in scaledbz) |> fockspaceunion

    # TODO: Use of latticeoff requires revisit, since latticeoff only truncates the decimals.
    unscaledblockedlatticeoffsets::Region = Subset(pbc(crystal, p) |> latticeoff for p in unscaledblockedregion)
    unscaledblockedunitcellfock::RegionFock = spanoffset(basismodes, unscaledblockedlatticeoffsets)|>RegionFock

    restrictedfourier::FockMap = fourier(crystalfock, unscaledblockedunitcellfock)
    volumeratio::Number = (crystal |> vol) / (scaledcrystal |> vol)
    permutedfourier::FockMap = Zipper.permute(restrictedfourier, outspace=scaledfock) / sqrt(volumeratio)

    function repackfourierblocks(source::FockMap, kscaled::Momentum, partition::Subset{Mode})::FockMap
        partitionrows::FockMap = rows(source, partition |> FockSpace)
        inspace::FockSpace = FockSpace(mode|>setattr(:k=>kscaled, :b=>*(scale, mode|>getpos))|>removeattr(:r) for mode in partitionrows|>getinspace)
        return FockMap(partitionrows.outspace, inspace, partitionrows |> rep)
    end

    repackedblocks::Base.Generator = (
        repackfourierblocks(permutedfourier, kscaled, partition)
        for (kscaled, partition) in Iterators.zip(scaledbz, scaledfock |> rep))
    blocking::FockMap = directsum(repackedblocks)
    return FockMap(blocking, inspace=FockSpace(blocking |> getinspace, reflected=scaledcrystal), outspace=crystalfock)'
end

function Base.:*(transformation::AffineTransform, regionfock::FockSpace{Region})::FockMap
    # This is used to correct the :b attribute, since the :b as a Point will be symmetrized,
    # which the basis point set might not include the symmetrized :b. Thus we would like to set
    # the :b to its corresponding basis point, and offload the difference to :r.
    function correctsymmetrizedmode(mode::Mode)::Mode
        actualposition::Offset = mode|>getattr(:R)
        basisposition::Offset = actualposition|>basispoint
        offset::Offset = actualposition - basisposition
        return mode|>setattr(:b=>basisposition)|>setattr(:r=>offset)|>removeattr(:R)
    end

    modemapping::Dict{Mode, Mode} = Dict()

    mergepositions(mode::Mode)::Mode = mode|>setattr(:R=>getpos(mode))|>removeattr(:r, :b)

    function modesymmetrize(mode::Mode)::Mode
        fixedmode = mode|>mergepositions
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
    
    return FockMap(outmodes |> FockSpace{Region}, regionfock, connections)
end

function Base.:*(transformation::AffineTransform, crystalfock::FockSpace{Crystal})::FockMap
    homefock::FockSpace = crystalfock |> unitcellfock
    homefocktransform::FockMap = transformation * homefock
    ksubspaces::Dict{Momentum, FockSpace} = crystalfock |> crystalsubspaces |> Dict
    fouriertransform::FockMap = fourier(crystalfock, homefock)
    transformedfourier::FockMap = fourier(crystalfock, homefocktransform.outspace)
    transform::FockMap = directsum(
        rows(transformedfourier, ksubspaces[transformation * k |> basispoint]) * homefocktransform * rows(fouriertransform, fockspace)'
        for (k, fockspace) in ksubspaces)
    crystal::Crystal = crystalfock |> getcrystal
    return FockMap(transform, outspace=FockSpace(transform.outspace, reflected=crystal), inspace=FockSpace(transform.inspace, reflected=crystal), performpermute=false)
end

"""
    spatialmap(fockmap::FockMap)::FockMap

Given `fockmap` with a `outspace` of `FockMap{Region}`, determine the center position of the column function and generate a identity map that transforms
the `inspace` of `fockmap` to include the actual physical attribute of `:r` and `:b`.

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
        return inmode|>setattr(:r=>offset)|>setattr(:b=>basis)
    end

    spatialinspace::FockSpace{Region} = FockSpace{Region}(fockmap[:, m] |> spatialinmode for m in fockmap |> getinspace)
    return idmap(spatialinspace, fockmap |> getinspace)
end
export spatialmap

function Base.:*(symmetry::AffineTransform, fockmap::FockMap)::FockMap
    inspacerep::FockMap = *(symmetry, fockmap |> getinspace)
    hassamespan(inspacerep |> getoutspace, fockmap |> getinspace) || error("The symmetry action on the inspace in not closed!")
    outspacerep::FockMap = fockmap' * inspacerep * fockmap

    phasespectrum::EigenSpectrum = outspacerep |> eigspec
    outspace::FockSpace = FockSpace(
        m |> setattr(:orbital => findeigenfunction(symmetry, eigenvalue=(phasespectrum |> geteigenvalues)[m]))
          |> removeattr(:eigenindex) # The :orbital can subsitute the :eigenindex.
        for m in phasespectrum |> geteigenvectors |> getinspace)s
    return FockMap(phasespectrum |> geteigenvectors, inspace=outspace, performpermute=false)'
end

function Base.:*(fockmap::FockMap, symmetry::AffineTransform)::FockMap
    outspacerep::FockMap = *(symmetry, fockmap |> getoutspace)
    hassamespan(outspacerep |> getoutspace, fockmap |> getoutspace) || error("The symmetry action on the outspace in not closed!")
    inspacerep::FockMap = fockmap' * outspacerep * fockmap

    phasespectrum::EigenSpectrum = inspacerep |> eigspec
    inspace::FockSpace = FockSpace(
        m |> setattr(:orbital => findeigenfunction(symmetry, eigenvalue=(phasespectrum |> geteigenvalues)[m]))
        # TODO: Requires revisit.  
        # |> removeattr(:eigenindex) # The :orbital can subsitute the :eigenindex.
        for m in phasespectrum |> geteigenvectors |> getinspace)
    return FockMap(phasespectrum |> geteigenvectors, inspace=inspace, performpermute=false)
end
