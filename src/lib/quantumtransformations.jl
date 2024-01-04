function Base.:*(scale::Scale, crystalfock::FockSpace{Crystal})::FockMap
    crystal::Crystal = crystalfock |> getcrystal
    scaledcrystal::Crystal = scale * crystal
    unscaledblockedregion::Subset{Offset} = (scale |> inv) * scaledcrystal.unitcell
    bz::Subset{Momentum} = crystal |> brillouinzone
    basismodes::Subset{Mode} = crystalfock |> rep |> first

    momentummappings::Base.Generator = (basispoint(scale * p) => p for p in bz)

    unscaledblockedlatticeoffsets::Subset{Offset} = Subset(pbc(crystal, p) |> latticeoff for p in unscaledblockedregion)
    unscaledblockedunitcellfock::FockSpace = spanoffset(basismodes, unscaledblockedlatticeoffsets)

    volumeratio::Number = (crystal |> vol) / (scaledcrystal |> vol)

    crystalfocksubspaces::Dict{Momentum, FockSpace} = crystalfock|>crystalsubspaces|>Dict
    restrictedfourier::FockMap = fourier(bz, unscaledblockedunitcellfock)'
    blocks::Dict = Dict()

    scaledksubspaces::Dict{Momentum, FockSpace} = Dict()
    for (scaledk, k) in momentummappings
        kfourier::FockMap = columns(restrictedfourier, crystalfocksubspaces[k]) / sqrt(volumeratio)
        if !haskey(scaledksubspaces, scaledk)
            scaledksubspaces[scaledk] = FockSpace(
                setattr(mode, :offset=>scaledk, :pos=>(scale * getpos(mode))) for mode in kfourier|>getoutspace)
        end
        blocks[(scaledk, k)] = FockMap(scaledksubspaces[scaledk], kfourier|>getinspace, kfourier|>rep)
    end

    return CrystalFockMap(scaledcrystal, crystal, blocks)
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
          |> removeattr(:eigenindex) # The :orbital can subsitute the :eigenindex.
        for m in phasespectrum |> geteigenvectors |> getinspace)
    return FockMap(phasespectrum |> geteigenvectors, inspace=inspace, performpermute=false)
end
