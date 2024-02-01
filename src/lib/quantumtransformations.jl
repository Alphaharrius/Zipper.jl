function Base.:*(scale::Scale, crystalfock::FockSpace{Crystal})::FockMap
    crystal::Crystal = crystalfock |> getcrystal
    scaledcrystal::Crystal = scale * crystal
    unscaledblockedregion::Subset{Offset} = (scale |> inv) * scaledcrystal.unitcell
    bz::Subset{Momentum} = crystal |> brillouinzone
    basismodes::Subset{Mode} = crystalfock|>unitcellfock|>orderedmodes

    momentummappings::Base.Generator = (basispoint(scale * p) => p for p in bz)

    unscaledblockedlatticeoffsets::Subset{Offset} = Subset(pbc(crystal, p) |> latticeoff for p in unscaledblockedregion)
    unscaledblockedunitcellfock::FockSpace = spanoffset(basismodes, unscaledblockedlatticeoffsets)

    volumeratio::Number = (crystal |> vol) / (scaledcrystal |> vol)

    crystalfocksubspaces::Dict{Momentum, FockSpace} = crystalfock|>crystalsubspaces|>Dict
    restrictedfourier::FockMap = fourier(crystalfock, unscaledblockedunitcellfock|>RegionFock)'
    blocks::Dict = Dict()

    scaledksubspaces::Dict{Momentum, FockSpace} = Dict()
    for (scaledk, k) in momentummappings
        kfourier::FockMap = columns(restrictedfourier, crystalfocksubspaces[k]) / sqrt(volumeratio)
        if !haskey(scaledksubspaces, scaledk)
            scaledksubspaces[scaledk] = FockSpace(
                setattr(mode, :k=>scaledk, :b=>(scale * getpos(mode)))|>removeattr(:r) for mode in kfourier|>getoutspace)
        end
        blocks[(scaledk, k)] = FockMap(scaledksubspaces[scaledk], kfourier|>getinspace, kfourier|>rep)
    end

    return CrystalFockMap(scaledcrystal, crystal, blocks)
end

function Base.:*(transformation::AffineTransform, regionfock::FockSpace{Region})::FockMap
    # This is used to correct the :b attribute, since the :b as a Point will be symmetrized,
    # which the basis point set might not include the symmetrized :b. Thus we would like to set
    # the :b to its corresponding basis point, and offload the difference to :r.
    function correctsymmetrizedmode(mode::Mode)::Mode
        actualposition::Offset = mode |> getattr(:b)
        basisposition::Offset = actualposition |> basispoint
        offset::Offset = actualposition - basisposition
        return mode |> setattr(:b => basisposition) |> setattr(:r => offset)
    end

    modemapping::Dict{Mode, Mode} = Dict()

    function mergepositions(mode::Mode)::Mode
        latticeoffset::Point = mode |> getattr(:r)
        latticeoffset isa Offset || error("Transforming a mode based on momentum must be done with crystal fockspace!")
        actualposition::Offset = getattr(mode, :b) + latticeoffset
        return mode |> setattr(:b => actualposition) |> removeattr(:r)
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
    fouriertransform::FockMap = fourier(crystalfock, homefock|>RegionFock)
    transformedfourier::FockMap = fourier(crystalfock, homefocktransform|>getoutspace|>RegionFock)
    transform::FockMap = directsum(
        rows(transformedfourier, ksubspaces[transformation * k |> basispoint]) * homefocktransform * rows(fouriertransform, fockspace)'
        for (k, fockspace) in ksubspaces)
    crystal::Crystal = crystalfock |> getcrystal
    return FockMap(transform, outspace=FockSpace(transform.outspace, reflected=crystal), inspace=FockSpace(transform.inspace, reflected=crystal), performpermute=false)
end

function Base.:*(symmetry::AffineTransform, state::RegionState)
    statemap::FockMap = state|>FockMap
    stateregionfock::RegionFock = statemap|>getoutspace
    localsymmetry::AffineTransform = symmetry|>recenter(stateregionfock|>getregion|>getcenter)
    outspacerep::FockMap = localsymmetry * stateregionfock
    hassamespan(outspacerep|>getoutspace, stateregionfock) || error("The symmetry action on the state region fockspace in not closed!")
    inspacerep::FockMap = statemap' * outspacerep * statemap
    phasespectrum::EigenSpectrum = inspacerep|>eigspec

    function mapper(mode::Mode)::Mode
        eigenvalue::Complex = (phasespectrum|>geteigenvalues)[mode]
        basisfunction::BasisFunction = findeigenfunction(localsymmetry, eigenvalue=eigenvalue)
        # The existing :flavor attribute is not needed since we will use mapmodes to determine the new one;
        # the :eigenindex is not needed since it does not represent any physical attributes.
        return mode|>setattr(:orbital=>basisfunction)|>removeattr(:eigenindex, :flavor)
    end

    symmetricfock::FockSpace = phasespectrum|>geteigenvectors|>getinspace|>mapmodes(mapper)|>FockSpace
    symmetrizer::FockMap = FockMap(phasespectrum|>geteigenvectors, inspace=symmetricfock, performpermute=false)

    return RegionState(statemap * symmetrizer)
end

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

function Base.:*(fockmap::FockMap, symmetry::AffineTransform)
    outspacerep::FockMap = *(symmetry, fockmap |> getoutspace)
    hassamespan(outspacerep |> getoutspace, fockmap |> getoutspace) || error("The symmetry action on the outspace in not closed!")
    inspacerep::FockMap = fockmap' * outspacerep * fockmap
    phasespectrum::EigenSpectrum = inspacerep |> eigspec
    # This is used to correct the :flavor attribute when there are multiple mode with the same orbital is found.
    processedmodes::Dict{Mode, Integer} = Dict()
    newmodes::Vector{Mode} = []
    # Previously this is using a function to generate the new mode, but it is not working with the FockSpace constructor 
    # when creating the inspace since it will call the constructor of Subset(::Base.Generator), which it uses Base.first 
    # to determine the type of the generator elements, which will invoke the function and cause the first mode to be processed 
    # once before the main iteration. This will cause one of the orbital group start with index 2 instead of 1 since it 
    # has been processed once and loaded into the dictionary processedmodes.
    for mode in phasespectrum|>geteigenvectors|>getinspace
        eigenvalue::Complex = (phasespectrum|>geteigenvalues)[mode]
        basisfunction::BasisFunction = findeigenfunction(symmetry, eigenvalue=eigenvalue)
        newmode::Mode = mode|>setattr(:orbital=>basisfunction)|>removeattr(:eigenindex, :flavor)
        if !haskey(processedmodes, newmode)
            processedmodes[newmode] = 1
        else
            processedmodes[newmode] += 1
        end
        push!(newmodes, newmode|>setattr(:flavor=>processedmodes[newmode]))
    end

    inspace::FockSpace = FockSpace(newmodes)
    return FockMap(phasespectrum |> geteigenvectors, inspace=inspace, performpermute=false)
end
