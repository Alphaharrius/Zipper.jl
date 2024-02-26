function extendedcrystalscale(; crystal::Crystal, generatingvector::Offset)::Scale
    snfinput::Matrix{Integer} = map(Integer, vcat(generatingvector|>vec|>transpose, crystal|>size|>diagm)) 
    _, S, Vd = snfinput|>dosnf
    return Scale((S|>diag|>diagm) * Vd |>transpose, crystal|>getspace)
end
export extendedcrystalscale

struct ExtendedRestrict <: Transformation{Matrix{Float64}}
    scale::Scale
    normalvector::Offset
    stripradius::Real
end
export ExtendedRestrict

Base.:convert(::Type{Matrix{Float64}}, restrict::ExtendedRestrict) = restrict.scale|>rep
Zipper.getspace(restrict::ExtendedRestrict) = restrict.scale|>getspace

function extendedcrystalrestrict(; crystal::Crystal, normalvector::Offset, stripradius::Real)::ExtendedRestrict
    snfinput::Matrix{Integer} = map(Integer, vcat(normalvector|>vec|>transpose, crystal|>size|>diagm)) 
    _, S, Vd = snfinput|>dosnf
    scale = Scale((S|>diag|>diagm) * Vd |>transpose, crystal|>getspace)
    return ExtendedRestrict(scale, normalvector, stripradius)
end
export extendedcrystalrestrict

function Base.:*(restrict::ExtendedRestrict, crystalfock::CrystalFock)
    crystal::Crystal = crystalfock|>getcrystal
    scaledcrystal::Crystal = (restrict.scale) * crystal
    scaledspace::RealSpace = scaledcrystal|>getspace
    scaledkspace::MomentumSpace = convert(MomentumSpace, scaledspace)
    # This is the unit cell of the 1-D embedded crystal defining the strip region in the blocked crystal metric space.
    unscaledextendedunitcell::Region = getcrosssection(crystal=crystal, normalvector=restrict.normalvector, radius=restrict.stripradius)
    unscaledextendedhomefock = matchreduce(
        quantize(unscaledextendedunitcell, 1), crystalfock|>unitcellfock,
        leftmatch=m->m|>getattr(:b)|>euclidean,
        rightmatch=m->m|>getattr(:b)|>euclidean,
        reducer=(a, b)->merge(a, b))|>RegionFock

    bz::Subset{Momentum} = crystal|>brillouinzone
    momentummappings::Base.Generator = (basispoint(scaledkspace * k) => k for k in bz)

    restrictedfourier = fourier(crystalfock, unscaledextendedhomefock)'
    blocks::Dict = Dict()

    stripunitcell::Region = Subset(scaledspace * r for r in unscaledextendedunitcell)
    stripcrystal::Crystal = Crystal(stripunitcell, scaledcrystal|>size)
    volumeratio::Real = vol(crystal) / vol(stripcrystal)

    scaledksubspaces::Dict{Momentum, FockSpace} = Dict()
    for (scaledk, k) in momentummappings
        kfourier::FockMap = restrictedfourier[:, k] / sqrt(volumeratio)
        if !haskey(scaledksubspaces, scaledk)
            scaledksubspaces[scaledk] = FockSpace(
                setattr(mode, :k=>scaledk, :b=>(scaledspace * getpos(mode)))|>removeattr(:r) for mode in kfourier|>getoutspace)
        end
        blocks[(scaledk, k)] = FockMap(scaledksubspaces[scaledk], kfourier|>getinspace, kfourier|>rep)
    end

    return CrystalFockMap(stripcrystal, crystal, blocks)
end
