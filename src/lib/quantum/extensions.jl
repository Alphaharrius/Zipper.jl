# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ extensions.jl ◆

# This file contains the extensions APIs that integrates cross-type operations within the folder /quantum.
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Geometrical extension APIs ◆
"""
    weightedspatialmap(localisometry::FockMap)::FockMap{RegionFock, RegionFock}

Given `localisometry` with an `outspace` of `FockMap{Region}`, determine the center position of the column 
function by the weights of the eigenfunction and generate a identity map that transforms the `inspace` of 
`localisometry` to include the actual physical attribute of `:r` and `:b`.

### Output
The transformer `FockMap` with `inspace` of `localisometry` and the spatially decorated version as the `outspace`.
"""
function weightedspatialmap(localisometry::FockMap)::FockMap{RegionFock, RegionFock}
    function mapper(mode::Mode)::Mode
        modeisometry::FockMap = localisometry[:, mode]
        absmap::FockMap = modeisometry|>abs
        weights::FockMap = absmap / (absmap|>rep|>collect|>real|>sum)
        pos::Offset = reduce(+, (outmode|>getpos) * (weights[outmode, mode]|>real) for outmode in weights|>getoutspace)
        b::Offset = pos|>basispoint
        r::Offset = pos - b
        # removeattr(:eigenindex) is added since localisometry is coming from a EigenSpectrum.
        return mode|>setattr(:r=>r, :b=>b)|>removeattr(:eigenindex)
    end
    return idmap(localisometry|>getinspace, localisometry|>getinspace|>mapmodes(mapper)|>RegionFock)
end
export weightedspatialmap

function spatialmap(localisometry::FockMap)::FockMap
    function mapper(mode::Mode)::Mode
        pos::Offset = localisometry|>getoutspace|>getregion|>getcenter
        b::Offset = pos|>basispoint
        r::Offset = pos - b
        # removeattr(:eigenindex) is added since localisometry is coming from a EigenSpectrum.
        return mode|>setattr(:r=>r, :b=>b)|>removeattr(:eigenindex)
    end
    return idmap(localisometry|>getinspace, localisometry|>getinspace|>mapmodes(mapper)|>RegionFock)
end
export spatialmap

"""
    spatialremapper(regionfock::RegionFock; offsets::Region, unitcell::Region)

This function remaps the positions of modes in a `RegionFock` object based on provided offsets and a unit cell.
This function is used when the unit-cell defined in the original `regionfock` does not match with what you desire, 
most of the case the `unitcell` is not defined in the positive parallelogram spanned by the lattice basis vectors.

### Input
- `regionfock`  The `RegionFock` object to be remapped.

### Output
It creates a new `RegionFock` object with the remapped modes and returns an identity map between the new and original 
`RegionFock` objects.
"""
function spatialremapper(regionfock::RegionFock; offsets::Region, unitcell::Region)
    positions::Dict{Offset, Tuple} = Dict((r + b)=>(r, b) for (r, b) in Iterators.product(offsets, unitcell))
    remappingdata::Base.Generator = ((positions[mode|>getpos], mode) for mode in regionfock)
    remappedfock::RegionFock = RegionFock(mode|>setattr(:r=>r, :b=>b) for ((r, b), mode) in remappingdata)
    return idmap(remappedfock, regionfock)
end
export spatialremapper
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
