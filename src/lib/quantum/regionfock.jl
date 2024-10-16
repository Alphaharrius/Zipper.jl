# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ RegionFock APIs ◆
# RegionFock is an extension of NormalFock which provide 
# region information via the reflected attribute.

"""
    RegionFock(input)

From the input modes, reformat the mode format to have both `:r` and `:b` 
attribute and output the set as a `RegionFock`.

### Error
- If the input mode already has a `:k` attribute, an error will be thrown.
"""
function Zipper.RegionFock(input)
    function preprocess(mode::Mode)::Mode
        if hasattr(mode, :k)
            error("The mode has a momentum space position attribute!")
        end
        if !hasattr(mode, :r)
            return mode|>setattr(:r=>mode|>getattr(:b)|>getspace|>getorigin)
        end
        return mode
    end

    processedinput = (m|>preprocess for m in input)
    region::Region = Subset(m|>getpos for m in processedinput)
    return FockSpace(processedinput, reflected=region)
end

""" Get the reflected region of the regional `FockSpace`. """
Zipper.:getregion(regionfock::RegionFock)::Region = regionfock |> getreflected

"""
Translate a `RegionFock` by a `Offset` and recompute the 
corresponding `:r` and `:b` attributes.
"""
function Base.:+(regionfock::RegionFock, offset::Offset)::RegionFock
    positiontomodes = ((getpos(mode) + offset)=>mode for mode in regionfock)
    return RegionFock(
        mode|>setattr(:r=>(pos - basispoint(pos)), :b=>basispoint(pos)) 
        for (pos, mode) in positiontomodes)
end

Base.:-(regionfock::RegionFock, offset::Offset)::RegionFock = regionfock + (-offset)

"""
    snap2unitcell(regionfock::RegionFock)

Snap the `RegionFock` to the unit cell by forcing all `:b` to the positive 
parallelogram spanned by the spatial basis vectors, and relocate the differences 
to the `:r` attribute.
"""
function snap2unitcell(regionfock::RegionFock)
    spatials = ((mode|>getpos, mode) for mode in regionfock)
    spatials = ((r, r|>basispoint, mode) for (r, mode) in spatials)
    spatials = ((r-b, b, mode) for (r, b, mode) in spatials)
    return (
        m|>setattr(:r=>r, :b=>b) 
        for (r, b, m) in spatials)|>mapmodes(m->m)|>RegionFock
end
export snap2unitcell
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ RegionFock implementation of FockSpace interfaces ◆
unitcellfock(regionfock::RegionFock) = regionfock|>removeattr(:r)|>FockSpace
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ RegionFock display extension ◆
Base.:show(io::IO, ::Type{RegionFock}) = print(io, string("RegionFock"))
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
