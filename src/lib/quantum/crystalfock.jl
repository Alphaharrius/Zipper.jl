# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalFock definition ◆
struct CrystalFock <: FockSpace
    crystal::Crystal
    korderings::Dict{Momentum, Integer}
    homefock::NormalFock
end
export CrystalFock
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Implementation of FockSpace APIs ◆
orderedmodes(crystalfock::CrystalFock) = (
    # brillouinzone is the default ordering of the momentums since korderings is generated from it too.
    mode|>setattr(:k=>k) for k in crystalfock|>getcrystal|>brillouinzone for mode in crystalfock|>unitcellfock)

getmodes(crystalfock::CrystalFock) = Set(crystalfock|>orderedmodes)

"""
    unitcellfock(crystalfock::CrystalFock)::FockSpace

Retrieve the unit cell fockspace of the system from a `CrystalFock`, 
positioned at the origin of the parent `AffineSpace`.
"""
unitcellfock(crystalfock::CrystalFock) = crystalfock.homefock
export unitcellfock

"""
    subspacecount(fockspace::FockSpace)::Integer

Get the number of sub-fockspaces of this `fockspace`.
"""
subspacecount(crystalfock::CrystalFock)::Integer = vol(crystalfock.crystal)
export subspacecount

function Base.:getindex(crystalfock::CrystalFock, mode::Mode)
    k = mode|>getattr(:k)
    kordering = crystalfock.korderings[k]
    homefock = crystalfock|>unitcellfock
    return (kordering - 1) * dimension(homefock) + (crystalfock|>unitcellfock)[mode|>removeattr(:k)]
end

Base.:getindex(crystalfock::CrystalFock, range::UnitRange) = (crystalfock|>collect)[range]

""" Check whether a `Mode` is in the `FockSpace`. """
Base.:in(mode::Mode, crystalfock::CrystalFock)::Bool = (
    haskey(mode|>getattr(:k), crystalfock.korderings) && mode|>removeattr(:k) in crystalfock|>unitcellfock)
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Implementation of spatial/geometric APIs ◆
"""
    dimension(fockspace::FockSpace)

Returns the number of unique member modes within the `fockspace`, each of those represents a vector 
from the Hilbert space.
"""
Zipper.:dimension(crystalfock::CrystalFock) = vol(crystalfock.crystal) * dimension(crystalfock.homefock)

"""
    getcrystal(crystalfock::CrystalFock)::Crystal

Shorthand for retrieving the `Crystal` of a `CrystalFock`.
"""
Zipper.:getcrystal(crystalfock::CrystalFock)::Crystal = crystalfock.crystal
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ CrystalFock APIs ◆
"""
    getcrystalfock(homefock::NormalFock, crystal::Crystal)

A short hand to build the crystal fockspace, which is the fockspace containing all modes spanned from 
`homefock` by the brillouin zone of the `crystal`.
"""
@memoize function getcrystalfock(homefock::NormalFock, crystal::Crystal)
    korderings::Dict{Momentum, Integer} = Dict(k=>n for (n, k) in crystal|>brillouinzone|>enumerate)
    return CrystalFock(crystal, korderings, homefock)
end

getcrystalfock(homemodes::Subset{Mode}, crystal::Crystal) = getcrystalfock(FockSpace(homemodes), crystal)

export getcrystalfock

"""
    getmomentum(subspace::MomentumFock)::Momentum

Shorthand for retrieving the `Momentum` of a `MomentumFock`.
"""
getmomentum(subspace::MomentumFock)::Momentum = subspace.reflected
export getmomentum

"""
    getsubspace(crystalfock::CrystalFock, k::Momentum)::MomentumFock

Retrieve the momentum subspace of the `CrystalFock` at the given `k` momentum.
"""
getsubspace(crystalfock::CrystalFock, k::Momentum)::MomentumFock = FockSpace(
    crystalfock.homefock|>setattr(:k=>k), reflected=k)
export getsubspace

"""
    crystalsubsets(crystalfock::CrystalFock)::Dict{Momentum, Subset{Mode}}

    Retrieve mappings from the crystal momentums to the corresponding `Subset{Mode}`.
"""
crystalsubsets(crystalfock::CrystalFock)::Dict{Momentum, Subset{Mode}} = Dict(
    k=>crystalfock.homefock|>setattr(:k=>k) for (k, _) in crystalfock.korderings)
export crystalsubsets

"""
    crystalsubspaces(crystalfock::CrystalFock)::Base.Generator

Retrieve mappings from the crystal momentums to the corresponding fockspaces.
"""
crystalsubspaces(crystalfock::CrystalFock)::Base.Generator = (
    k=>FockSpace(
        crystalfock.homefock|>setattr(:k=>k), reflected=k) for (k, _) in crystalfock.korderings)
export crystalsubspaces

""" Check if two `CrystalFock` objects have the same span. """
function Zipper.:hassamespan(a::CrystalFock, b::CrystalFock)::Bool
    return getcrystal(a) == getcrystal(b) && hassamespan(a|>unitcellfock, b|>unitcellfock)
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
