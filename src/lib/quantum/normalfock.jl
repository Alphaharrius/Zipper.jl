# ================================================================================================
# NormalFock implementation
struct NormalFock{T} <: Zipper.FockSpace
    reflected::T
    subset::Subset
end
export NormalFock

""" Create an empty fockspace. """
NormalFock() = NormalFock(nothing, Subset{Mode}())
"""
    NormalFock(subset::Subset{Mode}; reflected=nothing)

Create a `FockSpace` with the given `subset` of `Mode` which follows the ordering of the `subset`.

### Input
- `subset` The subset of `Mode` to use.
- `reflected` The physical object this `FockSpace` is reflected to, defaults to `nothing`.
"""
NormalFock(subset::Subset{Mode}; reflected=nothing) = NormalFock(reflected, subset)
""" Create a `FockSpace` from an iterator of `Mode`. """
NormalFock(iter; reflected=nothing) = NormalFock(reflected, Subset(iter))
""" Adding a reflected object to a `FockSpace`."""
NormalFock(fock::NormalFock; reflected=nothing) = NormalFock(reflected, fock|>rep)
""" Create a `FockSpace` with a single `Mode`. """
NormalFock(mode::Mode; reflected=nothing) = NormalFock(reflected, Subset(mode))
# ================================================================================================

# ================================================================
# NormalFock logicals
""" Check whether a `Mode` is in the `FockSpace`. """
Base.:in(mode::Mode, fock::NormalFock)::Bool = in(mode, fock|>rep)
# ================================================================

# =============================================
# NormalFock APIs
getreflected(fock::NormalFock) = fock.reflected
export getreflected
# =============================================

# ===============================================================================================================
# NormalFock essentials
Base.:eltype(::NormalFock) = Mode
# We will not include a hash function for NormalFock since it is handled by the default hash function of Element.
# ===============================================================================================================

# =========================================================================================================
# FockSpace interface implementations
Base.:convert(::Type{Subset}, fock::NormalFock) = fock.subset

""" Get the dimension of the `FockSpace` which is essentially the number of `Mode` objects it contains. """
Zipper.:dimension(fock::NormalFock) = fock|>rep|>length

"""
    orderedmodes(fock::NormalFock)

Get the `Mode` objects sequentially according to their order in the `FockSpace`.
"""
orderedmodes(fock::NormalFock) = fock|>rep

# We don't have to include other iterator methods like `iterate`, `length` and `lastindex` since they are 
# supported via the definition of `orderedmodes` and `dimension`.

"""
    getmodes(fockspace::FockSpace)::Set{Mode}

Returns an unordered set of modes of `fockspace`, this is a more efficient way 
to retrieve the underlying modes for a fockspace with more than one partitions.
"""
getmodes(fock::NormalFock)::Set{Mode} = Set(mode for mode in fock)

""" There are no concrete subspace in `NormalFock`, but it is a subspace itself. """
subspacecount(fock::NormalFock) = 1

"""
This have multiple implications, if the input is an Integer, it will return the `Mode` at 
the `i`th order of the `FockSpace`; if the input is a `Mode`, it will return the order of 
the `Mode` in the `FockSpace`.
"""
Base.:getindex(fock::NormalFock, v) = (fock|>rep)[v]
""" Get a slice of the `FockSpace`, this will not preserve the `reflected` attribute. """
Base.:getindex(fock::NormalFock, range::UnitRange) = NormalFock((fock|>rep)[range])
# =======================================================================================

# =================================================================================================
# NormalFock display
""" Displays the fock type, subspace count and dimension information of a `FockSpace`. """
Base.:show(io::IO, fock::NormalFock) = print(io, string("$(fock|>typeof)(dim=$(fock|>dimension))"))
# =================================================================================================

# =================================================
# NormalFock subtypes
""" Shorthand alias for `NormalFock{Region}`. """
RegionFock = NormalFock{Region}
export RegionFock

""" Shorthand alias for `NormalFock{Momentum}`. """
MomentumFock = NormalFock{Momentum}
export MomentumFock

""" Shorthand alias for `NormalFock{Offset}`. """
SiteFock = NormalFock{Offset}
export SiteFock
# =================================================

# ==========================================================================================
# FockSpace default implementation

# Since the name FockSpace is used in many legacy locations, 
# we have to maintain backward compatabilities.
FockSpace(subset::Subset{Mode}; reflected=nothing) = NormalFock(subset, reflected=reflected)
FockSpace(iter; reflected=nothing) = NormalFock(iter, reflected=reflected)
FockSpace(fock::NormalFock; reflected=nothing) = NormalFock(fock, reflected=reflected)
FockSpace(mode::Mode; reflected=nothing) = NormalFock(mode, reflected=reflected)
# ==========================================================================================
