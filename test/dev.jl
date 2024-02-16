using Base.Iterators
using Zipper

square = euclidean(RealSpace, 2)
point = [1/2, 1/2] ∈ square
spatialsnappingcalibration([point])

kspace = convert(MomentumSpace, square)

c4 = pointgrouptransform([0 -1; 1 0])

unitcell = Subset(point)
crystal = Crystal(unitcell, [384, 384])
reciprocalhashcalibration(crystal.sizes)

bz = brillouinzone(crystal)

# ================================================================================================
# NormalFock implementation
struct NormalFock{T} <: Zipper.FockSpace
    reflected::T
    subset::Subset
end

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

# ===============================================================================================================
# NormalFock essentials
Base.:convert(::Type{Subset}, fock::NormalFock) = fock.subset
Base.:eltype(::NormalFock) = Mode
# We will not include a hash function for NormalFock since it is handled by the default hash function of Element.
# ===============================================================================================================

# =========================================================================================================
# NormalFock iterator methods
""" Get the dimension of the `FockSpace` which is essentially the number of `Mode` objects it contains. """
Zipper.:dimension(fock::NormalFock) = fock|>rep|>length

"""
    orderedmodes(fock::NormalFock)

Get the `Mode` objects sequentially according to their order in the `FockSpace`.
"""
orderedmodes(fock::NormalFock) = fock|>rep

# We don't have to include other iterator methods like `iterate`, `length` and `lastindex` since they are 
# supported via the definition of `orderedmodes` and `dimension`.
# =========================================================================================================

# =======================================================================================
# NormalFock indexing
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

# ================================================================================================
# NormalFock common APIs for FockSpace

""" There are no concrete subspace in `NormalFock`, but it is a subspace itself. """
subspacecount(fock::NormalFock) = 1
# ================================================================================================

mode = quantize(point, 1)|>first
fock = NormalFock(Subset(mode))

fock[mode]

getindexk(crystal::Crystal, k::Momentum) = vec(k) .* size(crystal)

kindexoperator(crystal::Crystal) = [1, (size(crystal)[1:end-1]|>cumprod)...]

Base.getindex(crystal::Crystal, k::Momentum)::Integer = kindexoperator(crystal)' * getindexk(crystal, k) + 1

Base.:in(k::Momentum, crystal::Crystal) = getspace(k) == convert(MomentumSpace, getspace(crystal)) && getindexk(crystal, k)|>sum|>isinteger

k = bz[129]
k ∈ crystal
kindex(crystal, bz[123214])
