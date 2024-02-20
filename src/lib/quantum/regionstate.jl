# ===========
# RegionState definition
struct RegionState{Dim} <: Element{FockMap}
    spstates::Dict{Mode, FockMap}
end
export RegionState
# ===========

# ===========
# RegionState essentials
Base.:convert(::Type{FockMap}, source::RegionState)::FockMap = reduce(+, spstate for (_, spstate) in spstates)
# ===========

# ===========
# RegionState iterator methods
Base.:iterate(state::RegionState, i...) = iterate(state.spstates, i...)
Base.:length(state::RegionState) = state.spstates |> length
# ===========

# ===========
# RegionState arithmetics
function Base.:+(a::RegionState, b::RegionState)
    dim::Integer = a|>getoutspace|>getregion|>getspace|>dimension
    combinedstates::Vector = [a.spstates..., b.spstates...]
    combinedinmodes::Base.Generator = (mode for (mode, _) in combinedstates)
    # This step ensures that any duplicated modes from a and b will be mapped to different :flavor.
    mergedmodes::Subset{Mode} = combinedinmodes|>mapmodes(m -> m)
    mappedstates = (mode=>FockMap(state, inspace=mode|>FockSpace, performpermute=false) for (mode, (_, state)) in zip(mergedmodes, combinedstates))
    return RegionState{dim}(mappedstates|>Dict)
end
# ===========

# ===========
# RegionState APIs
function regionalrestriction(crystalstate::FockMap, regionfock::RegionFock)::RegionState
    eigenmodes::Subset{Mode} = crystalstate |> getinspace |> unitcellfock |> orderedmodes

    function extractregionstate(mode::Mode)
        rightfourier::FockMap = fourier(crystalstate |> getinspace, mode |> RegionFock)
        leftfourier::FockMap = fourier(crystalstate |> getoutspace, regionfock)
        return leftfourier' * crystalstate * rightfourier
    end

    return Dict(mode => mode |> extractregionstate for mode in eigenmodes) |> RegionState{crystalstate |> getoutspace |> getcrystal |> dimension}
end
export regionalrestriction
# ===========

# ===========
# RegionState conversions
"""" Convert a state like `FockMap` with outspace type of `RegionFock` to a `RegionState`. """
function RegionState(localstates::SparseFockMap{RegionFock, <:FockSpace})
    decorated::FockMap = localstates * spatialmap(localstates)
    dim::Integer = decorated|>getoutspace|>getregion|>getspace|>dimension
    return Dict(mode=>decorated[:, mode] for mode in decorated|>getinspace)|>RegionState{dim}
end

function FockMap(regionstate::RegionState)
    statemap::FockMap = reduce(+, state for (_, state) in regionstate.spstates)
    return FockMap(statemap, outspace=statemap|>getoutspace|>RegionFock, inspace=statemap|>getinspace|>RegionFock)
end
# ===========

# ===========
# FockMap interfaces extensions
getinspace(state::RegionState) = FockSpace(m for (m, _) in state.spstates)
getoutspace(state::RegionState) = state.spstates |> first |> last |> getoutspace
# ===========

# ===========
# RegionState display
Base.:show(io::IO, state::RegionState) = print(io, string("$(typeof(state))(count=$(state|>getinspace|>dimension))"))
# ===========
