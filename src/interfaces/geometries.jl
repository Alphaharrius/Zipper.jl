""" Alias for a subset of `Offset` positions which is region like. """
Region = Subset{Offset}
export Region

struct Crystal <: AbstractSubset{Crystal}
    unitcell::Subset{Offset}
    sizes::Vector{Int64}
end
export Crystal

""" Interface to get the underlying region of an `Element`. """
getregion(::Element)::Region = notimplemented()
export getregion

""" Interface to get the underlying crystal of an `Element`. """
getcrystal(::Element) = notimplemented()
export getcrystal