""" Interface to get the dimension of an `Element`. """
dimension(::Element)::Integer = notimplemented()
export dimension

""" Interface to check if two `Element` objects have the same span. """
hassamespan(::Element, ::Element)::Bool = notimplemented()
export hassamespan

""" Interface to retrieve the space of an `Element`. """
getspace(::Element)::AffineSpace = notimplemented()
export getspace

""" Interface to retrieve the position of an `Element`. """
getpos(::Element) = notimplemented()
export getpos
