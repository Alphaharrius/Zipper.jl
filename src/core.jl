"""
Simply means an distinct type of object that can be represented by a concrete type `R`.
"""
abstract type Element{R <: Any} end
export Element

""" Restrict all print out of elements to just printing the type. """
Base.:show(io::IO, element::Element) = print(io, string(typeof(element)))

""" Reflexive relation. """
Base.:convert(::Type{T}, source::T) where {T <: Element} = source
""" The set representation of an `Element` is a set containing itself. """
Base.:convert(::Type{OrderedSet{T}}, source::T) where {T <: Element} = OrderedSet{T}([source])

Base.:*(a::A, b::B) where {A <: Element, B <: Element} = error("Composition algebra not defined for `$(a |> typeof)` with `$(b |> typeof)`!")
Base.:+(a::A, b::B) where {A <: Element, B <: Element} = error("Addition algebra not defined for `$(a |> typeof)` with `$(b |> typeof)`!")
Base.:-(::A, ::B) where {A <: Element, B <: Element} = error("Subtraction algebra not defined for `$(a |> typeof)` with `$(b |> typeof)`!")
Base.:/(::A, ::B) where {A <: Element, B <: Element} = error("Division algebra not defined for `$(a |> typeof)` with `$(b |> typeof)`!")

"""
    rep(source::Element{R}) where {R <: Any}

Retrieve the representation of the `source`, for each subtype of `Element{R}` are ensured to have a representation of type `R`.
This method will throw error if the subtype `T <: Element{R}` does not overload `Base.convert` to convert `T` to `R`.
"""
rep(source::Element{R}) where {R <: Any} = convert(R, source)
export rep
