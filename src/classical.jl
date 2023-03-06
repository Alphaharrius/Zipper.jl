import Base.+

module classical

struct ClassicalSpace
    basis::Matrix{Number}
    dimension::Int8
end

"""
A general point in a [`ClassicalSpace`](@ref).
"""
struct Point
    """
    The position in spacial basis in the `parent_space`.
    """
    position::Vector{Number}
    """
    The parent space that includes this point.
    """
    parent_space::ClassicalSpace
end

"""
    point(position[, parent_space])

Create a new [`Point`](@ref) of `position` in the [`ClassicalSpace`](@ref) of `parent_space`.
"""
function point(position::Vector, parent_space::ClassicalSpace)::Point
    @assert(sizeof(position) == parent_space.dimension)
    return Point(position, parent_space)
end

function +(a::Point, b::Point)::Point
    @assert(typeof(a.parent_space) == typeof(b.parent_space))
    return Point(a.position + b.position, a.parent_space)
end

function __init__()
end

end
