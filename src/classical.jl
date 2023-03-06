import Base.+

module classical

struct ClassicalSpace
end

"""
A general point in a [`ClassicalSpace`](@ref)
"""
struct Point
    position::Vector{Number}
    parent_space::ClassicalSpace
    dimension::Int8
end

function point(position::Vector, parent_space::ClassicalSpace)::Point
    return Point(position, parent_space, sizeof(position))
end

function +(a::Point, b::Point)::Point
    @assert(typeof(a.parent_space) == typeof(b.parent_space))
    return Point(a.position + b.position, a.parent_space, a.dimension)
end

function __init__()
end

end
