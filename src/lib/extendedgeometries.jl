function getsphericalregion(;from::Offset, generators::Subset{Offset}, symmetry::AffineTransform, radius::Real, metricspace::RealSpace)
    from |> latticeoff == from || error("`from` must be a lattice offset!")
    
    getspancount(generator::Offset)::Integer = radius / (generator |> euclidean |> norm) |> ceil |> Integer

    generated::Region = Subset(from |> getspace |> getorigin)
    println(generated |> collect)
    for vec in generators
        expanded::Region = Subset(p + vec * n for p in generated for n in 0:getspancount(vec))
        generated = filter(p -> lineartransform(metricspace, p) |> norm <= radius, expanded)
    end

    for transform in symmetry |> pointgroupelements
        generated += transform * generated
    end

    return generated + from
end
export getsphericalregion
