function groupmodesbydist(;
    region::Subset{Offset}
    regionfock::FockSpace,
    center::Point,
    samedistancethreshold::Int = 8)::Dict{Integer, Mode}

    visualspace = region |> getspace |> euclidean
    distancewithmode = ((norm(lineartransform(visualspace, (mode |> getpos)-center) |> vec),mode) for mode in regionfock)

    df = DataFrame()
    df.distance = [round(dist; digits=10) for (dist,_) in save]
    df.mode = [mode for (_,mode) in save]
    grouped_df = groupby(df, :distance)
    store = Dict()
    for (ind,group) in enumerate(grouped_df)
        store[ind] = []
        for (distance,mode) in zip(group.distance,group.mode)
            push!(store[ind],(distance,mode))
        end
    end

    return store
end


function localwannierization(localbasis::FockMap, localseeds::FockMap, svdorthothreshold::Number = 1e-1)::FockMap
    U, Σ, Vt = svd(localbasis'*localseeds)
    minsvdvalue::Number = minimum(v for (_, v) in Σ)
    println("min svdvalue", minsvdvalue)
    precarioussvdvalues::Vector = []
    if minsvdvalue < svdorthothreshold
        push!(precarioussvdvalues, minsvdvalue)
    end
    if (precarioussvdvalues |> length) > 0
        @warn "Precarious wannier projection with minimum svdvalue of $(precarioussvdvalues |> minimum)"
    end
    unitary::FockMap = U * Vt
    wannierizedbasis = localbasis*unitary
    wannierizedbasis = wannierizedbasis*_spatialmap(wannierizedbasis)'
    return wannierizedbasis
end