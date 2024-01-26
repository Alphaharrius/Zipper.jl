"""
    modeselectionbythreshold(threshold::Float64)

A selection strategy for `localisometries` that selects the modes with correlation eigenvalues γ
that is localized within the local region ∀ γ ≈ 0 or γ ≈ 1 up to a `threshold`. Extension of "frozenselectionbythreshold" in renormalization by including courier modes

### Input
- `threshold::Float64`: The threshold for the correlation eigenvalues.

### Output
A function representing the selection strategy with an input of the local correlations 𝐶ᵣ `FockMap` and outputs a `Dict{Symbol, FockMap}`
keyed by three grouping symbols with their associated local isometries selected from the unitary that diagonalizes 𝐶ᵣ. The grouping symbol
`:filled` represents the eigenmodes that corresponds to γ ≈ 0; `:empty` represents the eigenmodes that corresponds to γ ≈ 1; 
`:courier` represents the eigenmodes that corresponds to 0 ≈< γ <≈ 1
"""
function modeselectionbythreshold(threshold::Float64)::Function
    function modesfockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
        spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
        filledmodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second < threshold, spectrum |> geteigenvalues)))
        emptymodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second > 1.0 - threshold, spectrum |> geteigenvalues)))
        couriermodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> threshold <= p.second <= 1.0 - threshold, spectrum |> geteigenvalues)))
        return Dict(:filled => columns(spectrum |> geteigenvectors, FockSpace(filledmodes)), :empty => columns(spectrum |> geteigenvectors, FockSpace(emptymodes)),
        :courier => columns(spectrum |> geteigenvectors, FockSpace(couriermodes)))
    end
    return modesfockmaps
end
export modeselectionbythreshold

function modeselectionbycount(count::Integer)
    function modefockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
        𝑈ᵣ::FockMap = eigvecsh(𝐶ᵣ)
        modes::Vector{Mode} = [orderedmodes(𝑈ᵣ.inspace)...]
        return Dict(:filled => columns(𝑈ᵣ, FockSpace(Subset(modes[1:count]))), :empty => columns(𝑈ᵣ, FockSpace(Subset(modes[(end - count + 1):end]))),
        :courier => columns(𝑈ᵣ, FockSpace(Subset(modes[count+1:(end - count)]))))
    end
    return modefockmaps
end
export modeselectionbycount

localisometries(
    correlations::FockMap, regionfock::FockSpace;
    selectionstrategy::Function = modeselectionbythreshold(1e-3))::Dict{Symbol, FockMap} = (
    regioncorrelations(correlations, regionfock) |> selectionstrategy)
export localisometries

# function circularregionmodes(correlations::FockMap, origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
#     currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthospace
#     physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
#     return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
# end
# export circularregionmodes

# function groupmodebydist

# function localdistillerhamiltonian(;
#     correlations::FockMap, restrictspace::FockSpace,
#     localisometryselectionstrategy::Function, manualeigenenergies::Dict{Symbol, <:Number} = Dict(:filled => -1, :empty => 1))::FockMap

#     localisometries::Dict{Symbol} = localfrozenisometries(correlations, restrictspace, selectionstrategy=localisometryselectionstrategy)
#     crystalprojectors::Dict{Symbol, FockMap} = Dict(
#         name => crystalprojector(localisometry=localisometries[name], crystalfock=correlations.inspace)
#         for (name, isometry) in localisometries)
#     return reduce(+, manualeigenenergies[name] * crystalprojector for (name, crystalprojector) in crystalprojectors)
# end
# export globaldistillerhamiltonian


# function groupmodesbydist(;
#     region::Subset{Offset}
#     regionfock::FockSpace,
#     center::Point,
#     samedistancethreshold::Int = 8)::Dict{Integer, Mode}

#     visualspace = region |> getspace |> euclidean
#     distancewithmode = ((norm(lineartransform(visualspace, (mode |> getpos)-center) |> vec),mode) for mode in regionfock)

#     df = DataFrame()
#     df.distance = [round(dist; digits=10) for (dist,_) in save]
#     df.mode = [mode for (_,mode) in save]
#     grouped_df = groupby(df, :distance)
#     store = Dict()
#     for (ind,group) in enumerate(grouped_df)
#         store[ind] = []
#         for (distance,mode) in zip(group.distance,group.mode)
#             push!(store[ind],(distance,mode))
#         end
#     end

#     return store
# end


# function localwannierization(localbasis::FockMap, localseeds::FockMap, svdorthothreshold::Number = 1e-1)::FockMap
#     U, Σ, Vt = svd(localbasis'*localseeds)
#     minsvdvalue::Number = minimum(v for (_, v) in Σ)
#     println("min svdvalue", minsvdvalue)
#     precarioussvdvalues::Vector = []
#     if minsvdvalue < svdorthothreshold
#         push!(precarioussvdvalues, minsvdvalue)
#     end
#     if (precarioussvdvalues |> length) > 0
#         @warn "Precarious wannier projection with minimum svdvalue of $(precarioussvdvalues |> minimum)"
#     end
#     unitary::FockMap = U * Vt
#     wannierizedbasis = localbasis*unitary
#     wannierizedbasis = wannierizedbasis*_spatialmap(wannierizedbasis)'
#     return wannierizedbasis
# end