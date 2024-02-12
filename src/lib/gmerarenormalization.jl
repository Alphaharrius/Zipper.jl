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
        filledmodes::Subset{Mode} = Subset(pchosen.first for pchosen in filter(p -> p.second < threshold, spectrum |> geteigenvalues))
        emptymodes::Subset{Mode} = Subset(pchosen.first for pchosen in filter(p -> p.second > 1.0 - threshold, spectrum |> geteigenvalues))
        couriermodes::Subset{Mode} = Subset(pchosen.first for pchosen in filter(p -> threshold <= p.second <= 1.0 - threshold, spectrum |> geteigenvalues))
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

function gmeracrystalisometries(; localisometry::FockMap, crystalfock::CrystalFock,
    addinspacemomentuminfo::Bool = false)

    crystal::Crystal = getcrystal(crystalfock)
    transform::FockMap = fourier(crystalfock, localisometry|>getoutspace|>RegionFock)

    function preprocesslocalisometry(k::Momentum)::FockMap
        if !addinspacemomentuminfo
            return localisometry
        end
        inspace::FockSpace = localisometry.inspace |> orderedmodes |> setattr(:k => k) |> removeattr(:r) |> FockSpace
        return FockMap(localisometry, inspace=inspace, performpermute=false)
    end

    isometries = paralleltasks(
        name="crystalisometries",
        tasks=(()->(k=>transform[getsubspace(crystalfock, k), :]*preprocesslocalisometry(k)) for k in crystal|>brillouinzone),
        count=crystal|>vol)|>parallel

    return isometries
end
export crystalisometries

function gmeraspatialmap(fockmap::FockMap)::FockMap
    function spatialinmode(colmap::FockMap, ind::Integer)
        inmode::Mode = colmap |> getinspace |> first
        absmap::FockMap = colmap |> abs 
        modecenter::Offset = sort(absmap |> Zipper.columnspec, by=p->p.second |> real) |> last |> first |> getpos
        basis::Offset = modecenter |> basispoint
        offset::Offset = modecenter - basis
        return inmode |> setattr(:r => offset) |> setattr(:b => basis) 
    end

    spatialinspace::RegionFock = RegionFock(spatialinmode(fockmap[:, m],i) for (i,m) in fockmap |> getinspace |> enumerate)
    return idmap(spatialinspace, fockmap |> getinspace)
end

# Need further improvement
function groupmodesbydist(;
    region::Subset{Point{RealSpace}},
    regionfock::FockSpace,
    center::Point,
    samedistancethreshold::Int = 3)

    visualspace = region |> getspace |> euclidean
    physicalnorm = m -> lineartransform(visualspace, m |> getpos) |> norm
    distancewithmode = sort([(physicalnorm(mode-center),mode) for mode in regionfock], by = first, rev = true)
    df = DataFrame()
    df.distance = [round(dist; digits=samedistancethreshold) for (dist,_) in distancewithmode]
    df.mode = [mode for (_,mode) in distancewithmode]
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
export groupmodesbydist

function localwannierseedslists(modebydist,localiso::Dict{Symbol,FockMap})
    modeorderedbydist = [modewifdist[2] for r in range(1,length(modebydist)) for modewifdist in modebydist[length(modebydist)-r+1]]
    nooffilledmodes = localiso[:filled]|>getinspace|>dimension
    noofemptymodes = localiso[:empty]|>getinspace|>dimension
    emptyseedslist = modeorderedbydist[1:noofemptymodes]
    filledseedslist = modeorderedbydist[noofemptymodes+1:noofemptymodes+nooffilledmodes]
    courierseedslist = modeorderedbydist[noofemptymodes+nooffilledmodes+1:length(modeorderedbydist)]
    return Dict(:filled => filledseedslist, :empty => emptyseedslist, :courier => courierseedslist)
end
export localwannierseedslists


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
    # wannierizedbasis = wannierizedbasis*gmeraspatialmap(wannierizedbasis)'
    return wannierizedbasis
end
export localwannierization

function combinedlocalwannierization(localiso,localwannierseedslists, localseedingfock)
    emptyseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in localwannierseedslists[:empty])] 
    filledseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in localwannierseedslists[:filled])] 
    courierseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in localwannierseedslists[:courier])]

    wannierempty = localwannierization(localiso[:empty], emptyseeds)
    wannierfilled = localwannierization(localiso[:filled], filledseeds)
    wanniercourier = localwannierization(localiso[:courier], courierseeds)

    return Dict(:wannierempty => wannierempty, :wannierfilled => wannierfilled, :wanniercourier => wanniercourier,
    :emptyseeds => emptyseeds, :filledseeds => filledseeds, :courierseeds => courierseeds)
end
export combinedlocalwannierization

# function gmerastep1(correlations::CrystalFockMap,localcenterlist)
#     crystal::Crystal = getcrystal(correlations|>getinspace)
#     space::RealSpace = crystal|>getspace
#     function fulllocalwannierization(correlations::CrystalFockMap,localcenter::Offset,selectionstragedy::Function=modeselectionbycount(3))
#         shift = [1/2,1/2] ∈ space
#         localregion::Region =  Subset(pt+localcenter for pt in (crystal|>getunitcell))
#         localseedingfock::RegionFock = quantize(localregion,1)
#         modebydist = groupmodesbydist(region = localregion,regionfock = localseedingfock,center = localcenter+shift) 
#         localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstragedy)
#         localwannierlist = localwannierseedslists(modebydist,localiso)
#         return combinedlocalwannierization(localiso,localwannierlist,localseedingfock)
#     end
#     wannierinfos =  Dict(localcenter=>fulllocalwannierization(correlations,localcenter) for localcenter in localcenterlist)
    
#     extendedwannierizedempty =  sum(wannierinfos[center][:wannierempty] for center in localcenterlist)
#     extendedwannierizedfilled =  sum(wannierinfos[center][:wannierfilled] for center in localcenterlist)
#     extendedwannierizedcourier =  sum(wannierinfos[center][:wanniercourier] for center in localcenterlist)

#     # extendedemptyseeds = sum(wannierinfos[center][:emptyseeds] for center in localcenterlist)
#     # extendedfilledseeds = sum(wannierinfos[center][:filledseeds] for center in localcenterlist)
#     # extendedcourierseeds = sum(wannierinfos[center][:courierseeds] for center in localcenterlist)
    
#     origin = [0, 0] ∈ space
#     refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
#     refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
#     refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

#     # return wannierinfos

#     function globalwannierfunction(localisometry::FockMap)
#         wanniercrystalisos = gmeracrystalisometries(localisometry=localisometry, crystalfock=correlations|>getoutspace,addinspacemomentuminfo=true)

#         wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b)) for mode in localisometry|>getinspace |> orderedmodes)
#         wanniercrystall::Crystal = Crystal(wannierunitcell, (correlations|>getoutspace|> getcrystal).sizes)
#         globalwannierizedfunction::FockMap = crystaldirectsum(outcrystal = correlations|>getoutspace|>getcrystal, incrystal=wanniercrystall,wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)

#         return globalwannierizedfunction
#     end

#     wannieremptyisometry = globalwannierfunction(extendedwannierizedempty[:,refunictcellfockempty])
#     wannierfilledisometry = globalwannierfunction(extendedwannierizedfilled[:,refunictcellfockfilled])
#     wanniercourierisometry = globalwannierfunction(extendedwannierizedcourier[:,refunictcellfockcourier])

#     couriercorrelations = wanniercourierisometry' * correlations * wanniercourierisometry
#     couriercorrelationspectrum = couriercorrelations |> crystalspectrum
#     purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
#     purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

#     return Dict(
#         :emptyisometry => wannieremptyisometry,
#         :filledisometry => wannierfilledisometry,
#         :courierisometry => wanniercourierisometry,
#         :couriercorrelations => couriercorrelations,
#         :correlations => purifiedcouriercorrelations)
# end
# export gmerastep1

# function gmerastep2(rgblockedcorrelations::CrystalFockMap,correlations::CrystalFockMap,localcenterlist)
#     origcrystal::Crystal = getcrystal(rgblockedcorrelations|>getinspace)
#     crystal::Crystal = getcrystal(correlations|>getinspace)
#     space::RealSpace = crystal|>getspace
#     function fulllocalwannierization(correlations::CrystalFockMap,localcenter::Offset,selectionstragedy::Function=modeselectionbycount(3))
#         shift = [1/2,1/2] ∈ space
#         shiftedunitcell::Region =  Subset(pt+localcenter*2 for pt in (crystal|>getunitcell))
#         doubleunitcell::Region = shiftedunitcell+(crystal|>getunitcell)
#         shiftedorigunitcell::Region = Subset(pt+localcenter for pt in (origcrystal|>getunitcell))
#         localregion = intersect(doubleunitcell,shiftedorigunitcell)
#         localseedingfock::RegionFock = quantize(localregion,1)
#         modebydist = groupmodesbydist(region = localregion,regionfock = localseedingfock,center = localcenter+shift) 
#         localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstragedy)
#         localwannierlist = localwannierseedslists(modebydist,localiso)
#         return combinedlocalwannierization(localiso,localwannierlist,localseedingfock)
#     end
#     wannierinfos =  Dict(localcenter=>fulllocalwannierization(correlations,localcenter) for localcenter in localcenterlist)
    
#     extendedwannierizedempty =  sum(wannierinfos[center][:wannierempty] for center in localcenterlist)
#     extendedwannierizedfilled =  sum(wannierinfos[center][:wannierfilled] for center in localcenterlist)
#     extendedwannierizedcourier =  sum(wannierinfos[center][:wanniercourier] for center in localcenterlist)
    
#     origin = [0, 0] ∈ space
#     refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
#     refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
#     refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

#     # return wannierinfos

#     function globalwannierfunction(localisometry::FockMap)
#         wanniercrystalisos = gmeracrystalisometries(localisometry=localisometry, crystalfock=correlations|>getoutspace,addinspacemomentuminfo=true)

#         wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b)) for mode in localisometry|>getinspace |> orderedmodes)
#         wanniercrystall::Crystal = Crystal(wannierunitcell, (correlations|>getoutspace|> getcrystal).sizes)
#         globalwannierizedfunction::FockMap = crystaldirectsum(outcrystal = correlations|>getoutspace|>getcrystal, incrystal=wanniercrystall,wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)

#         return globalwannierizedfunction
#     end

#     wannieremptyisometry = globalwannierfunction(extendedwannierizedempty[:,refunictcellfockempty])
#     wannierfilledisometry = globalwannierfunction(extendedwannierizedfilled[:,refunictcellfockfilled])
#     wanniercourierisometry = globalwannierfunction(extendedwannierizedcourier[:,refunictcellfockcourier])

#     couriercorrelations = wanniercourierisometry' * correlations * wanniercourierisometry
#     couriercorrelationspectrum = couriercorrelations |> crystalspectrum
#     purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
#     purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

#     return Dict(
#         :emptyisometry => wannieremptyisometry,
#         :filledisometry => wannierfilledisometry,
#         :courierisometry => wanniercourierisometry,
#         :couriercorrelations => couriercorrelations,
#         :correlations => purifiedcouriercorrelations)
# end
# export gmerastep2

function gmerastep(rgblockedcorrelations::CrystalFockMap,correlations::CrystalFockMap,offsetlist)
    origcrystal::Crystal = getcrystal(rgblockedcorrelations|>getinspace)
    crystal::Crystal = getcrystal(correlations|>getinspace)
    space::RealSpace = crystal|>getspace
    function fulllocalwannierization(correlations::CrystalFockMap,offset::Offset,selectionstragedy::Function=modeselectionbythreshold(3))
        shift = [1/2,1/2] ∈ space
        shiftedorigunitcell::Region = Subset(pt+offset for pt in (origcrystal|>getunitcell))
        if offset == ([1/2,0] ∈ space) || offset == ([-1/2,0] ∈ space) || offset == ([0,1/2] ∈ space) || offset == ([0,-1/2] ∈ space)
            @info ("gmera step for vertical or horizontal offset...")
            shiftedunitcell::Region =  Subset(pt+offset*2 for pt in (crystal|>getunitcell))
            doubleunitcell::Region = shiftedunitcell+(crystal|>getunitcell)
            localregion = intersect(doubleunitcell,shiftedorigunitcell)
        elseif offset == ([1/2,1/2] ∈ space) || offset == ([-1/2,-1/2] ∈ space) || offset == ([1/2,-1/2] ∈ space) || offset == ([-1/2,1/2] ∈ space)
            @info ("gmera step for diagonal offset...")
            if offset == ([1/2,1/2] ∈ space)
                horizontal = [1,0] ∈ space
                vertical = [0,1] ∈ space
            elseif offset == ([-1/2,-1/2] ∈ space)
                horizontal = [-1,0] ∈ space
                vertical = [0,-1] ∈ space
            elseif offset == ([-1/2,1/2] ∈ space)
                horizontal = [-1,0] ∈ space
                vertical = [0,1] ∈ space
            elseif offset == ([1/2,-1/2] ∈ space)
                horizontal = [1,0] ∈ space
                vertical = [0,-1] ∈ space
            end
            shiftedhorizontalunitcell::Region =  Subset(pt+horizontal for pt in (crystal|>getunitcell))
            shiftedverticalunitcell::Region =  Subset(pt+vertical for pt in (crystal|>getunitcell))
            shifteddiagunitcell::Region =  Subset(pt+offset*2 for pt in (crystal|>getunitcell))
            allunitcell::Region = shiftedhorizontalunitcell+shiftedverticalunitcell+shifteddiagunitcell+(crystal|>getunitcell)
            localregion = intersect(allunitcell,shiftedorigunitcell)
        else
            @info ("gmera step without offset...")
            localregion = crystal|>getunitcell
        end
        localseedingfock::RegionFock = quantize(localregion,1)
        modebydist = groupmodesbydist(region = localregion,regionfock = localseedingfock,center = offset+shift) 
        localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstragedy)
        localwannierlist = localwannierseedslists(modebydist,localiso)
        return combinedlocalwannierization(localiso,localwannierlist,localseedingfock)
    end
    wannierinfos =  Dict(offset=>fulllocalwannierization(correlations,offset) for offset in offsetlist)
    
    extendedwannierizedempty =  sum(wannierinfos[offset][:wannierempty] for offset in offsetlist)
    extendedwannierizedfilled =  sum(wannierinfos[offset][:wannierfilled] for offset in offsetlist)
    extendedwannierizedcourier =  sum(wannierinfos[offset][:wanniercourier] for offset in offsetlist)
    
    origin = [0, 0] ∈ space
    refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

    # return wannierinfos

    function globalwannierfunction(localisometry::FockMap)
        wanniercrystalisos = gmeracrystalisometries(localisometry=localisometry, crystalfock=correlations|>getoutspace,addinspacemomentuminfo=true)

        wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b)) for mode in localisometry|>getinspace |> orderedmodes)
        wanniercrystall::Crystal = Crystal(wannierunitcell, (correlations|>getoutspace|> getcrystal).sizes)
        globalwannierizedfunction::FockMap = crystaldirectsum(outcrystal = correlations|>getoutspace|>getcrystal, incrystal=wanniercrystall,wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)

        return globalwannierizedfunction
    end

    wannieremptyisometry = globalwannierfunction(extendedwannierizedempty[:,refunictcellfockempty])
    wannierfilledisometry = globalwannierfunction(extendedwannierizedfilled[:,refunictcellfockfilled])
    wanniercourierisometry = globalwannierfunction(extendedwannierizedcourier[:,refunictcellfockcourier])

    couriercorrelations = wanniercourierisometry' * correlations * wanniercourierisometry
    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

    return Dict(
        :emptyisometry => wannieremptyisometry,
        :filledisometry => wannierfilledisometry,
        :courierisometry => wanniercourierisometry,
        :couriercorrelations => couriercorrelations,
        :correlations => purifiedcouriercorrelations)
end
export gmerastep