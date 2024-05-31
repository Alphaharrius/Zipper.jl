function modeselectionbythreshold(threshold::Float64)::Function
    @info ("modeselection by threshold")
    function modesfockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
        spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
        modesdict = Dict()

        filteredfilled = filter(p -> p.second < threshold, spectrum |> geteigenvalues)
        if isempty(filteredfilled)
            println("no filledmodes within this threshold!")
        else
            filledevals = [pchosen.second for pchosen in filteredfilled]
            println("last filled eigenvalues ", filledevals[end])
            filledmodes::Subset{Mode} = Subset(pchosen.first for pchosen in filteredfilled)
            println("no of distillable filledmodes ", length(filledmodes))
            modesdict[:filled] = columns(spectrum |> geteigenvectors, FockSpace(filledmodes))
        end

        filteredempty = filter(p -> p.second > 1.0 - threshold, spectrum |> geteigenvalues)
        if isempty(filteredfilled)
            println("no emptymodes within this threshold!")
        else
            emptyevals = [pchosen.second for pchosen in filteredempty]
            println("last empty eigenvalues ", emptyevals[end])
            emptymodes::Subset{Mode} = Subset(pchosen.first for pchosen in filteredempty)
            println("no of distillable filledmodes ", length(emptymodes))
            modesdict[:empty] = columns(spectrum |> geteigenvectors, FockSpace(emptymodes))
        end

        filteredcourier = filter(p -> threshold <= p.second <= 1.0 - threshold, spectrum |> geteigenvalues)
        if isempty(filteredcourier)
            println("no couriermodes remain!")
        else
            couriermodes::Subset{Mode} = Subset(pchosen.first for pchosen in filteredcourier)
            println("no of couriermodes ", length(couriermodes))
            modesdict[:courier] = columns(spectrum |> geteigenvectors, FockSpace(couriermodes))
        end
        return modesdict
    end
    return modesfockmaps
end
export modeselectionbythreshold

function modeselectionbycount(count::Integer)::Function
    @info ("modeselection by count")
    function modefockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
        spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
        evals = spectrum |> geteigenvalues
        sortedmodeandevalpairs = sort!(collect(evals), by = x->x.second)
        filledmodes::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[1:count])
        reffilledeval = sortedmodeandevalpairs[count][2]
        println("ref filled eigenvalue ", reffilledeval)
        emptymodes::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[(end - count + 1):end])
        refemptyeval = sortedmodeandevalpairs[end - count + 1][2]
        println("ref empty eigenvalue ", refemptyeval)
        if length(emptymodes)+length(filledmodes) == length(sortedmodeandevalpairs)
            println("no couriermodes ")
            return Dict(:filled => columns(spectrum |> geteigenvectors, FockSpace(filledmodes)), :empty => columns(spectrum |> geteigenvectors, FockSpace(emptymodes)))
        else
            couriermodes::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[count+1:end-count])
            println("no of couriermodes ", length(couriermodes))
            return Dict(:filled => columns(spectrum |> geteigenvectors, FockSpace(filledmodes)), :empty => columns(spectrum |> geteigenvectors, FockSpace(emptymodes)),
            :courier => columns(spectrum |> geteigenvectors, FockSpace(couriermodes)))
        end
    end
    return modefockmaps
end
export modeselectionbycount

localisometries(
    correlations::FockMap, regionfock::RegionFock;
    selectionstrategy::Function = modeselectionbythreshold(1e-3))::Dict{Symbol, FockMap} = (
    regioncorrelations(correlations, regionfock) |> selectionstrategy)
export localisometries


function gmeracrystalisometries(; localisometry::FockMap, crystalfock::CrystalFock,
    addinspacemomentuminfo::Bool = false)

    crystal::Crystal = getcrystal(crystalfock)
    fouriermap::FockMap = fourier(crystalfock, localisometry|>getoutspace|>RegionFock) 
    # momentumfouriers::Base.Generator = rowsubmaps(fouriermap)
    bz::Subset{Momentum} = brillouinzone(crystal)

    function preprocesslocalisometry(k::Momentum)::FockMap
        if !addinspacemomentuminfo
            return localisometry
        end
        inspace::FockSpace = localisometry.inspace |> orderedmodes |> setattr(:k => k) |> removeattr(:r) |> FockSpace
        return FockMap(localisometry, inspace=inspace, performpermute=false)
    end

    # return Dict(k => kfourier * preprocesslocalisometry(k) for (k, kfourier) in zip(bz, momentumfouriers))
    # return Dict(k => fouriermap[getsubspace(crystalfock, k), :] * preprocesslocalisometry(k) for k in bz)
    isometries = paralleltasks(
        name="crystalisometries",
        tasks=(()->(k=>fouriermap[getsubspace(crystalfock, k), :]*preprocesslocalisometry(k)) for k in crystal|>brillouinzone),
        count=crystal|>vol)|>parallel

    return isometries
end
export gmeracrystalisometries

function modewifdiagpairs(;projector::FockMap,modes::Subset{Mode}=Subset(projector|>getinspace),noofchosen::Number=modes|>length)
    # modewifdiag = Dict(mode=>norm(projector[mode,mode]) for mode in modes)
    modewifdiag = Dict(mode=>norm(projector[mode,:]) for mode in modes)
    return(Dict([pair[1]=>pair[2] for pair in sort!(collect(modewifdiag), by = x->x.second)][end-noofchosen+1:end]))
end
export modewifdiagpairs


function localwannierseedslists(localiso::Dict{Symbol, FockMap})
    localisocourier = localiso[:courier]
    localisoempty = localiso[:empty]
    localisofilled = localiso[:filled]

    noofcourierseeds = length(localisocourier|>getinspace|>orderedmodes)
    noofemptyseeds = length(localisoempty|>getinspace|>orderedmodes)
    nooffilledseeds = length(localisofilled|>getinspace|>orderedmodes)

    projcourier = localisocourier*localisocourier'
    projfilled = localisofilled*localisofilled'
    projempty = localisoempty*localisoempty'

    chosencouriermodes = modewifdiagpairs(projector=projcourier,noofchosen=noofcourierseeds)
    realspacemodes = Subset(localisocourier|>getoutspace)
    couriermodes = Subset([pair[1] for pair in chosencouriermodes])

    frozenmodes = realspacemodes-couriermodes
    filledmodewifdiagpairs = modewifdiagpairs(projector=projfilled,modes=frozenmodes)
    emptymodewifdiagpairs = modewifdiagpairs(projector=projempty,modes=frozenmodes)

    sortedfilledmodes = [pair[1] for pair in sort!(collect(filledmodewifdiagpairs), by = x->x.second)]
    sortedemptymodes = [pair[1] for pair in sort!(collect(emptymodewifdiagpairs), by = x->x.second)]
    println([pair[2] for pair in sort!(collect(filledmodewifdiagpairs), by = x->x.second)])
    println([pair[2] for pair in sort!(collect(filledmodewifdiagpairs), by = x->x.second)])

    filledmodes = []
    emptymodes = []

    while noofemptyseeds>0 || nooffilledseeds>0
        if nooffilledseeds==0
            for mode in sortedemptymodes
                push!(emptymodes,mode)
                deleteat!(sortedemptymodes, findall(x->x==mode, sortedemptymodes))
                deleteat!(sortedfilledmodes, findall(x->x==mode, sortedfilledmodes))
                noofemptyseeds-=1
            end
        elseif noofemptyseeds==0
            for mode in sortedfilledmodes
                push!(filledmodes,mode)
                deleteat!(sortedemptymodes, findall(x->x==mode, sortedemptymodes))
                deleteat!(sortedfilledmodes, findall(x->x==mode, sortedfilledmodes))
                nooffilledseeds-=1
            end
        else
            if sortedfilledmodes[end] != sortedemptymodes[end]
                push!(filledmodes,sortedfilledmodes[end])
                push!(emptymodes,sortedemptymodes[end])
                deleteat!(sortedemptymodes, findall(x->x==sortedfilledmodes[end], sortedemptymodes))
                pop!(sortedfilledmodes)
                deleteat!(sortedfilledmodes, findall(x->x==sortedemptymodes[end], sortedfilledmodes))
                pop!(sortedemptymodes)
                nooffilledseeds-=1
                noofemptyseeds-=1
            else
                if round(filledmodewifdiagpairs[sortedfilledmodes[end]],digits=5)>round(emptymodewifdiagpairs[sortedemptymodes[end]],digits=5)
                    push!(filledmodes,sortedfilledmodes[end])
                    pop!(sortedfilledmodes)
                    pop!(sortedemptymodes)
                    nooffilledseeds-=1
                elseif round(filledmodewifdiagpairs[sortedfilledmodes[end]],digits=5)<round(emptymodewifdiagpairs[sortedemptymodes[end]],digits=5)
                    push!(emptymodes,sortedemptymodes[end])
                    pop!(sortedemptymodes)
                    pop!(sortedfilledmodes)
                    noofemptyseeds-=1
                else
                    if nooffilledseeds>noofemptyseeds
                        push!(filledmodes,sortedfilledmodes[end])
                        pop!(sortedfilledmodes)
                        pop!(sortedemptymodes)
                        nooffilledseeds-=1
                    elseif nooffilledseeds<noofemptyseeds
                        push!(emptymodes,sortedemptymodes[end])
                        pop!(sortedemptymodes)
                        pop!(sortedfilledmodes)
                        noofemptyseeds-=1
                    else
                        push!(filledmodes,sortedfilledmodes[end])
                        pop!(sortedfilledmodes)
                        pop!(sortedemptymodes)
                        nooffilledseeds-=1
                    end
                end
            end
        end
    end

    if isempty(sortedfilledmodes) & isempty(sortedemptymodes)
        @info("use up all the frozenmodes")
    end

    filledmodes = Subset(filledmodes)
    emptymodes = Subset(emptymodes)

    courierseedslist = [mode for mode in couriermodes]
    emptyseedslist = [mode for mode in emptymodes]
    filledseedslist = [mode for mode in filledmodes]

    emptyrslist = Subset([(mode |> getattr(:b)) + (mode |> getattr(:r)) for mode in emptyseedslist])
    filledrslist = Subset([(mode |> getattr(:b)) + (mode |> getattr(:r)) for mode in filledseedslist])
    courierrslist = Subset([(mode |> getattr(:b)) + (mode |> getattr(:r)) for mode in courierseedslist])
    display(visualize(emptyrslist))
    display(visualize(filledrslist))
    display(visualize(courierrslist))

    emptybslist = [mode |> getattr(:b) for mode in emptyseedslist]
    filledbslist = [mode |> getattr(:b) for mode in filledseedslist]
    courierbslist = [mode |> getattr(:b) for mode in courierseedslist]
    return Dict(:filled => filledseedslist, :empty => emptyseedslist, :courier => courierseedslist,
    :bempty=>emptybslist,:bfilled=>filledbslist,:bcourier=>courierbslist)
end
export localwannierseedslists

function localwannierseedslistswifoutcourier(localiso::Dict{Symbol, Any})
    localisoempty = localiso[:empty]
    localisofilled = localiso[:filled]

    noofemptyseeds = length(localisoempty|>getinspace|>orderedmodes)
    nooffilledseeds = length(localisofilled|>getinspace|>orderedmodes)

    projfilled = localisofilled*localisofilled'
    projempty = localisoempty*localisoempty'

    filledmodewifdiagpairs = modewifdiagpairs(projector=projfilled)
    emptymodewifdiagpairs = modewifdiagpairs(projector=projempty)

    sortedfilledmodes = [pair[1] for pair in sort!(collect(filledmodewifdiagpairs), by = x->x.second)]
    sortedemptymodes = [pair[1] for pair in sort!(collect(emptymodewifdiagpairs), by = x->x.second)]

    filledmodes = []
    emptymodes = []

    while noofemptyseeds>0 || nooffilledseeds>0
        if nooffilledseeds==0
            for mode in sortedemptymodes
                push!(emptymodes,mode)
                deleteat!(sortedemptymodes, findall(x->x==mode, sortedemptymodes))
                deleteat!(sortedfilledmodes, findall(x->x==mode, sortedfilledmodes))
                noofemptyseeds-=1
            end
        elseif noofemptyseeds==0
            for mode in sortedfilledmodes
                push!(filledmodes,mode)
                deleteat!(sortedemptymodes, findall(x->x==mode, sortedemptymodes))
                deleteat!(sortedfilledmodes, findall(x->x==mode, sortedfilledmodes))
                nooffilledseeds-=1
            end
        else
            if sortedfilledmodes[end] != sortedemptymodes[end]
                push!(filledmodes,sortedfilledmodes[end])
                push!(emptymodes,sortedemptymodes[end])
                deleteat!(sortedemptymodes, findall(x->x==sortedfilledmodes[end], sortedemptymodes))
                pop!(sortedfilledmodes)
                deleteat!(sortedfilledmodes, findall(x->x==sortedemptymodes[end], sortedfilledmodes))
                pop!(sortedemptymodes)
                nooffilledseeds-=1
                noofemptyseeds-=1
            else
                if round(filledmodewifdiagpairs[sortedfilledmodes[end]],digits=5)>round(emptymodewifdiagpairs[sortedemptymodes[end]],digits=5)
                    push!(filledmodes,sortedfilledmodes[end])
                    pop!(sortedfilledmodes)
                    pop!(sortedemptymodes)
                    nooffilledseeds-=1
                elseif round(filledmodewifdiagpairs[sortedfilledmodes[end]],digits=5)<round(emptymodewifdiagpairs[sortedemptymodes[end]],digits=5)
                    push!(emptymodes,sortedemptymodes[end])
                    pop!(sortedemptymodes)
                    pop!(sortedfilledmodes)
                    noofemptyseeds-=1
                else
                    if nooffilledseeds>noofemptyseeds
                        push!(filledmodes,sortedfilledmodes[end])
                        pop!(sortedfilledmodes)
                        pop!(sortedemptymodes)
                        nooffilledseeds-=1
                    elseif nooffilledseeds<noofemptyseeds
                        push!(emptymodes,sortedemptymodes[end])
                        pop!(sortedemptymodes)
                        pop!(sortedfilledmodes)
                        noofemptyseeds-=1
                    else
                        push!(filledmodes,sortedfilledmodes[end])
                        pop!(sortedfilledmodes)
                        pop!(sortedemptymodes)
                        nooffilledseeds-=1
                    end
                end
            end
        end
    end

    if isempty(sortedfilledmodes) & isempty(sortedemptymodes)
        @info("use up all the frozenmodes")
    end

    filledmodes = Subset(filledmodes)
    emptymodes = Subset(emptymodes)

    emptyseedslist = [mode for mode in emptymodes]
    filledseedslist = [mode for mode in filledmodes]

    emptybslist = [mode |> getattr(:b) for mode in emptyseedslist]
    filledbslist = [mode |> getattr(:b) for mode in filledseedslist]
    return Dict(:filled => filledseedslist, :empty => emptyseedslist,
    :bempty=>emptybslist,:bfilled=>filledbslist)
end
export localwannierseedslistswifoutcourier

function localwannierseedslistsfromref(modekeybyb,refselectionstragedydict)
    emptyseedslist = [modekeybyb[b] for b in refselectionstragedydict[:empty]]
    filledseedslist = [modekeybyb[b] for b in refselectionstragedydict[:filled]]
    courierseedslist = [modekeybyb[b] for b in refselectionstragedydict[:courier]]

    return Dict(:filled => filledseedslist, :empty => emptyseedslist, :courier => courierseedslist)
end
export localwannierseedslistsfromref

function localwannierseedslistsfromrefwifoutcourier(modekeybyb,refselectionstragedydict)
    filledseedslist = [modekeybyb[b] for b in refselectionstragedydict[:filled]]
    emptyseedslist = [modekeybyb[b] for b in refselectionstragedydict[:empty]]

    return Dict(:filled => filledseedslist, :empty => emptyseedslist)
end
export localwannierseedslistsfromrefwifoutcourier

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
    return wannierizedbasis
end
export localwannierization

function startingcombinedlocalwannierization(localiso,localwannierseedslists, localseedingfock)
    emptyseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in localwannierseedslists[:empty])] 
    filledseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in localwannierseedslists[:filled])] 
    courierseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in localwannierseedslists[:courier])]

    wannierempty = localwannierization(localiso[:empty], emptyseeds)
    wannierfilled = localwannierization(localiso[:filled], filledseeds)
    wanniercourier = localwannierization(localiso[:courier], courierseeds)

    return Dict(:wannierempty => wannierempty, :wannierfilled => wannierfilled, :wanniercourier => wanniercourier,
    :emptyseeds => emptyseeds, :filledseeds => filledseeds, :courierseeds => courierseeds,
    :bempty=>localwannierseedslists[:bempty],:bfilled=>localwannierseedslists[:bfilled],:bcourier=>localwannierseedslists[:bcourier])
end
export startingcombinedlocalwannierization

function startingcombinedlocalwannierizationwifoutcourier(localiso,localwannierseedslists, localseedingfock)
    emptyseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in localwannierseedslists[:empty])] 
    filledseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in localwannierseedslists[:filled])] 

    wannierempty = localwannierization(localiso[:empty], emptyseeds)
    wannierfilled = localwannierization(localiso[:filled], filledseeds)

    return Dict(:wannierempty => wannierempty, :wannierfilled => wannierfilled,
    :emptyseeds => emptyseeds, :filledseeds => filledseeds,
    :bempty=>localwannierseedslists[:bempty],:bfilled=>localwannierseedslists[:bfilled])
end
export startingcombinedlocalwannierizationwifoutcourier

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

function combinedlocalwannierizationwifoutcourier(localiso,localwannierseedslists, localseedingfock)
    emptyseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in localwannierseedslists[:empty])] 
    filledseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in localwannierseedslists[:filled])] 

    wannierempty = localwannierization(localiso[:empty], emptyseeds)
    wannierfilled = localwannierization(localiso[:filled], filledseeds)

    return Dict(:wannierempty => wannierempty, :wannierfilled => wannierfilled,
    :emptyseeds => emptyseeds, :filledseeds => filledseeds)
end
export combinedlocalwannierizationwifoutcourier

function gmerastep(rgblockedcorrelations::CrystalFockMap,correlations::CrystalFockMap,offsetlist,selectionstragedy)
    origcrystal::Crystal = getcrystal(rgblockedcorrelations|>getinspace)
    crystal::Crystal = getcrystal(correlations|>getinspace)
    space::RealSpace = crystal|>getspace
    hexagonalregion = gethexagonalregion(crystal=crystal)

    function startinglocalisometries(correlations::CrystalFockMap,offset::Offset,selectionstragedy::Function=modeselectionbycount(6))
        @info ("starting step for one layer of gmera...")
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
        localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstragedy)
        return localiso
    end

    function startingfulllocalwannierization(correlations::CrystalFockMap,offset::Offset,selectionstragedy::Function=modeselectionbycount(6))
        @info ("starting step for one layer of gmera...")
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
        # modebydistwifb = groupmodesbydistwifb(region = localregion,regionfock = localseedingfock,center = offset+shift) 
        localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstragedy)
        # localwannierlist = localwannierseedslists(modebydistwifb,localiso)
        if haskey(localiso,:courier)
            localwannierlist = localwannierseedslists(localiso)
            return startingcombinedlocalwannierization(localiso,localwannierlist,localseedingfock)
        else
            @info("starting gmera without courier")
            localwannierlist = localwannierseedslistswifoutcourier(localiso)
            return startingcombinedlocalwannierizationwifoutcourier(localiso,localwannierlist,localseedingfock)
        end
    end

    function fulllocalwannierization(correlations::CrystalFockMap,offset::Offset,selectionstrategy::Function,refselectionstrategydict)
        @info ("other steps for one layer of gmera...")
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
        modekeybyb = Dict(m|>getattr(:b)=>m for m in localseedingfock)
        localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstrategy)
        if haskey(localiso,:courier)
            localwannierlist = localwannierseedslistsfromref(modekeybyb,refselectionstrategydict)
            return combinedlocalwannierization(localiso,localwannierlist,localseedingfock)
        else
            @info("continuing gmera without courier")
            localwannierlist = localwannierseedslistsfromrefwifoutcourier(modekeybyb,refselectionstrategydict)
            return combinedlocalwannierizationwifoutcourier(localiso,localwannierlist,localseedingfock)
        end
    end

    # startinglocalisometries = startinglocalisometries(correlations,offsetlist[1],selectionstragedy)
    wannierinfos =  Dict(offsetlist[1]=>startingfulllocalwannierization(correlations,offsetlist[1],selectionstragedy))
    if haskey(wannierinfos[offsetlist[1]],:bcourier)
        refselectionstragedydict = Dict(:empty=>(wannierinfos[offsetlist[1]][:bempty]),:filled=>(wannierinfos[offsetlist[1]][:bfilled]),:courier=>(wannierinfos[offsetlist[1]][:bcourier]))
    else
        refselectionstragedydict = Dict(:empty=>(wannierinfos[offsetlist[1]][:bempty]),:filled=>(wannierinfos[offsetlist[1]][:bfilled]))
    end
    
    for offset in offsetlist[2:length(offsetlist)]
        wannierinfos[offset] = fulllocalwannierization(correlations,offset,selectionstragedy,refselectionstragedydict)
    end

    function globalwannierfunction(localisometry::FockMap)
        wanniercrystalisos = gmeracrystalisometries(localisometry=localisometry, crystalfock=correlations|>getoutspace,addinspacemomentuminfo=true)

        wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b)) for mode in localisometry|>getinspace |> orderedmodes)
        wanniercrystal::Crystal = Crystal(wannierunitcell, (correlations|>getoutspace|> getcrystal).sizes)

        blocks = Dict(((k, k)=>isometry) for (k, isometry) in wanniercrystalisos|>Dict)
        
        # blocks = paralleltasks(
        #     name="construct global isometry",
        #     tasks=(()->((k, k)=>isometry) for (k, isometry) in wanniercrystalisos|>Dict),
        #     count=wanniercrystalisos|>Dict|>length)|>parallel|>Dict
        
        globalwannierizedfunction::FockMap = crystalfockmap(correlations|>getoutspace|>getcrystal, wanniercrystal, blocks)
        # globalwannierizedfunction::FockMap = crystaldirectsum(outcrystal = correlations|>getoutspace|>getcrystal, incrystal=wanniercrystall,wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)

        return globalwannierizedfunction
    end

    if haskey(wannierinfos[offsetlist[1]],:bcourier)
        extendedwannierizedempty =  sum(wannierinfos[offset][:wannierempty] for offset in offsetlist)
        extendedwannierizedfilled =  sum(wannierinfos[offset][:wannierfilled] for offset in offsetlist)
        extendedwannierizedcourier =  sum(wannierinfos[offset][:wanniercourier] for offset in offsetlist)
        
        origin = [0, 0] ∈ space
        refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
        refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
        refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))


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
    else
        extendedwannierizedempty =  sum(wannierinfos[offset][:wannierempty] for offset in offsetlist)
        extendedwannierizedfilled =  sum(wannierinfos[offset][:wannierfilled] for offset in offsetlist)
        
        origin = [0, 0] ∈ space
        refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
        refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    
        wannieremptyisometry = globalwannierfunction(extendedwannierizedempty[:,refunictcellfockempty])
        wannierfilledisometry = globalwannierfunction(extendedwannierizedfilled[:,refunictcellfockfilled])


        return Dict(
            :emptyisometry => wannieremptyisometry,
            :filledisometry => wannierfilledisometry)
    end
end
export gmerastep


function globalwannierfunction(correlations::FockMap, localisometry::FockMap)
    wanniercrystalisos = gmeracrystalisometries(localisometry=localisometry, crystalfock=correlations|>getoutspace,addinspacemomentuminfo=true)

    wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b)) for mode in localisometry|>getinspace |> orderedmodes)
    wanniercrystal::Crystal = Crystal(wannierunitcell, (correlations|>getoutspace|> getcrystal).sizes)

    blocks = Dict(((k, k)=>isometry) for (k, isometry) in wanniercrystalisos|>Dict)
    
    globalwannierizedfunction::FockMap = crystalfockmap(correlations|>getoutspace|>getcrystal, wanniercrystal, blocks)
    
    return globalwannierizedfunction
end
export globalwannierfunction

