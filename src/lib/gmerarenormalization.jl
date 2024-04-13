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
# function modeselectionbythreshold(threshold::Float64)::Function
#     function modesfockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, Dict{Symbol, Any}}
#         spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
#         evals = spectrum |> geteigenvalues
#         filledmodes::Subset{Mode} = Subset(pchosen.first for pchosen in filter(p -> p.second < threshold, spectrum |> geteigenvalues))
#         filledmodewifevaldict = Dict(mode=> evals[mode] for mode in filledmodes)
#         println("no of distillable filledmodes ", length(filledmodes))
#         emptymodes::Subset{Mode} = Subset(pchosen.first for pchosen in filter(p -> p.second > 1.0 - threshold, spectrum |> geteigenvalues))
#         emptymodewifevaldict =  Dict(mode=> evals[mode] for mode in emptymodes)
#         println("no of distillable emptymodes ", length(emptymodes))
#         couriermodes::Subset{Mode} = Subset(pchosen.first for pchosen in filter(p -> threshold <= p.second <= 1.0 - threshold, spectrum |> geteigenvalues))
#         couriermodewifevaldict =  Dict(mode=> evals[mode] for mode in couriermodes)
#         println("no of couriermodes ", length(couriermodes))
#         # return Dict(:filled => columns(spectrum |> geteigenvectors, FockSpace(filledmodes)), :empty => columns(spectrum |> geteigenvectors, FockSpace(emptymodes)),
#         # :courier => columns(spectrum |> geteigenvectors, FockSpace(couriermodes)))
#         return Dict(:filled => Dict(:isometry =>columns(spectrum |> geteigenvectors, FockSpace(filledmodes)),:modewifevaldict =>filledmodewifevaldict), 
#         :empty => Dict(:isometry =>columns(spectrum |> geteigenvectors, FockSpace(emptymodes)),:modewifevaldict =>emptymodewifevaldict),
#         :courier => Dict(:isometry =>columns(spectrum |> geteigenvectors, FockSpace(couriermodes)),:modewifevaldict =>couriermodewifevaldict))
#     end
#     return modesfockmaps
# end
# export modeselectionbythreshold

function modeselectionbythreshold(threshold::Float64)::Function
    @info ("modeselection by threshold")
    function modesfockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, Any}
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

# function modeselectionbycount(count::Integer)
#     function modefockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
#         𝑈ᵣ::FockMap = eigvecsh(𝐶ᵣ)
#         modes::Vector{Mode} = [orderedmodes(𝑈ᵣ.inspace)...]
#         return Dict(:filled => columns(𝑈ᵣ, FockSpace(Subset(modes[1:count]))), :empty => columns(𝑈ᵣ, FockSpace(Subset(modes[(end - count + 1):end]))),
#         :courier => columns(𝑈ᵣ, FockSpace(Subset(modes[count+1:(end - count)]))))
#     end
#     return modefockmaps
# end

# function modeselectionbycount(count::Integer)::Function
#     function modefockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, Dict{Symbol, Any}}
#         spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
#         evals = spectrum |> geteigenvalues
#         sortedmodeandevalpairs = sort!(collect(evals), by = x->x.second)
#         filledmodes::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[1:count])
#         filledmodewifevaldict = Dict(mode=> evals[mode] for mode in filledmodes)
#         reffilledeval = sortedmodeandevalpairs[count][2]
#         println("ref filled eigenvalue ", reffilledeval)
#         emptymodes::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[(end - count + 1):end])
#         emptymodewifevaldict =  Dict(mode=> evals[mode] for mode in emptymodes)
#         refemptyeval = sortedmodeandevalpairs[end - count + 1][2]
#         println("ref empty eigenvalue ", refemptyeval)
#         couriermodes::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[count+1:end-count])
#         couriermodewifevaldict =  Dict(mode=> evals[mode] for mode in couriermodes)
#         println("no of couriermodes ", length(couriermodes))
#         # return Dict(:filled => columns(spectrum |> geteigenvectors, FockSpace(filledmodes)), :empty => columns(spectrum |> geteigenvectors, FockSpace(emptymodes)),
#         # :courier => columns(spectrum |> geteigenvectors, FockSpace(couriermodes)))
#         return Dict(:filled => Dict(:isometry =>columns(spectrum |> geteigenvectors, FockSpace(filledmodes)),:modewifevaldict =>filledmodewifevaldict), 
#         :empty => Dict(:isometry =>columns(spectrum |> geteigenvectors, FockSpace(emptymodes)),:modewifevaldict =>emptymodewifevaldict),
#         :courier => Dict(:isometry =>columns(spectrum |> geteigenvectors, FockSpace(couriermodes)),:modewifevaldict =>couriermodewifevaldict))
#     end
#     return modefockmaps
# end
# export modeselectionbycount

function modeselectionbycount(count::Integer)::Function
    @info ("modeselection by count")
    function modefockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, Any}
        spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
        evals = spectrum |> geteigenvalues
        sortedmodeandevalpairs = sort!(collect(evals), by = x->x.second)
        filledmodes::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[1:count])
        reffilledeval = sortedmodeandevalpairs[count][2]
        println("ref filled eigenvalue ", reffilledeval)
        emptymodes::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[(end - count + 1):end])
        refemptyeval = sortedmodeandevalpairs[end - count + 1][2]
        println("ref empty eigenvalue ", refemptyeval)
        couriermodes::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[count+1:end-count])
        println("no of couriermodes ", length(couriermodes))
        return Dict(:filled => columns(spectrum |> geteigenvectors, FockSpace(filledmodes)), :empty => columns(spectrum |> geteigenvectors, FockSpace(emptymodes)),
        :courier => columns(spectrum |> geteigenvectors, FockSpace(couriermodes)))
    end
    return modefockmaps
end
export modeselectionbycount

# function modeselection1stbycountthenbythreshold(count::Integer,threshold::Float64)::Function
#     function modefockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, Dict{Symbol, Any}}
#         spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
#         evals = spectrum |> geteigenvalues
#         sortedmodeandevalpairs = sort!(collect(evals), by = x->x.second)
#         filledmodesbycount::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[1:count])
#         reffilledeval = sortedmodeandevalpairs[count][2]
#         println("ref filled eigenvalue ", reffilledeval)
#         emptymodesbycount::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[(end - count + 1):end])
#         refemptyeval = sortedmodeandevalpairs[end - count + 1][2]
#         println("ref empty eigenvalue ", refemptyeval)
#         filledmodesbythreshold::Subset{Mode} = Subset(pchosen.first for pchosen in filter(p -> p.second-reffilledeval < threshold, spectrum |> geteigenvalues))
#         emptymodesbythreshold::Subset{Mode} = Subset(pchosen.first for pchosen in filter(p -> 1-threshold < p.second+1-refemptyeval, spectrum |> geteigenvalues))
#         filledmodes = filledmodesbycount + filledmodesbythreshold
#         filledmodewifevaldict = Dict(mode=> evals[mode] for mode in filledmodes)
#         println("no of distillable filledmodes ", length(filledmodes))
#         emptymodes = emptymodesbycount + emptymodesbythreshold
#         emptymodewifevaldict =  Dict(mode=> evals[mode] for mode in emptymodes)
#         println("no of distillable emptymodes ", length(emptymodes))
#         couriermodes::Subset{Mode} = Subset(pchosen.first for pchosen in filter(p -> threshold + reffilledeval <= p.second <= refemptyeval - threshold, spectrum |> geteigenvalues))
#         couriermodewifevaldict =  Dict(mode=> evals[mode] for mode in couriermodes)
#         println("no of couriermodes ", length(couriermodes))
#         # return Dict(:filled => columns(spectrum |> geteigenvectors, FockSpace(filledmodes)), :empty => columns(spectrum |> geteigenvectors, FockSpace(emptymodes)),
#         # :courier => columns(spectrum |> geteigenvectors, FockSpace(couriermodes)))
#         return Dict(:filled => Dict(:isometry =>columns(spectrum |> geteigenvectors, FockSpace(filledmodes)),:modewifevaldict =>filledmodewifevaldict), 
#         :empty => Dict(:isometry =>columns(spectrum |> geteigenvectors, FockSpace(emptymodes)),:modewifevaldict =>emptymodewifevaldict),
#         :courier => Dict(:isometry =>columns(spectrum |> geteigenvectors, FockSpace(couriermodes)),:modewifevaldict =>couriermodewifevaldict))
#     end
#     return modefockmaps
# end
# export modeselection1stbycountthenbythreshold

function modeselection1stbycountthenbythreshold(count::Integer,threshold::Float64)::Function
    @info ("modeselection 1st by count then by threshold")
    function modefockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, Any}
        spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
        evals = spectrum |> geteigenvalues
        modesdict = Dict()
        sortedmodeandevalpairs = sort!(collect(evals), by = x->x.second)
        filledmodesbycount::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[1:count])
        reffilledeval = sortedmodeandevalpairs[count][2]
        println("ref filled eigenvalue ", reffilledeval)
        emptymodesbycount::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[(end - count + 1):end])
        refemptyeval = sortedmodeandevalpairs[end - count + 1][2]
        println("ref empty eigenvalue ", refemptyeval)
        filledmodesbythreshold::Subset{Mode} = Subset(pchosen.first for pchosen in filter(p -> p.second-reffilledeval < threshold, spectrum |> geteigenvalues))
        emptymodesbythreshold::Subset{Mode} = Subset(pchosen.first for pchosen in filter(p -> 1-threshold < p.second+1-refemptyeval, spectrum |> geteigenvalues))
        filledmodes = filledmodesbycount + filledmodesbythreshold
        println("no of distillable filledmodes ", length(filledmodes))
        modesdict[:filled] = columns(spectrum |> geteigenvectors, FockSpace(filledmodes))
        emptymodes = emptymodesbycount + emptymodesbythreshold
        println("no of distillable emptymodes ", length(emptymodes))
        modesdict[:empty] = columns(spectrum |> geteigenvectors, FockSpace(emptymodes))
        filteredcourier = filter(p -> threshold + reffilledeval <= p.second <= refemptyeval - threshold, spectrum |> geteigenvalues)
        if isempty(filteredcourier)
            println("no couriermodes remain!")
        else
            couriermodes::Subset{Mode} = Subset(pchosen.first for pchosen in filteredcourier)
            println("no of couriermodes ", length(couriermodes))
            modesdict[:courier] = columns(spectrum |> geteigenvectors, FockSpace(couriermodes))
        end
        return modesdict
    end
    return modefockmaps
end
export modeselection1stbycountthenbythreshold

function modeselectiononemodeplusepilon(threshold::Float64,epsilon::Float64)::Function
    @info ("modeselection 1 mode plus epsilon")
    function modefockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, Any}
        spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
        evals = spectrum |> geteigenvalues
        modesdict = Dict()
        sortedmodeandevalpairs = sort!(collect(evals), by = x->x.second)
        reffilledeval = sortedmodeandevalpairs[1][2]
        println("ref filled eigenvalue ", reffilledeval)
        refemptyeval = sortedmodeandevalpairs[end][2]
        println("ref empty eigenvalue ", refemptyeval)
        if reffilledeval<threshold
            firstfilledmode::Subset{Mode} = Subset(sortedmodeandevalpairs[1][1])
            epsilonfilledmodes::Subset{Mode} = Subset(pchosen.first for pchosen in filter(p -> p.second-reffilledeval < epsilon, spectrum |> geteigenvalues))
            filledchosenevals = [pchosen.second for pchosen in filter(p -> p.second-reffilledeval < epsilon, spectrum |> geteigenvalues)]
            println("chosen empty evals", filledchosenevals)
            filledmodes = firstfilledmode + epsilonfilledmodes
            println("no of distillable filledmodes ", length(filledmodes))
            modesdict[:filled] = columns(spectrum |> geteigenvectors, FockSpace(filledmodes))
            filledepsilon = epsilon
        else
            reffilledeval = 0
            filledepsilon = 0
            @info("the first filled mode is not within threshold => no filled mode")
        end
        if refemptyeval>1-threshold
            firstemptymode::Subset{Mode} = Subset(sortedmodeandevalpairs[end][1])
            epsilonemptymodes::Subset{Mode} = Subset(pchosen.first for pchosen in filter(p -> 1-epsilon < p.second+1-refemptyeval, spectrum |> geteigenvalues))
            emptychosenevals = [pchosen.second for pchosen in filter(p -> 1-epsilon < p.second+1-refemptyeval, spectrum |> geteigenvalues)]
            println("chosen empty evals", emptychosenevals)
            emptymodes = firstemptymode + epsilonemptymodes
            println("no of distillable emptymodes ", length(emptymodes))
            modesdict[:empty] = columns(spectrum |> geteigenvectors, FockSpace(emptymodes))
            emptyepsilon = epsilon
        else
            refemptyeval=1
            emptyepsilon = 0
            @info("the first empty mode is not within threshold => no empty mode")
        end
        
        
        filteredcourier = filter(p -> reffilledeval + filledepsilon <= p.second <= refemptyeval - emptyepsilon, spectrum |> geteigenvalues)
        if isempty(filteredcourier)
            println("no couriermodes remain!")
        else
            couriermodes::Subset{Mode} = Subset(pchosen.first for pchosen in filteredcourier)
            println("no of couriermodes ", length(couriermodes))
            modesdict[:courier] = columns(spectrum |> geteigenvectors, FockSpace(couriermodes))
        end
        return modesdict
    end
    return modefockmaps
end
export modeselectiononemodeplusepilon

localisometries(
    correlations::FockMap, regionfock::FockSpace;
    selectionstrategy::Function = modeselectionbythreshold(1e-3))::Dict{Symbol, Any} = (
    regioncorrelations(correlations, regionfock) |> selectionstrategy)
export localisometries


function fullftmap(correlations::FockMap)
    crystal = correlations|>getinspace|>getcrystal
    crystalpoints::Subset{Offset} = latticepoints(crystal)
    crystalfock::CrystalFock = correlations|>getinspace
    realspacemodes::RegionFock = RegionFock(spanoffset(correlations|>getinspace|>unitcellfock|>orderedmodes, crystalpoints))
    fouriermap::FockMap = fourier(crystalfock, realspacemodes)/ (crystalfock|>subspacecount|>sqrt)
    return fouriermap
end
export fullftmap

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
        return FockMap(localisometry, inspace=inspace, permute=false)
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
function groupmodesbydistwifb(;
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
    df.b = [mode |> getattr(:b) for (_,mode) in distancewithmode]
    grouped_df = groupby(df, :distance)
    store = Dict()
    for (ind,group) in enumerate(grouped_df)
        store[ind] = []
        for (distance,mode,b) in zip(group.distance,group.mode,group.b)
            push!(store[ind],(distance,mode,b))
        end
    end
    return store
end
export groupmodesbydistwifb

function modewifdiagpairs(;projector::FockMap,modes::Subset{Mode}=Subset(projector|>getinspace),noofchosen::Number=modes|>length)
    modewifdiag = Dict(mode=>norm(projector[mode,mode]) for mode in modes)
    return(Dict([pair[1]=>pair[2] for pair in sort!(collect(modewifdiag), by = x->x.second)][end-noofchosen+1:end]))
end
export modewifdiagpairs

function groupmodesbydiag(modewifdiagpairs::Dict{Mode,Float64},samedistancethreshold::Int = 5)
    df = DataFrame()
    df.diag = [round(diag; digits=samedistancethreshold) for (_,diag) in modewifdiagpairs]
    df.mode = [mode for (mode,_) in modewifdiagpairs]
    grouped_df = groupby(df, :diag)
    store = Dict()
    for (ind,group) in enumerate(grouped_df)
        store[ind] = []
        for (diag,mode) in zip(group.diag,group.mode)
            push!(store[ind],(diag,mode))
        end
    end
    return store
end
export groupmodesbydiag

# function genposX(regioncorrelations::FockMap)
#     modes = regioncorrelations|>getinspace
#     modewithpos = Dict(m=>((m|>getattr(:b))+(m|>getattr(:r))).pos[1] for m in modes)
#     sortedmodewithpos = [pair for pair in sort!(collect(modewithpos), by = x->x.second)]
#     diagentry = ComplexF64[]
#     xpos = ComplexF64(0)
#     ref = sortedmodewithpos[1][2]
#     for pair in sortedmodewithpos
#         if round(pair[2], digits=5) == round(ref, digits=5)
#             push!(diagentry,xpos)
#         else
#             xpos = xpos+1
#             push!(diagentry,xpos)
#             ref = pair[2]
#         end
#     end
#     rep = spdiagm(diagentry)
#     regionfock = RegionFock(Subset(pair[1] for pair in sortedmodewithpos))
#     return FockMap(regionfock,regionfock,rep)
# end
# export genposX

# function genposY(regioncorrelations::FockMap)
#     modes = regioncorrelations|>getinspace
#     modewithpos = Dict(m=>((m|>getattr(:b))+(m|>getattr(:r))).pos[2] for m in modes)
#     sortedmodewithpos = [pair for pair in sort!(collect(modewithpos), by = x->x.second)]
#     diagentry = ComplexF64[]
#     xpos = ComplexF64(0)
#     ref = sortedmodewithpos[1][2]
#     for pair in sortedmodewithpos
#         if round(pair[2], digits=5) == round(ref, digits=5)
#             push!(diagentry,xpos)
#         else
#             xpos = xpos+1
#             push!(diagentry,xpos)
#             ref = pair[2]
#         end
#     end
#     rep = spdiagm(diagentry)
#     regionfock = RegionFock(Subset(pair[1] for pair in sortedmodewithpos))
#     return FockMap(regionfock,regionfock,rep)
# end
# export genposY

# function localwannierseedslists(modebydistwifb,localiso::Dict{Symbol,FockMap})
#     modeorderedbydist = [modewifdistandb[2] for r in range(1,length(modebydistwifb)) for modewifdistandb in modebydistwifb[length(modebydistwifb)-r+1]]
#     borderedbydist = [modewifdistandb[3] for r in range(1,length(modebydistwifb)) for modewifdistandb in modebydistwifb[length(modebydistwifb)-r+1]]
#     nooffilledmodes = localiso[:filled]|>getinspace|>dimension
#     noofemptymodes = localiso[:empty]|>getinspace|>dimension
#     emptyseedslist = modeorderedbydist[1:noofemptymodes]
#     filledseedslist = modeorderedbydist[noofemptymodes+1:noofemptymodes+nooffilledmodes]
#     courierseedslist = modeorderedbydist[noofemptymodes+nooffilledmodes+1:length(modeorderedbydist)]
#     emptybslist = borderedbydist[1:noofemptymodes]
#     filledbslist = borderedbydist[noofemptymodes+1:noofemptymodes+nooffilledmodes]
#     courierbslist = borderedbydist[noofemptymodes+nooffilledmodes+1:length(modeorderedbydist)]
#     return Dict(:filled => filledseedslist, :empty => emptyseedslist, :courier => courierseedslist,
#     :bempty=>emptybslist,:bfilled=>filledbslist,:bcourier=>courierbslist)
# end

# function localwannierseedslists(localiso::Dict{Symbol,FockMap})
#     refmode = Subset(localiso[:empty]|>getoutspace)
#     filledseedslist = []
#     emptyseedslist = []
#     localisoempty = localiso[:empty]
#     localisofilled = localiso[:filled]
#     for inmode in localisofilled|>getinspace
#         chosenmode = sort!(collect(Dict((mode,norm(sum(localisofilled[mode,inmode]|>rep))) for mode in refmode)), by = x->x.second)[end].first
#         push!(filledseedslist,chosenmode)
#         refmode = refmode-Subset(chosenmode)
#     end
#     for inmode in localisoempty|>getinspace
#         chosenmode = sort!(collect(Dict((mode,norm(sum(localisoempty[mode,inmode]|>rep))) for mode in refmode)), by = x->x.second)[end].first
#         push!(emptyseedslist,chosenmode)
#         refmode = refmode-Subset(chosenmode)
#     end
#     courierseedslist = [mode for mode in refmode]
#     emptybslist = [mode |> getattr(:b) for mode in emptyseedslist]
#     filledbslist = [mode |> getattr(:b) for mode in filledseedslist]
#     courierbslist = [mode |> getattr(:b) for mode in courierseedslist]
#     return Dict(:filled => filledseedslist, :empty => emptyseedslist, :courier => courierseedslist,
#     :bempty=>emptybslist,:bfilled=>filledbslist,:bcourier=>courierbslist)
# end
# export localwannierseedslists

# function localwannierseedslists(localiso::Dict{Symbol, Dict{Symbol, Any}})
#     localisoempty = localiso[:empty][:isometry]
#     emptymodewitheval = localiso[:empty][:modewifevaldict]
#     noofemptyseeds = length(emptymodewitheval)
#     localisofilled = localiso[:filled][:isometry]
#     filledmodewitheval = localiso[:filled][:modewifevaldict]
#     nooffilledseeds = length(filledmodewitheval)
#     localisowannier = localiso[:courier][:isometry]
#     couriermodewitheval = localiso[:courier][:modewifevaldict]
#     noofcourierseeds = length(couriermodewitheval)
   
#     realspacemodes = Subset(localisoempty|>getoutspace)

#     # filledcountour = Dict(rmode=>sum([abs2(sum(localisofilled[rmode,inmode]|>rep)) for (inmode,eval) in filledmodewitheval]) for rmode in realspacemodes)
#     # filledseedslist = [pair.first for pair in sort!(collect(filledcountour), by = x->x.second)[end-nooffilledseeds+1:end]]
#     # realspacemodes = realspacemodes-Subset(filledseedslist)

#     courierecountour = Dict(rmode=>sum([abs2(sum(localisowannier[rmode,inmode]|>rep)) for (inmode,eval) in couriermodewitheval]) for rmode in realspacemodes)
#     courierseedslist = [pair.first for pair in sort!(collect(courierecountour), by = x->x.second)[end-noofcourierseeds+1:end]]
#     realspacemodes = realspacemodes-Subset(courierseedslist)

#     emptycountour = Dict(rmode=>sum([abs2(sum(localisoempty[rmode,inmode]|>rep)) for (inmode,eval) in emptymodewitheval]) for rmode in realspacemodes)
#     emptyseedslist = [pair.first for pair in sort!(collect(emptycountour), by = x->x.second)[end-noofemptyseeds+1:end]]
#     realspacemodes = realspacemodes-Subset(emptyseedslist)

#     filledseedslist = [rmode for rmode in realspacemodes]

#     emptybslist = [mode |> getattr(:b) for mode in emptyseedslist]
#     filledbslist = [mode |> getattr(:b) for mode in filledseedslist]
#     courierbslist = [mode |> getattr(:b) for mode in courierseedslist]
#     return Dict(:filled => filledseedslist, :empty => emptyseedslist, :courier => courierseedslist,
#     :bempty=>emptybslist,:bfilled=>filledbslist,:bcourier=>courierbslist)
# end
# export localwannierseedslists

function localwannierseedslists(localiso::Dict{Symbol, Any})
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

    # println(sortedfilledmodes[end]==sortedemptymodes[end])

    if isempty(sortedfilledmodes) & isempty(sortedemptymodes)
        @info("use up all the frozenmodes")
    end

    # filledmodes = frozenmodes[1:nooffilledseeds]
    # emptymodes = frozenmodes-filledmodes

    filledmodes = Subset(filledmodes)
    emptymodes = Subset(emptymodes)

    courierseedslist = [mode for mode in couriermodes]
    emptyseedslist = [mode for mode in emptymodes]
    filledseedslist = [mode for mode in filledmodes]

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

    # println(sortedfilledmodes[end]==sortedemptymodes[end])

    if isempty(sortedfilledmodes) & isempty(sortedemptymodes)
        @info("use up all the frozenmodes")
    end

    # filledmodes = frozenmodes[1:nooffilledseeds]
    # emptymodes = frozenmodes-filledmodes

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

# function entanglemententropy(eval::Number)
#     return -eval*log(eval)-(1-eval)*log(1-eval)
# end
# export entanglemententropy

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
    # wannierizedbasis = wannierizedbasis*gmeraspatialmap(wannierizedbasis)'
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

# function gmerastep1(correlations::CrystalFockMap,localcenterlist)
#     crystal::Crystal = getcrystal(correlations|>getinspace)
#     space::RealSpace = crystal|>getspace
#     function startingfulllocalwannierization(correlations::CrystalFockMap,localcenter::Offset,selectionstrategy::Function=modeselectionbycount(3))
#         shift = [1/2,1/2] ∈ space
#         localregion::Region =  Subset(pt+localcenter for pt in (crystal|>getunitcell))
#         localseedingfock::RegionFock = quantize(localregion,1)
#         modebydistwifb = groupmodesbydistwifb(region = localregion,regionfock = localseedingfock,center = localcenter+shift) 
#         localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstrategy)
#         localwannierlist = localwannierseedslists(localiso)
#         return startingcombinedlocalwannierization(localiso,localwannierlist,localseedingfock)
#     end

#     function fulllocalwannierization(correlations::CrystalFockMap,localcenter::Offset,selectionstrategy::Function,refselectionstrategydict)
#         localregion::Region =  Subset(pt+localcenter for pt in (crystal|>getunitcell))
#         localseedingfock::RegionFock = quantize(localregion,1)
#         modekeybyb = Dict(m|>getattr(:b)=>m for m in localseedingfock)
#         localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstrategy)
#         localwannierlist = localwannierseedslistsfromref(modekeybyb,refselectionstrategydict)
#         return combinedlocalwannierization(localiso,localwannierlist,localseedingfock)
#     end
#     wannierinfos =  Dict(localcenterlist[1]=>startingfulllocalwannierization(correlations,localcenterlist[1]))
#     refselectionstragedydict = Dict(:empty=>(wannierinfos[localcenterlist[1]][:bempty]),:filled=>(wannierinfos[localcenterlist[1]][:bfilled]),:courier=>(wannierinfos[localcenterlist[1]][:bcourier]))
#     for localcenter in localcenterlist[2:length(localcenterlist)]
#         wannierinfos[localcenter] = fulllocalwannierization(correlations,localcenter,modeselectionbycount(3),refselectionstragedydict)
#     end
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
#         :wannierinfos => wannierinfos,
#         :emptyisometry => wannieremptyisometry,
#         :filledisometry => wannierfilledisometry,
#         :courierisometry => wanniercourierisometry,
#         :couriercorrelations => couriercorrelations,
#         :correlations => purifiedcouriercorrelations,
#         :localemptyisometry => extendedwannierizedempty,
#         :localfilledisometry => extendedwannierizedfilled,
#         :localcourierisometry => extendedwannierizedcourier,
#         :refunictcellfockempty => refunictcellfockempty,
#         :refunictcellfockfilled => refunictcellfockfilled,
#         :refunictcellfockcourier => refunictcellfockcourier)
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

function gmerastep(rgblockedcorrelations::CrystalFockMap,correlations::CrystalFockMap,offsetlist,selectionstragedy)
    origcrystal::Crystal = getcrystal(rgblockedcorrelations|>getinspace)
    crystal::Crystal = getcrystal(correlations|>getinspace)
    space::RealSpace = crystal|>getspace

    function startinglocalisometries(correlations::CrystalFockMap,offset::Offset,selectionstragedy::Function=modeselectionbycount(6))
        @info ("starting step for one layer of gmera...")
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
        localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstragedy)
        return localiso
    end

    function startingfulllocalwannierization(correlations::CrystalFockMap,offset::Offset,selectionstragedy::Function=modeselectionbycount(6))
        @info ("starting step for one layer of gmera...")
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

        # return wannierinfos


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