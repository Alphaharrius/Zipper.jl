# function modeselectionbythreshold(threshold::Float64)::Function
#     @info ("modeselection by threshold")
#     function modesfockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
#         spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
#         modesdict = Dict()

#         filteredfilled = filter(p -> p.second < threshold, spectrum |> geteigenvalues)
#         if isempty(filteredfilled)
#             println("no filledmodes within this threshold!")
#         else
#             filledevals = [pchosen.second for pchosen in filteredfilled]
#             println("last filled eigenvalues ", filledevals[end])
#             filledmodes::Subset{Mode} = Subset(pchosen.first for pchosen in filteredfilled)
#             println("no of distillable filledmodes ", length(filledmodes))
#             modesdict[:filled] = columns(spectrum |> geteigenvectors, FockSpace(filledmodes))
#         end

#         filteredempty = filter(p -> p.second > 1.0 - threshold, spectrum |> geteigenvalues)
#         if isempty(filteredfilled)
#             println("no emptymodes within this threshold!")
#         else
#             emptyevals = [pchosen.second for pchosen in filteredempty]
#             println("last empty eigenvalues ", emptyevals[end])
#             emptymodes::Subset{Mode} = Subset(pchosen.first for pchosen in filteredempty)
#             println("no of distillable filledmodes ", length(emptymodes))
#             modesdict[:empty] = columns(spectrum |> geteigenvectors, FockSpace(emptymodes))
#         end

#         filteredcourier = filter(p -> threshold <= p.second <= 1.0 - threshold, spectrum |> geteigenvalues)
#         if isempty(filteredcourier)
#             println("no couriermodes remain!")
#         else
#             couriermodes::Subset{Mode} = Subset(pchosen.first for pchosen in filteredcourier)
#             println("no of couriermodes ", length(couriermodes))
#             modesdict[:courier] = columns(spectrum |> geteigenvectors, FockSpace(couriermodes))
#         end
#         return modesdict
#     end
#     return modesfockmaps
# end
# export modeselectionbythreshold

# function modeselectionbycount(count::Integer)::Function
#     @info ("modeselection by count")
#     function modefockmaps(𝐶ᵣ::FockMap)::Dict{Symbol, FockMap}
#         spectrum::EigenSpectrum = 𝐶ᵣ |> eigspech
#         evals = spectrum |> geteigenvalues
#         sortedmodeandevalpairs = sort!(collect(evals), by = x->x.second)
#         filledmodes::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[1:count])
#         reffilledeval = sortedmodeandevalpairs[count][2]
#         println("ref filled eigenvalue ", reffilledeval)
#         emptymodes::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[(end - count + 1):end])
#         refemptyeval = sortedmodeandevalpairs[end - count + 1][2]
#         println("ref empty eigenvalue ", refemptyeval)
#         if length(emptymodes)+length(filledmodes) == length(sortedmodeandevalpairs)
#             println("no couriermodes ")
#             return Dict(:filled => columns(spectrum |> geteigenvectors, FockSpace(filledmodes)), :empty => columns(spectrum |> geteigenvectors, FockSpace(emptymodes)))
#         else
#             couriermodes::Subset{Mode} = Subset(pair[1] for pair in sortedmodeandevalpairs[count+1:end-count])
#             println("no of couriermodes ", length(couriermodes))
#             return Dict(:filled => columns(spectrum |> geteigenvectors, FockSpace(filledmodes)), :empty => columns(spectrum |> geteigenvectors, FockSpace(emptymodes)),
#             :courier => columns(spectrum |> geteigenvectors, FockSpace(couriermodes)))
#         end
#     end
#     return modefockmaps
# end
# export modeselectionbycount

function generatesystem(t_a,t_b,tₙ,tₕ,systemsize::Number)
    triangular = RealSpace([1 0; 1/2 sqrt(3)/2]')
    kspace = convert(MomentumSpace, triangular)

    pa = [1/3, 0] ∈ triangular
    pb = [0, 1/3] ∈ triangular
    pc = [0, 2/3] ∈ triangular
    pd = [1/3, 2/3] ∈ triangular
    pe = [2/3, 1/3] ∈ triangular
    pf = [2/3, 0] ∈ triangular

    pg = (pa + pb + pc + pd + pe + pf) / 6
    spatialsnappingcalibration((pa, pb, pc, pd, pe, pf, pg))

    unitcell = Subset(pa, pb, pc, pd, pe, pf)
    crystal = Crystal(unitcell, [systemsize, systemsize])
    modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
    m0, m1, m2, m3, m4, m5 = members(modes)

    onsite = [
    (m0, m0) => t_a,
    (m1, m1) => t_b,
    (m2, m2) => t_a,
    (m3, m3) => t_b,
    (m4, m4) => t_a,
    (m5, m5) => t_b
]

    nearestneighbor = [
        (m1, m0) => tₙ,
        (m1, m4|> setattr(:r => [-1, 0] ∈ triangular)) => tₙ,
        (m1, m2) => tₙ,
        (m3, m2) => tₙ,
        (m3, m0|> setattr(:r => [0, 1] ∈ triangular)) => tₙ,
        (m3, m4) => tₙ,
        (m5, m4) => tₙ,
        (m5, m2|> setattr(:r => [1, -1] ∈ triangular)) => tₙ,
        (m5, m0) => tₙ]

    haldaneterm = [
        (m1, m3|> setattr(:r => [-1, 0] ∈ triangular)) => tₕ,
        (m1, m3) => tₕ,
        (m1, m3|> setattr(:r => [0, -1] ∈ triangular)) => tₕ,
        (m3, m5|> setattr(:r => [-1, 1] ∈ triangular)) => tₕ,
        (m3, m5|> setattr(:r => [0, 1] ∈ triangular)) => tₕ,
        (m3, m5) => tₕ,
        (m5, m1) => tₕ,
        (m5, m1|> setattr(:r => [1, 0] ∈ triangular)) => tₕ,
        (m5, m1|> setattr(:r => [1, -1] ∈ triangular)) => tₕ,
        (m0, m4|> setattr(:r => [-1, 0] ∈ triangular)) => -tₕ,
        (m0, m4) => -tₕ,
        (m0, m4|> setattr(:r => [0, -1] ∈ triangular)) => -tₕ,
        (m2, m0|> setattr(:r => [-1, 1] ∈ triangular)) => -tₕ,
        (m2, m0|> setattr(:r => [0, 1] ∈ triangular)) => -tₕ,
        (m2, m0) => -tₕ,
        (m4, m2) => -tₕ,
        (m4, m2|> setattr(:r => [1, 0] ∈ triangular)) => -tₕ,
        (m4, m2|> setattr(:r => [1, -1] ∈ triangular)) => -tₕ]

    bonds::FockMap = bondmap([onsite..., nearestneighbor...,haldaneterm...])
    energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)

    groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=3)
    groundstateprojector = groundstates|>crystalprojector
    correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector
    H = CrystalFockMap(energyspectrum)

    return correlations,H
end
export generatesystem

localisometries(
    correlations::FockMap, regionfock::RegionFock;
    selectionstrategy::Function = modeselectionbythreshold(1e-3))::Dict{Symbol, FockMap} = (
    regioncorrelations(correlations, regionfock) |> selectionstrategy)
export localisometries

function separatefilledandemptymodes(spectrum::EigenSpectrum)
    evaldict = spectrum |> geteigenvalues
    filleddict = Dict(key=>value for (key,value) in evaldict if value<0.5)
    emptydict = Dict(key=>value for (key,value) in evaldict if value>0.5)
    return filleddict,emptydict
end
export separatefilledandemptymodes

function calculateavgoverlapofmodes(truncatedeigenvec,modes::Subset{Mode})
    overlaps = [(norm(columns(truncatedeigenvec, FockSpace(mode))))^2 for mode in modes]
    avgoverlap = sum(overlaps)/length(overlaps)
    return tuple(avgoverlap,modes)
end
export calculateavgoverlapofmodes

function sortgroupdictwifdist(dict::Dict{Mode, Float64},rev::Bool)
    sortedtupledata = sort([(value,key) for (key,value) in dict],rev=rev,by=first)
    refvalue = round(sortedtupledata[1][1],digits=8)
    result = []
    subresult = Subset(sortedtupledata[1][2])
    for pair in sortedtupledata
        if round(pair[1],digits=8)==refvalue
            subresult = subresult+Subset(pair[2])
        else
            append!(result,tuple([refvalue,subresult]))
            refvalue = round(pair[1],digits=8)
            subresult = Subset(pair[2])
        end
    end
    append!(result,tuple([refvalue,subresult]))

    return result
end
export sortgroupdictwifdist

function generatehexagonalregionref(systemsize::Number,scale::Scale)
    triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')

    pa = [1/3, 2/3] ∈ triangular
    pb = [2/3, 1/3] ∈ triangular
    pc = (pa + pb) / 2
    spatialsnappingcalibration((pa, pb, pc))

    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    unitcell = Subset(pa, pb)
    crystal = Crystal(unitcell, [systemsize, systemsize])
    scaledcrystal = scale*crystal
    scaledunitcell = scaledcrystal|>getunitcell
    hexagonalregionatorig = scaledunitcell+c3*scaledunitcell+c3^2*scaledunitcell
    return hexagonalregionatorig
end
export generatehexagonalregionref

function generaterefAandBsublatticeunitcellmodes(crystal::Crystal)
    Asublatticeunitcellmodes = Subset([m|>removeattr(:r) for (ind,m) in enumerate(members(quantize(crystal|>getunitcell,1)|>orderedmodes)) if isodd(ind)])
    Bsublatticeunitcellmodes = Subset([m|>removeattr(:r) for (ind,m) in enumerate(members(quantize(crystal|>getunitcell,1)|>orderedmodes)) if iseven(ind)])
    return Asublatticeunitcellmodes,Bsublatticeunitcellmodes
end
export generaterefAandBsublatticeunitcellmodes

function distinguishABsublatticemodesforhexagon(modes::Subset{Mode},Arefmodes::Subset{Mode},Brefmodes::Subset{Mode})
    modesdict = Dict(m|>removeattr(:r)=> m for m in modes )
    Amodes = sum([Subset(modesdict[am]) for am in Arefmodes if haskey(modesdict,am)])
    Bmodes = sum([Subset(modesdict[bm]) for bm in Brefmodes if haskey(modesdict,bm)]) 
    return Amodes,Bmodes
end
export distinguishABsublatticemodesforhexagon

function checkingdiffwifchosen(truncatedeigenvec, chosenmodes::Subset{Mode},newmodes::Subset{Mode})
    chosencolumns = columns(truncatedeigenvec, FockSpace(chosenmodes))
    U, Σ, Vt = svd(chosencolumns)
    # minsvdvalue::Number = minimum(v for (_, v) in Σ)
    # println("min svdvalue during checking before", minsvdvalue)
    projector = U*U'
    candidate = columns(truncatedeigenvec, FockSpace(newmodes))
    W, Σ, Xt = svd(candidate)
    diff = projector*W-W
    normlist = [norm(columns(diff,FockSpace(m))) for m in diff|>getinspace]
    return sum(normlist)/length(normlist)
end
export checkingdiffwifchosen

function choosemodesforwannierization(rankedandgroupedmodes,truncatedeigenvec,threshold::Float64,noofmodes::Number)
    chosenmodes = rankedandgroupedmodes[1]
    for modes in rankedandgroupedmodes[2:end]
        if length(chosenmodes)<noofmodes
            if length(chosenmodes)+length(modes)>noofmodes
                @info("more than need")
            else
                diff = checkingdiffwifchosen(truncatedeigenvec, chosenmodes,modes)
                if diff>threshold
                    chosenmodes = chosenmodes+modes
                    chosencolumns = columns(truncatedeigenvec, FockSpace(chosenmodes))
                    # U, Σ, Vt = svd(chosencolumns)
                    # minsvdvalue::Number = minimum(v for (_, v) in Σ)
                    # println("min svdvalue during checking after", minsvdvalue)
                else
                    @info("reject this candidate due to small diff", diff)
                end
            end
        else
            return chosenmodes
        end
    end
    @error("threshold too large")
end
export choosemodesforwannierization

# function findmaxofoverlap(truncatedeigenvec, chosenmodes,modes)
#     return maximum(abs.((columns(truncatedeigenvec,FockSpace(chosenmodes))'*columns(truncatedeigenvec,FockSpace(modes))).rep)) 
# end 

# function choosemodesforwannierization(rankedandgroupedmodes,truncatedeigenvec,threshold::Float64,noofmodes::Number)
#     chosenmodes = rankedandgroupedmodes[1]
#     for modes in rankedandgroupedmodes[2:end]
#         if length(chosenmodes)<noofmodes
#             if length(chosenmodes)+length(modes)>noofmodes
#                 @info("more than need")
#             else
#                 max = findmaxofoverlap(truncatedeigenvec, chosenmodes,modes)
#                 if max<threshold
#                     chosenmodes = chosenmodes+modes
#                     chosencolumns = columns(truncatedeigenvec, FockSpace(chosenmodes))
#                     U, Σ, Vt = svd(chosencolumns)
#                     minsvdvalue::Number = minimum(v for (_, v) in Σ)
#                     @info("min svdvalue during checking after", minsvdvalue)
#                 else
#                     @info("reject this candidate due to small max", max)
#                 end
#             end
#         else
#             return chosenmodes
#         end
#     end
#     @error("threshold too large")
# end
# export choosemodesforwannierization

function findfrozenandcourierrsmode(sortedrsmode,noofmodes)
    allrsmodes = sum([pair[2] for pair in sortedrsmode])
    frozenrsmodes = sortedrsmode[1][2]
    for pair in sortedrsmode[2:end]
        if length(frozenrsmodes)<noofmodes
            frozenrsmodes = frozenrsmodes+pair[2]
        elseif length(frozenrsmodes)>noofmodes
            @error("given noofmodes cannot be achieved")
        else
            courierrsmodes = allrsmodes-frozenrsmodes
            return frozenrsmodes,courierrsmodes
        end
    end
end
export findfrozenandcourierrsmode

function sortgroupdictwifvalue(dict::Dict{Mode, Float64},rev::Bool)
    sortedtupledata = sort([(value,key) for (key,value) in dict],rev=rev,by=first)
    refvalue = sortedtupledata[1][2]|>getattr(:eigenindex)
    result = []
    evals = []
    subresult = Subset(sortedtupledata[1][2])
    for pair in sortedtupledata
        if pair[2]|>getattr(:eigenindex)==refvalue
            subresult = subresult + Subset(pair[2])
            append!(evals,pair[1])
        else
            append!(result,[subresult])
            refvalue = pair[2]|>getattr(:eigenindex)
            evals = []
            subresult = Subset(pair[2])
        end
    end
    append!(result,[subresult])

    return result
end
export sortgroupdictwifvalue

function sortgroupdictwifvaluefilled(dict::Dict{Mode, Float64}, threshold::Float64)
    tupledata = sort([(value,key) for (key,value) in dict],by=first)
    refvalue = tupledata[1][2]|>getattr(:eigenindex)
    result = []
    evals = []
    subresult = Subset(tupledata[1][2])
    for pair in tupledata
        if pair[1]<threshold
            if pair[2]|>getattr(:eigenindex)==refvalue
                subresult = subresult + Subset(pair[2])
                append!(evals,pair[1])
            else
                append!(result,[subresult])
                refvalue = pair[2]|>getattr(:eigenindex)
                evals = []
                subresult = Subset(pair[2])
            end
        else 
            return result
        end
    end
    append!(result,[subresult])
    return result
end
export sortgroupdictwifvaluefilled

function sortgroupdictwifvalueempty(dict::Dict{Mode, Float64}, threshold::Float64)
    tupledata = sort([(value,key) for (key,value) in dict],by=first,rev=true)
    refvalue = tupledata[1][2]|>getattr(:eigenindex)
    result = []
    evals = []
    subresult = Subset(tupledata[1][2])
    for pair in tupledata
        # println(1-pair[1],threshold)
        if (1-pair[1])<threshold
            if pair[2]|>getattr(:eigenindex)==refvalue
                subresult = subresult + Subset(pair[2])
                append!(evals,pair[1])
            else
                append!(result,[subresult])
                refvalue = pair[2]|>getattr(:eigenindex)
                evals = []
                subresult = Subset(pair[2])
            end
        else 
            return result
        end
    end
    append!(result,[subresult])
    return result
end
export sortgroupdictwifvalueempty

function localwannierization(localbasis::FockMap, localseeds::FockMap, svdorthothreshold::Number = 1e-1)::FockMap
    U, Σ, Vt = svd(localbasis'*localseeds)
    minsvdvalue::Number = minimum(v for (_, v) in Σ)
    @info("min svdvalue", minsvdvalue)
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

function globalwannierfunction(correlations::FockMap, localisometry::FockMap)
    wanniercrystalisos = gmeracrystalisometries(localisometry=localisometry, crystalfock=correlations|>getoutspace,addinspacemomentuminfo=true)

    wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b)) for mode in localisometry|>getinspace |> orderedmodes)
    wanniercrystal::Crystal = Crystal(wannierunitcell, (correlations|>getoutspace|> getcrystal).sizes)

    blocks = Dict(((k, k)=>isometry) for (k, isometry) in wanniercrystalisos|>Dict)
    
    globalwannierizedfunction::FockMap = crystalfockmap(correlations|>getoutspace|>getcrystal, wanniercrystal, blocks)
    
    return globalwannierizedfunction
end
export globalwannierfunction

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

    isometries = paralleltasks(
        name="crystalisometries",
        tasks=(()->(k=>fouriermap[getsubspace(crystalfock, k), :]*preprocesslocalisometry(k)) for k in crystal|>brillouinzone),
        count=crystal|>vol)|>parallel

    return isometries
end
export gmeracrystalisometries

function firstgmerastep(correlations,threshold::Float64,nooffrozenrsmodes::Number,Asublatticeunitcellmodes::Subset{Mode},Bsublatticeunitcellmodes::Subset{Mode})
    crystalfock = correlations|>getoutspace
    crystal::Crystal = crystalfock|>getcrystal
    space::RealSpace = crystal|>getspace

    @info("Computing local correlations...")
    firstcenter = [0,0] ∈ space
    firsthexagonalregion = gethexagonalregion(crystal=crystal, center=firstcenter, metricspace=space)
    firsthexagonalregionfock = quantize(firsthexagonalregion,1)
    
    shiftedfirstcenter1 = [1,0] ∈ space
    shiftedfirsthexagonalregion1 = firsthexagonalregion.+shiftedfirstcenter1
    shiftedfirsthexagonalregion1fock = quantize(shiftedfirsthexagonalregion1,1)
    
    shiftedfirstcenter2 = [0,1] ∈ space
    shiftedfirsthexagonalregion2 = firsthexagonalregion.+shiftedfirstcenter2
    shiftedfirsthexagonalregion2fock = quantize(shiftedfirsthexagonalregion2,1)
    
    shiftedfirstcenter3 = [1,1] ∈ space
    shiftedfirsthexagonalregion3 = firsthexagonalregion.+shiftedfirstcenter3
    shiftedfirsthexagonalregion3fock = quantize(shiftedfirsthexagonalregion3,1)
    
    allregion = firsthexagonalregion+shiftedfirsthexagonalregion1+shiftedfirsthexagonalregion2+shiftedfirsthexagonalregion3
    if intersect(allregion,crystal|>getunitcell) == crystal|>getunitcell
        @info("cover the whole unitcell")
    else
        @error("the allregion cannot cover the whole unitcell ")
    end

    localcorrelations = regioncorrelations(correlations,firsthexagonalregionfock)
    localspectrum = localcorrelations|>eigspech
    display(localspectrum|>visualize)
    filleddict,emptydict = separatefilledandemptymodes(localspectrum)
    rankedandgroupedfilledmodes = sortgroupdictwifvalue(filleddict,false)
    rankedandgroupedemptymodes = sortgroupdictwifvalue(emptydict,true)

    rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-firstcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
    sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
    frozenrsmodes, courierrsmodes = findfrozenandcourierrsmode(sortedrsmode,nooffrozenrsmodes)
    filledrsmodes,emptyrsmodes = distinguishABsublatticemodesforhexagon(frozenrsmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)

    filledseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in filledrsmodes)]
    emptyseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in emptyrsmodes)]
    frozenseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
    courierseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]

    truncatedeigenvectofilled = rows((localspectrum|>geteigenvectors),FockSpace(filledrsmodes))
    truncatedeigenvectoempty = rows((localspectrum|>geteigenvectors),FockSpace(emptyrsmodes))

    # rankedandgroupedfilledmodeswifoverlap = [pair[2] for pair in sort([calculateavgoverlapofmodes(truncatedeigenvectofilled,modes) for modes in rankedandgroupedfilledmodes],by=first,rev=true)]
    # rankedandgroupedemptymodeswifoverlap = [pair[2] for pair in sort([calculateavgoverlapofmodes(truncatedeigenvectoempty,modes) for modes in rankedandgroupedemptymodes],by=first,rev=true)]

    chosenfilledmodes = choosemodesforwannierization(rankedandgroupedfilledmodes,truncatedeigenvectofilled,threshold,div(nooffrozenrsmodes,2))
    chosenemptymodes = choosemodesforwannierization(rankedandgroupedemptymodes,truncatedeigenvectoempty,threshold,div(nooffrozenrsmodes,2))
    chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
    chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

    localisofilled = columns(localspectrum |> geteigenvectors, FockSpace(chosenfilledmodes))
    localisoempty = columns(localspectrum |> geteigenvectors, FockSpace(chosenemptymodes))
    localisofrozen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))

    wannierfilled = localwannierization(localisofilled, filledseeds)
    wannierempty = localwannierization(localisoempty, emptyseeds)
    wannierfrozen = localwannierization(localisofrozen, frozenseeds)
    wanniercourier = localwannierization(localisocourier, courierseeds)

    (wannierfrozen'*localcorrelations*wannierfrozen)|>eigspec|>visualize|>display

    localwannierresults =  Dict(:wannierfilled => wannierfilled, :wannierempty => wannierempty,
                                :wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                                :filledseeds => filledseeds,:emptyseeds => emptyseeds,
                                :frozenseeds => frozenseeds,:courierseeds => courierseeds)

    wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

    firstshiftedhexagonalregionfocklist = [(shiftedfirstcenter1,shiftedfirsthexagonalregion1fock), (shiftedfirstcenter2,shiftedfirsthexagonalregion2fock), (shiftedfirstcenter3,shiftedfirsthexagonalregion3fock)]

    for (shifted,hexagonalregionfock) in firstshiftedhexagonalregionfocklist
        localcorrelations = regioncorrelations(correlations,hexagonalregionfock)
        localspectrum = localcorrelations|>eigspech
        filleddict,emptydict = separatefilledandemptymodes(localspectrum)
        rankedandgroupedfilledmodes = sortgroupdictwifvalue(filleddict,false)
        rankedandgroupedemptymodes = sortgroupdictwifvalue(emptydict,true)
        
        rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-(firstcenter+shifted))) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
        sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
        frozenrsmodes, courierrsmodes = findfrozenandcourierrsmode(sortedrsmode,nooffrozenrsmodes)
        filledrsmodes,emptyrsmodes = distinguishABsublatticemodesforhexagon(frozenrsmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)
    
        filledseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in filledrsmodes)]
        emptyseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in emptyrsmodes)]
        courierseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
        frozenseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
    
        truncatedeigenvectofilled = rows((localspectrum|>geteigenvectors),FockSpace(filledrsmodes))
        truncatedeigenvectoempty = rows((localspectrum|>geteigenvectors),FockSpace(emptyrsmodes))

        # rankedandgroupedfilledmodeswifoverlap = [pair[2] for pair in sort([calculateavgoverlapofmodes(truncatedeigenvectofilled,modes) for modes in rankedandgroupedfilledmodes],by=first,rev=true)]
        # rankedandgroupedemptymodeswifoverlap = [pair[2] for pair in sort([calculateavgoverlapofmodes(truncatedeigenvectoempty,modes) for modes in rankedandgroupedemptymodes],by=first,rev=true)]

        chosenfilledmodes = choosemodesforwannierization(rankedandgroupedfilledmodes,truncatedeigenvectofilled,threshold,div(nooffrozenrsmodes,2))
        chosenemptymodes = choosemodesforwannierization(rankedandgroupedemptymodes,truncatedeigenvectoempty,threshold,div(nooffrozenrsmodes,2))
        chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
        chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes
    
        localisofilled = columns(localspectrum |> geteigenvectors, FockSpace(chosenfilledmodes))
        localisoempty = columns(localspectrum |> geteigenvectors, FockSpace(chosenemptymodes))
        localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
        localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))
    
        wannierfilled = localwannierization(localisofilled, filledseeds)
        wannierempty = localwannierization(localisoempty, emptyseeds)
        wannierfrozen = localwannierization(localisoforzen, frozenseeds)
        wanniercourier = localwannierization(localisocourier, courierseeds)
        
        shiftedlocalwannierresults =  Dict(:wannierfilled => wannierfilled, :wannierempty => wannierempty,
                                    :wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                                    :filledseeds => filledseeds,:emptyseeds => emptyseeds,
                                    :frozenseeds => frozenseeds,:courierseeds => courierseeds)
    
        wannierinfos[hexagonalregionfock] =  shiftedlocalwannierresults
    end 

    firsthexagonalregionfocklist = [firsthexagonalregionfock, shiftedfirsthexagonalregion1fock, shiftedfirsthexagonalregion2fock, shiftedfirsthexagonalregion3fock]

    extendedwannierizedfrozen =  sum(wannierinfos[regionfock][:wannierfrozen] for regionfock in firsthexagonalregionfocklist)
    extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in firsthexagonalregionfocklist)
                
    origin = [0, 0] ∈ space
    refunictcellfockfrozen = FockSpace(Subset(mode for mode in extendedwannierizedfrozen |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

    wannierforzenisometry = globalwannierfunction(correlations,extendedwannierizedfrozen[:,refunictcellfockfrozen])
    wanniercourierisometry = globalwannierfunction(correlations,extendedwannierizedcourier[:,refunictcellfockcourier])

    couriercorrelations = wanniercourierisometry' * correlations * wanniercourierisometry
    couriercorrelations|>eigspech|>visualize|>display
    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

    firstrgedcorrelations = purifiedcouriercorrelations
    return firstrgedcorrelations
end
export firstgmerastep

function secondgmerastep(correlations,threshold::Float64,nooffrozenrsmodes::Number,Asublatticeunitcellmodes::Subset{Mode},Bsublatticeunitcellmodes::Subset{Mode})
    firstrgedcorrelations = correlations
    firstrgedcrystalfock = correlations|>getoutspace
    firstrgedcrystal::Crystal = firstrgedcrystalfock|>getcrystal
    firstrgedspace::RealSpace = firstrgedcrystal|>getspace

    @info("Computing local correlations...")
    secondcenter = [2/3,-1/3] ∈ firstrgedspace
    secondhexagonalregion = gethexagonalregion(crystal=firstrgedcrystal, center=secondcenter, metricspace=firstrgedspace)
    secondhexagonalregionfock = quantize(secondhexagonalregion,1)

    shiftedsecondcenter1 = [0,1] ∈ firstrgedspace
    shiftedsecondhexagonalregion1 = secondhexagonalregion.+shiftedsecondcenter1
    shiftedsecondhexagonalregion1fock = quantize(shiftedsecondhexagonalregion1,1)

    shiftedsecondcenter2 = [-1,1] ∈ firstrgedspace
    shiftedsecondhexagonalregion2 = secondhexagonalregion.+shiftedsecondcenter2
    shiftedsecondhexagonalregion2fock = quantize(shiftedsecondhexagonalregion2,1)

    allregion2 = secondhexagonalregion+shiftedsecondhexagonalregion1+shiftedsecondhexagonalregion2
    if intersect(allregion2,firstrgedcrystal|>getunitcell) == firstrgedcrystal|>getunitcell
        @info("cover the whole unitcell")
    else
        @error("the allregion cannot cover the whole unitcell ")
    end

    localcorrelations = regioncorrelations(firstrgedcorrelations,secondhexagonalregionfock)
    localspectrum = localcorrelations|>eigspech
    display(localspectrum|>visualize)
    filleddict,emptydict = separatefilledandemptymodes(localspectrum)
    rankedandgroupedfilledmodes = sortgroupdictwifvalue(filleddict,false)
    rankedandgroupedemptymodes = sortgroupdictwifvalue(emptydict,true)

    rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-secondcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
    sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
    frozenrsmodes, courierrsmodes = findfrozenandcourierrsmode(sortedrsmode,nooffrozenrsmodes)
    filledrsmodes,emptyrsmodes = distinguishABsublatticemodesforhexagon(frozenrsmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)

    filledseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in filledrsmodes)]
    emptyseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in emptyrsmodes)]
    frozenseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
    courierseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]

    truncatedeigenvectofilled = rows((localspectrum|>geteigenvectors),FockSpace(filledrsmodes))
    truncatedeigenvectoempty = rows((localspectrum|>geteigenvectors),FockSpace(emptyrsmodes))

    rankedandgroupedfilledmodeswifoverlap = [pair[2] for pair in sort([calculateavgoverlapofmodes(truncatedeigenvectofilled,modes) for modes in rankedandgroupedfilledmodes],by=first,rev=true)]
    rankedandgroupedemptymodeswifoverlap = [pair[2] for pair in sort([calculateavgoverlapofmodes(truncatedeigenvectoempty,modes) for modes in rankedandgroupedemptymodes],by=first,rev=true)]

    chosenfilledmodes = choosemodesforwannierization(rankedandgroupedfilledmodeswifoverlap,truncatedeigenvectofilled,threshold,div(nooffrozenrsmodes,2))
    chosenemptymodes = choosemodesforwannierization(rankedandgroupedemptymodeswifoverlap,truncatedeigenvectoempty,threshold,div(nooffrozenrsmodes,2))
    chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
    chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

    localisofilled = columns(localspectrum |> geteigenvectors, FockSpace(chosenfilledmodes))
    localisoempty = columns(localspectrum |> geteigenvectors, FockSpace(chosenemptymodes))
    localisofrozen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))

    wannierfilled = localwannierization(localisofilled, filledseeds)
    wannierempty = localwannierization(localisoempty, emptyseeds)
    wannierfrozen = localwannierization(localisofrozen, frozenseeds)
    wanniercourier = localwannierization(localisocourier, courierseeds)

    (wannierfrozen'*localcorrelations*wannierfrozen)|>eigspec|>visualize|>display

    localwannierresults =  Dict(:wannierfilled => wannierfilled, :wannierempty => wannierempty,
                                :wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                                :filledseeds => filledseeds,:emptyseeds => emptyseeds,
                                :frozenseeds => frozenseeds,:courierseeds => courierseeds)

    wannierinfos =  Dict(secondhexagonalregionfock=>localwannierresults)

    secondshiftedhexagonalregionfocklist = [(shiftedsecondcenter1,shiftedsecondhexagonalregion1fock), (shiftedsecondcenter2,shiftedsecondhexagonalregion2fock)]
                    
    for (shifted,hexagonalregionfock) in secondshiftedhexagonalregionfocklist
        localcorrelations = regioncorrelations(firstrgedcorrelations,hexagonalregionfock)
        localspectrum = localcorrelations|>eigspech
        filleddict,emptydict = separatefilledandemptymodes(localspectrum)
        rankedandgroupedfilledmodes = sortgroupdictwifvalue(filleddict,false)
        rankedandgroupedemptymodes = sortgroupdictwifvalue(emptydict,true)
                        
        rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-(secondcenter+shifted))) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
        sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
        frozenrsmodes, courierrsmodes = findfrozenandcourierrsmode(sortedrsmode,nooffrozenrsmodes)
        filledrsmodes,emptyrsmodes = distinguishABsublatticemodesforhexagon(frozenrsmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)

        filledseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in filledrsmodes)]
        emptyseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in emptyrsmodes)]
        courierseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
        frozenseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
                    
        truncatedeigenvectofilled = rows((localspectrum|>geteigenvectors),FockSpace(filledrsmodes))
        truncatedeigenvectoempty = rows((localspectrum|>geteigenvectors),FockSpace(emptyrsmodes))

        rankedandgroupedfilledmodeswifoverlap = [pair[2] for pair in sort([calculateavgoverlapofmodes(truncatedeigenvectofilled,modes) for modes in rankedandgroupedfilledmodes],by=first,rev=true)]
        rankedandgroupedemptymodeswifoverlap = [pair[2] for pair in sort([calculateavgoverlapofmodes(truncatedeigenvectoempty,modes) for modes in rankedandgroupedemptymodes],by=first,rev=true)]

        chosenfilledmodes = choosemodesforwannierization(rankedandgroupedfilledmodeswifoverlap,truncatedeigenvectofilled,threshold,div(nooffrozenrsmodes,2))
        chosenemptymodes = choosemodesforwannierization(rankedandgroupedemptymodeswifoverlap,truncatedeigenvectoempty,threshold,div(nooffrozenrsmodes,2))
        chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
        chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes
                        
        localisofilled = columns(localspectrum |> geteigenvectors, FockSpace(chosenfilledmodes))
        localisoempty = columns(localspectrum |> geteigenvectors, FockSpace(chosenemptymodes))
        localisofrozen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
        localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))

        wannierfilled = localwannierization(localisofilled, filledseeds)
        wannierempty = localwannierization(localisoempty, emptyseeds)
        wannierfrozen = localwannierization(localisofrozen, frozenseeds)
        wanniercourier = localwannierization(localisocourier, courierseeds)

        localwannierresults =  Dict(:wannierfilled => wannierfilled, :wannierempty => wannierempty,
                                    :wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                                    :filledseeds => filledseeds,:emptyseeds => emptyseeds,
                                    :frozenseeds => frozenseeds,:courierseeds => courierseeds)
                            
        shiftedlocalwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                                    :frozenseeds => frozenseeds,:courierseeds => courierseeds)
                    
        wannierinfos[hexagonalregionfock] =  shiftedlocalwannierresults
    end 
                    
    secondhexagonalregionfocklist = [secondhexagonalregionfock, shiftedsecondhexagonalregion1fock, shiftedsecondhexagonalregion2fock]
                    
    extendedwannierizedfrozen =  sum(wannierinfos[regionfock][:wannierfrozen] for regionfock in secondhexagonalregionfocklist)
    extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in secondhexagonalregionfocklist)
                                
    origin = [0, 0] ∈ firstrgedspace
    refunictcellfockfrozen = FockSpace(Subset(mode for mode in extendedwannierizedfrozen |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
                    
    wannierforzenisometry = globalwannierfunction(firstrgedcorrelations,extendedwannierizedfrozen[:,refunictcellfockfrozen])
    wanniercourierisometry = globalwannierfunction(firstrgedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
                    
    frozencorrelations = wannierforzenisometry' * firstrgedcorrelations * wannierforzenisometry              
    couriercorrelations2 = wanniercourierisometry' * firstrgedcorrelations * wanniercourierisometry
                    
    display(couriercorrelations2|>eigspech|>visualize)
    couriercorrelationspectrum2 = couriercorrelations2 |> crystalspectrum
    purifiedcorrelationspectrum2 = couriercorrelationspectrum2 |> roundingpurification
    purifiedcouriercorrelations2 = purifiedcorrelationspectrum2 |> CrystalFockMap
                    
    secondrgedcorrelations = purifiedcouriercorrelations2
    return secondrgedcorrelations
end 
export secondgmerastep

function thirdgmerastep(correlations,threshold::Float64,nooffrozenrsmodes::Number,Asublatticeunitcellmodes::Subset{Mode},Bsublatticeunitcellmodes::Subset{Mode})
    secondrgedcorrelations = correlations
    secondrgedcrystalfock = secondrgedcorrelations|>getoutspace
    secondrgedcrystal::Crystal = secondrgedcrystalfock|>getcrystal
    secondrgedspace::RealSpace = secondrgedcrystal|>getspace
    
    @info("Computing local correlations...")
    thirdcenter = [1/3,1/3] ∈ secondrgedspace
    thirdhexagonalregion = gethexagonalregion(crystal=secondrgedcrystal, center=thirdcenter, metricspace=secondrgedspace)
    thirdhexagonalregionfock = quantize(thirdhexagonalregion,1)

    allregion3 = thirdhexagonalregion

    if intersect(allregion3,secondrgedcrystal|>getunitcell) == secondrgedcrystal|>getunitcell
        @info("cover the whole unitcell")
    else
        @error("the allregion cannot cover the whole unitcell ")
    end

    localcorrelations = regioncorrelations(secondrgedcorrelations,thirdhexagonalregionfock)
    localspectrum = localcorrelations|>eigspech
    display(visualize(localspectrum))
    filleddict,emptydict = separatefilledandemptymodes(localspectrum)
    rankedandgroupedfilledmodes = sortgroupdictwifvalue(filleddict,false)
    rankedandgroupedemptymodes = sortgroupdictwifvalue(emptydict,true)

    rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-thirdcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
    sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
    frozenrsmodes, courierrsmodes = findfrozenandcourierrsmode(sortedrsmode,nooffrozenrsmodes)
    filledrsmodes,emptyrsmodes = distinguishABsublatticemodesforhexagon(frozenrsmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)
    
    filledseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in filledrsmodes)]
    emptyseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in emptyrsmodes)]
    courierseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
    frozenseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]

    truncatedeigenvectofilled = rows((localspectrum|>geteigenvectors),FockSpace(filledrsmodes))
    truncatedeigenvectoempty = rows((localspectrum|>geteigenvectors),FockSpace(emptyrsmodes))
                    
    chosenfilledmodes = choosemodesforwannierization(rankedandgroupedfilledmodes,truncatedeigenvectofilled,threshold,div(nooffrozenrsmodes,2))
    chosenemptymodes = choosemodesforwannierization(rankedandgroupedemptymodes,truncatedeigenvectoempty,threshold,div(nooffrozenrsmodes,2))
    chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
    chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

    localisofilled = columns(localspectrum |> geteigenvectors, FockSpace(chosenfilledmodes))
    localisoempty = columns(localspectrum |> geteigenvectors, FockSpace(chosenemptymodes))
    localisofrozen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))

    wannierfilled = localwannierization(localisofilled, filledseeds)
    wannierempty = localwannierization(localisoempty, emptyseeds)
    wannierfrozen = localwannierization(localisofrozen, frozenseeds)
    wanniercourier = localwannierization(localisocourier, courierseeds)

    (wannierfrozen'*localcorrelations*wannierfrozen)|>eigspec|>visualize|>display

    localwannierresults =  Dict(:wannierfilled => wannierfilled, :wannierempty => wannierempty,
                                :wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                                :filledseeds => filledseeds,:emptyseeds => emptyseeds,
                                :frozenseeds => frozenseeds,:courierseeds => courierseeds)

    wannierinfos =  Dict(thirdhexagonalregionfock=>localwannierresults)

    wannierforzenisometry = globalwannierfunction(secondrgedcorrelations,wannierinfos[thirdhexagonalregionfock][:wannierfrozen])
    wanniercourierisometry = globalwannierfunction(secondrgedcorrelations,wannierinfos[thirdhexagonalregionfock][:wanniercourier])

    frozencorrelations = wannierforzenisometry' * secondrgedcorrelations * wannierforzenisometry
    couriercorrelations3 = wanniercourierisometry' * secondrgedcorrelations * wanniercourierisometry

    couriercorrelationspectrum3 = couriercorrelations3 |> crystalspectrum
    display(couriercorrelations3|>eigspech|>visualize)
    purifiedcorrelationspectrum3 = couriercorrelationspectrum3 |> roundingpurification
    purifiedcouriercorrelations3 = purifiedcorrelationspectrum3 |> CrystalFockMap

    thirdrgedcorrelations = purifiedcouriercorrelations3
    return thirdrgedcorrelations
end 
export thirdgmerastep