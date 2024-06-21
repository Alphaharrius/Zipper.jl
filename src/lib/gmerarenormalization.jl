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

function rankedandgroupoffsets(hexagonalregion::Subset{Offset},chosen::Number)
    center = hexagonalregion|>getcenter
    offsetswifdist = sort([(offset,norm((offset|>euclidean)-(center|>euclidean))) for offset in hexagonalregion],by=x->x[2])
    ref = round(offsetswifdist[1][2],digits=5)
    result = []
    samedist = []
    for offsetwifdist in offsetswifdist
        if round(offsetwifdist[2],digits=5) == ref
            push!(samedist,offsetwifdist[1])
        else
            ref = round(offsetwifdist[2],digits=5)
            push!(result,samedist)
            samedist = []
            push!(samedist,offsetwifdist[1])
        end
    end
    return Subset([offset for offsets in result[1:chosen] for offset in offsets])
end
export rankedandgroupoffsets


function firstgmerastepbycount(blockedcorrelations,noofcouriermodes::Number)
    blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal
    blockedcrystalfock = blockedcorrelations|>getoutspace
    blockedspace::RealSpace = blockedcrystal|>getspace
    
    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    firstcenter = [0,0] ∈ blockedspace
    firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=firstcenter, metricspace=blockedspace)
    firsthexagonalregionfock = quantize(firsthexagonalregion,1)
    rgshiftedcenter1 = [2/3,-1/3] ∈ blockedspace
    firstrgshiftedhexagonalregion1 = firsthexagonalregion.+rgshiftedcenter1
    rgshiftedcenter2 = [1/3,1/3] ∈ blockedspace
    firstrgshiftedhexagonalregion2 = firsthexagonalregion.+rgshiftedcenter2

    siteAregion1 = intersect(intersect(firsthexagonalregion,firstrgshiftedhexagonalregion1 ),firstrgshiftedhexagonalregion2)
    siteAregionfock1 = quantize(siteAregion1,1)
    siteAregion2 = c3*siteAregion1
    siteAregionfock2 = quantize(siteAregion2,1)
    siteAregion3 = (c3)*siteAregion2
    siteAregionfock3 = quantize(siteAregion3,1)

    siteBregion1 = (c6)*siteAregion1
    siteBregionfock1 = quantize(siteBregion1,1)
    siteBregion2 = (c3)*siteBregion1
    siteBregionfock2 = quantize(siteBregion2,1)
    siteBregion3 = (c3)*siteBregion2
    siteBregionfock3 = quantize(siteBregion3,1)

    nooffrozenmodes = length(firsthexagonalregionfock|>orderedmodes)-noofcouriermodes
    nooffilledmodes = div(nooffrozenmodes,2)
    noofemptymodes = nooffilledmodes
    @info("no of filled modes = ",nooffilledmodes)
    @info("no of courier modes = ",noofcouriermodes)
    @info("no of empty modes = ",noofemptymodes)

    localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
    localspectrum = localcorrelations|>eigspec
    display(localspectrum|>visualize)
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofcouriermodes, noofemptymodes])
    localcourierisometry = localstates[2]|>FockMap
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[3]|>FockMap
    courierproj = localcourierisometry*localcourierisometry'

    nooflocalcourierseeds = div(noofcouriermodes,6)
    reduandancy = length((siteAregionfock1|>orderedmodes))-nooflocalcourierseeds

    localAcorrelations1 = columns(rows(courierproj,siteAregionfock1),siteAregionfock1)
    localAstates1 = getregionstates(localcorrelations=localAcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
    localAseeds1 = localAstates1[2]

    localAcorrelations2 = columns(rows(courierproj,siteAregionfock2),siteAregionfock2)
    localAstates2 = getregionstates(localcorrelations=localAcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
    localAseeds2 = localAstates2[2]


    localAcorrelations3 = columns(rows(courierproj,siteAregionfock3),siteAregionfock3)
    localAstates3 = getregionstates(localcorrelations=localAcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
    localAseeds3 = localAstates3[2]

    allAseedsstate = localAseeds1+localAseeds2+localAseeds3

    localBcorrelations1 = columns(rows(courierproj,siteBregionfock1),siteBregionfock1)
    localBstates1 = getregionstates(localcorrelations=localBcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
    localBseeds1 = localBstates1[2]

    localBcorrelations2 = columns(rows(courierproj,siteBregionfock2),siteBregionfock2)
    localBstates2 = getregionstates(localcorrelations=localBcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
    localBseeds2 = localBstates2[2]

    localBcorrelations3 = columns(rows(courierproj,siteBregionfock3),siteBregionfock3)
    localBstates3 = getregionstates(localcorrelations=localBcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
    localBseeds3 = localBstates3[2]

    allBseedsstate = localBseeds1+localBseeds2+localBseeds3

    allseedsstates = (allAseedsstate+allBseedsstate)
    courierseeds = allseedsstates|>FockMap

    wanniercourier = localwannierization(localcourierisometry, courierseeds)
    display((wanniercourier'*localcorrelations*wanniercourier)|>eigspec|>visualize)

    localwannierresults =  Dict(:wanniercourier => wanniercourier,:localfilledisometry => localfilledisometry,
                            :localemptyisometry => localemptyisometry, :courierseeds => courierseeds)

    wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

    shiftedfirstcenter1 = [1,0] ∈ blockedspace
    shiftedfirstcenter2 = [0,1] ∈ blockedspace
    shiftedfirstcenter3 = [1,1] ∈ blockedspace
    
    firstshiftedhexagonalcenters = [shiftedfirstcenter1, shiftedfirstcenter2, shiftedfirstcenter3]
    
    for hexagonalcenter in firstshiftedhexagonalcenters
        shiftedfirsthexagonalregion = firsthexagonalregion.+hexagonalcenter
        shiftedfirsthexagonalregionfock = quantize(shiftedfirsthexagonalregion,1)
        c6recenter = recenter(c6,hexagonalcenter)
        c3recenter = recenter(c3,hexagonalcenter)
        
        shiftedsiteAregion1 = siteAregion1.+hexagonalcenter
        shiftedsiteAregionfock1 = quantize(shiftedsiteAregion1,1)
        shiftedsiteAregion2 = (c3recenter)*shiftedsiteAregion1
        shiftedsiteAregionfock2 = quantize(shiftedsiteAregion2,1)
        shiftedsiteAregion3 = (c3recenter)*shiftedsiteAregion2
        shiftedsiteAregionfock3 = quantize(shiftedsiteAregion3,1)

        shiftedsiteBregion1 = (c6recenter)*shiftedsiteAregion1
        shiftedsiteBregionfock1 = quantize(shiftedsiteBregion1,1)
        shiftedsiteBregion2 = (c3recenter)*shiftedsiteBregion1
        shiftedsiteBregionfock2 = quantize(shiftedsiteBregion2,1)
        shiftedsiteBregion3 = (c3recenter)*shiftedsiteBregion2
        shiftedsiteBregionfock3 = quantize(shiftedsiteBregion3,1)

        shiftedlocalcorrelations = regioncorrelations(blockedcorrelations,shiftedfirsthexagonalregionfock)
        shiftedlocalstates = getregionstates(localcorrelations=shiftedlocalcorrelations, grouping=[nooffilledmodes, noofcouriermodes, noofemptymodes])
        shiftedlocalcourierisometry = shiftedlocalstates[2]|>FockMap
        shiftedlocalfilledisometry = shiftedlocalstates[1]|>FockMap
        shiftedlocalemptyisometry = shiftedlocalstates[3]|>FockMap
        shiftedcourierproj = shiftedlocalcourierisometry*shiftedlocalcourierisometry'

        shiftedlocalAcorrelations1 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock1),shiftedsiteAregionfock1)
        shiftedlocalAstates1 = getregionstates(localcorrelations=shiftedlocalAcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalAseeds1 = shiftedlocalAstates1[2]

        shiftedlocalAcorrelations2 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock2),shiftedsiteAregionfock2)
        shiftedlocalAstates2 = getregionstates(localcorrelations=shiftedlocalAcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalAseeds2 = shiftedlocalAstates2[2]

        shiftedlocalAcorrelations3 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock3),shiftedsiteAregionfock3)
        shiftedlocalAstates3 = getregionstates(localcorrelations=shiftedlocalAcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalAseeds3 = shiftedlocalAstates3[2]

        shiftedallAseedsstate = shiftedlocalAseeds1+shiftedlocalAseeds2+shiftedlocalAseeds3

        shiftedlocalBcorrelations1 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock1),shiftedsiteBregionfock1)
        shiftedlocalBstates1 = getregionstates(localcorrelations=shiftedlocalBcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalBseeds1 = shiftedlocalBstates1[2]

        shiftedlocalBcorrelations2 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock2),shiftedsiteBregionfock2)
        shiftedlocalBstates2 = getregionstates(localcorrelations=shiftedlocalBcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalBseeds2 = shiftedlocalBstates2[2]

        shiftedlocalBcorrelations3 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock3),shiftedsiteBregionfock3)
        shiftedlocalBstates3 = getregionstates(localcorrelations=shiftedlocalBcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalBseeds3 = shiftedlocalBstates3[2]

        shiftedallBseedsstate = shiftedlocalBseeds1+shiftedlocalBseeds2+shiftedlocalBseeds3

        shiftedallseedsstates = (shiftedallAseedsstate+shiftedallBseedsstate)
        shiftedcourierseeds = shiftedallseedsstates|>FockMap

        shiftedwanniercourier = localwannierization(shiftedlocalcourierisometry, shiftedcourierseeds)

        shiftedlocalwannierresults =  Dict(:wanniercourier => shiftedwanniercourier, :localfilledisometry => shiftedlocalfilledisometry,
                                        :localemptyisometry => shiftedlocalemptyisometry, :courierseeds => shiftedcourierseeds)
        wannierinfos[shiftedfirsthexagonalregionfock] =  shiftedlocalwannierresults
    end

    ref = [quantize(firsthexagonalregion.+hexagonalcenter,1) for hexagonalcenter in firstshiftedhexagonalcenters]
    firsthexagonalregionfocklist = [firsthexagonalregionfock,ref...]

    extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in firsthexagonalregionfocklist)
    extendedwannierizedfilled =  sum(wannierinfos[regionfock][:localfilledisometry] for regionfock in firsthexagonalregionfocklist)
    extendedwannierizedempty =  sum(wannierinfos[regionfock][:localemptyisometry] for regionfock in firsthexagonalregionfocklist)
 
                
    origin = [0, 0] ∈ blockedspace
    refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

    wanniercourierisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
    globalisometryfilled = globalwannierfunction(blockedcorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
    globalisometryempty = globalwannierfunction(blockedcorrelations,extendedwannierizedempty[:,refunictcellfockempty])

    couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
    # display(couriercorrelations|>eigspech|>visualize)

    filledcorrelations = globalisometryfilled' * blockedcorrelations * globalisometryfilled
    # display(filledcorrelations|>eigspech|>visualize)

    couriercorrelationspectrum = couriercorrelations |> crystalspectrum

    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

    firstgmeracorrelations = purifiedcouriercorrelations
    return firstgmeracorrelations,globalisometryfilled,globalisometryempty,wanniercourierisometry
end
export firstgmerastepbycount

function firstgmerastepbythreshold(blockedcorrelations,threshold::Float64)
    blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal
    blockedcrystalfock = blockedcorrelations|>getoutspace
    blockedspace::RealSpace = blockedcrystal|>getspace
    
    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    firstcenter = [0,0] ∈ blockedspace
    firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=firstcenter, metricspace=blockedspace)
    firsthexagonalregionfock = quantize(firsthexagonalregion,1)
    rgshiftedcenter1 = [2/3,-1/3] ∈ blockedspace
    firstrgshiftedhexagonalregion1 = firsthexagonalregion.+rgshiftedcenter1
    rgshiftedcenter2 = [1/3,1/3] ∈ blockedspace
    firstrgshiftedhexagonalregion2 = firsthexagonalregion.+rgshiftedcenter2

    siteAregion1 = intersect(intersect(firsthexagonalregion,firstrgshiftedhexagonalregion1 ),firstrgshiftedhexagonalregion2)
    siteAregionfock1 = quantize(siteAregion1,1)
    siteAregion2 = c3*siteAregion1
    siteAregionfock2 = quantize(siteAregion2,1)
    siteAregion3 = (c3)*siteAregion2
    siteAregionfock3 = quantize(siteAregion3,1)

    siteBregion1 = (c6)*siteAregion1
    siteBregionfock1 = quantize(siteBregion1,1)
    siteBregion2 = (c3)*siteBregion1
    siteBregionfock2 = quantize(siteBregion2,1)
    siteBregion3 = (c3)*siteBregion2
    siteBregionfock3 = quantize(siteBregion3,1)

    localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
    localspectrum = localcorrelations|>eigspec
    display(localspectrum|>visualize)

    reffilledmodes = Subset([m for (m, v) in localspectrum|>geteigenvalues if norm(v) < threshold])
    refemptymodes = Subset([m for (m, v) in localspectrum|>geteigenvalues if norm(v) > 1-threshold])

    if length(reffilledmodes)%3==0
        nooffilledmodes = length(reffilledmodes)
    else
        nooffilledmodes = length(reffilledmodes)-(length(reffilledmodes)%3)
    end
    # if length(refemptymodes)%3==0
    #     noofemptymodes = length(refemptymodes)
    # else
    #     noofemptymodes = length(refemptymodes)-(length(refemptymodes)%3)
    # end    
    noofemptymodes = nooffilledmodes
    nooffrozenmodes = nooffilledmodes + noofemptymodes
    noofcouriermodes = length(firsthexagonalregionfock|>orderedmodes)-nooffrozenmodes
    @info("no of filled modes = ",nooffilledmodes)
    @info("no of courier modes = ",noofcouriermodes)
    @info("no of empty modes = ",noofemptymodes)

    rankedmodeswifeval = sort([(m,norm(v)) for (m, v) in localspectrum|>geteigenvalues],by=x->x[2])
    pickedfilledeval = rankedmodeswifeval[nooffilledmodes][2]
    pickedemptyeval = rankedmodeswifeval[end-noofemptymodes+1][2]
        
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofcouriermodes, noofemptymodes])
    localcourierisometry = localstates[2]|>FockMap
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[3]|>FockMap
    courierproj = localcourierisometry*localcourierisometry'

    nooflocalcourierseeds = div(noofcouriermodes,6)
    reduandancy = length((siteAregionfock1|>orderedmodes))-nooflocalcourierseeds

    localAcorrelations1 = columns(rows(courierproj,siteAregionfock1),siteAregionfock1)
    localAstates1 = getregionstates(localcorrelations=localAcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
    localAseeds1 = localAstates1[2]

    localAcorrelations2 = columns(rows(courierproj,siteAregionfock2),siteAregionfock2)
    localAstates2 = getregionstates(localcorrelations=localAcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
    localAseeds2 = localAstates2[2]


    localAcorrelations3 = columns(rows(courierproj,siteAregionfock3),siteAregionfock3)
    localAstates3 = getregionstates(localcorrelations=localAcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
    localAseeds3 = localAstates3[2]

    allAseedsstate = localAseeds1+localAseeds2+localAseeds3

    localBcorrelations1 = columns(rows(courierproj,siteBregionfock1),siteBregionfock1)
    localBstates1 = getregionstates(localcorrelations=localBcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
    localBseeds1 = localBstates1[2]

    localBcorrelations2 = columns(rows(courierproj,siteBregionfock2),siteBregionfock2)
    localBstates2 = getregionstates(localcorrelations=localBcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
    localBseeds2 = localBstates2[2]

    localBcorrelations3 = columns(rows(courierproj,siteBregionfock3),siteBregionfock3)
    localBstates3 = getregionstates(localcorrelations=localBcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
    localBseeds3 = localBstates3[2]

    allBseedsstate = localBseeds1+localBseeds2+localBseeds3

    allseedsstates = (allAseedsstate+allBseedsstate)
    courierseeds = allseedsstates|>FockMap

    @info("wannierizing local courier state")
    wanniercourier = localwannierization(localcourierisometry, courierseeds)
    display((wanniercourier'*localcorrelations*wanniercourier)|>eigspec|>visualize)

    localwannierresults =  Dict(:wanniercourier => wanniercourier,:localfilledisometry => localfilledisometry,
                            :localemptyisometry => localemptyisometry, :courierseeds => courierseeds)

    wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

    shiftedfirstcenter1 = [1,0] ∈ blockedspace
    shiftedfirstcenter2 = [0,1] ∈ blockedspace
    shiftedfirstcenter3 = [1,1] ∈ blockedspace
    
    firstshiftedhexagonalcenters = [shiftedfirstcenter1, shiftedfirstcenter2, shiftedfirstcenter3]
    
    for hexagonalcenter in firstshiftedhexagonalcenters
        shiftedfirsthexagonalregion = firsthexagonalregion.+hexagonalcenter
        shiftedfirsthexagonalregionfock = quantize(shiftedfirsthexagonalregion,1)
        c6recenter = recenter(c6,hexagonalcenter)
        c3recenter = recenter(c3,hexagonalcenter)
        
        shiftedsiteAregion1 = siteAregion1.+hexagonalcenter
        shiftedsiteAregionfock1 = quantize(shiftedsiteAregion1,1)
        shiftedsiteAregion2 = (c3recenter)*shiftedsiteAregion1
        shiftedsiteAregionfock2 = quantize(shiftedsiteAregion2,1)
        shiftedsiteAregion3 = (c3recenter)*shiftedsiteAregion2
        shiftedsiteAregionfock3 = quantize(shiftedsiteAregion3,1)

        shiftedsiteBregion1 = (c6recenter)*shiftedsiteAregion1
        shiftedsiteBregionfock1 = quantize(shiftedsiteBregion1,1)
        shiftedsiteBregion2 = (c3recenter)*shiftedsiteBregion1
        shiftedsiteBregionfock2 = quantize(shiftedsiteBregion2,1)
        shiftedsiteBregion3 = (c3recenter)*shiftedsiteBregion2
        shiftedsiteBregionfock3 = quantize(shiftedsiteBregion3,1)

        shiftedlocalcorrelations = regioncorrelations(blockedcorrelations,shiftedfirsthexagonalregionfock)
        shiftedlocalstates = getregionstates(localcorrelations=shiftedlocalcorrelations, grouping=[nooffilledmodes, noofcouriermodes, noofemptymodes])
        shiftedlocalcourierisometry = shiftedlocalstates[2]|>FockMap
        shiftedlocalfilledisometry = shiftedlocalstates[1]|>FockMap
        shiftedlocalemptyisometry = shiftedlocalstates[3]|>FockMap
        shiftedcourierproj = shiftedlocalcourierisometry*shiftedlocalcourierisometry'

        shiftedlocalAcorrelations1 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock1),shiftedsiteAregionfock1)
        shiftedlocalAstates1 = getregionstates(localcorrelations=shiftedlocalAcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalAseeds1 = shiftedlocalAstates1[2]

        shiftedlocalAcorrelations2 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock2),shiftedsiteAregionfock2)
        shiftedlocalAstates2 = getregionstates(localcorrelations=shiftedlocalAcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalAseeds2 = shiftedlocalAstates2[2]

        shiftedlocalAcorrelations3 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock3),shiftedsiteAregionfock3)
        shiftedlocalAstates3 = getregionstates(localcorrelations=shiftedlocalAcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalAseeds3 = shiftedlocalAstates3[2]

        shiftedallAseedsstate = shiftedlocalAseeds1+shiftedlocalAseeds2+shiftedlocalAseeds3

        shiftedlocalBcorrelations1 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock1),shiftedsiteBregionfock1)
        shiftedlocalBstates1 = getregionstates(localcorrelations=shiftedlocalBcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalBseeds1 = shiftedlocalBstates1[2]

        shiftedlocalBcorrelations2 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock2),shiftedsiteBregionfock2)
        shiftedlocalBstates2 = getregionstates(localcorrelations=shiftedlocalBcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalBseeds2 = shiftedlocalBstates2[2]

        shiftedlocalBcorrelations3 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock3),shiftedsiteBregionfock3)
        shiftedlocalBstates3 = getregionstates(localcorrelations=shiftedlocalBcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalBseeds3 = shiftedlocalBstates3[2]

        shiftedallBseedsstate = shiftedlocalBseeds1+shiftedlocalBseeds2+shiftedlocalBseeds3

        shiftedallseedsstates = (shiftedallAseedsstate+shiftedallBseedsstate)
        shiftedcourierseeds = shiftedallseedsstates|>FockMap

        @info("wannierizing local courier state")
        shiftedwanniercourier = localwannierization(shiftedlocalcourierisometry, shiftedcourierseeds)

        shiftedlocalwannierresults =  Dict(:wanniercourier => shiftedwanniercourier, :localfilledisometry => shiftedlocalfilledisometry,
                                        :localemptyisometry => shiftedlocalemptyisometry, :courierseeds => shiftedcourierseeds)
        wannierinfos[shiftedfirsthexagonalregionfock] =  shiftedlocalwannierresults
    end

    ref = [quantize(firsthexagonalregion.+hexagonalcenter,1) for hexagonalcenter in firstshiftedhexagonalcenters]
    firsthexagonalregionfocklist = [firsthexagonalregionfock,ref...]

    extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in firsthexagonalregionfocklist)
    extendedwannierizedfilled =  sum(wannierinfos[regionfock][:localfilledisometry] for regionfock in firsthexagonalregionfocklist)
    extendedwannierizedempty =  sum(wannierinfos[regionfock][:localemptyisometry] for regionfock in firsthexagonalregionfocklist)
 
                
    origin = [0, 0] ∈ blockedspace
    refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

    wanniercourierisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
    globalfilledisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
    globalemptyisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedempty[:,refunictcellfockempty])


    @info "Computing local courier states..."
    leftrestrict = fourier(wanniercourierisometry|>getoutspace, firsthexagonalregionfock) / (blockedcrystal|>vol|>sqrt)
    rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
    wanniercourierstate = leftrestrict' * wanniercourierisometry * rightrestrict

    couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
    # display(couriercorrelations|>eigspech|>visualize)

    filledcorrelations = globalfilledisometry' * blockedcorrelations * globalfilledisometry
    # display(filledcorrelations|>eigspech|>visualize)

    emptycorrelations = globalemptyisometry' * blockedcorrelations * globalemptyisometry

    couriercorrelationspectrum = couriercorrelations |> crystalspectrum

    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

    return Dict(
        :couriercorrelations=>purifiedcouriercorrelations,
        :filledcorrelations=>filledcorrelations,
        :emptycorrelations=>emptycorrelations,
        :wanniercourierisometry=>wanniercourierisometry,
        :globalfilledisometry=>globalfilledisometry,
        :globalemptyisometry=>globalemptyisometry,
        :pickedfilledeval=>pickedfilledeval,
        :pickedemptyeval=>pickedemptyeval,
        :nooflocalcouriermodes=>noofcouriermodes,
        :wanniercourierstates=>wanniercourierstate|>RegionState,
        :rawcouriercorrelations=>couriercorrelations)
end
export firstgmerastepbythreshold


function secondgmerastepbycount(firstgmeracorrelations,noofcouriermodes::Number,noofflavourpermode::Number)
    firstgmeracrystalfock = firstgmeracorrelations|>getoutspace
    firstgmeracrystal::Crystal = firstgmeracrystalfock|>getcrystal
    firstgmeraspace::RealSpace = firstgmeracrystal|>getspace

    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    secondcenter = [2/3,-1/3] ∈ firstgmeraspace
    secondhexagonalregion = gethexagonalregion(rot=refrot,crystal=firstgmeracrystal, center=secondcenter, metricspace=firstgmeraspace)
    secondhexagonalregionfock = quantize(secondhexagonalregion,noofflavourpermode)
    rgshiftedcenter1 = [2/3,-1/3] ∈ firstgmeraspace
    secondrgshiftedhexagonalregion1 = secondhexagonalregion.+rgshiftedcenter1
    rgshiftedcenter2 = [1/3,1/3] ∈ firstgmeraspace
    secondrgshiftedhexagonalregion2 = secondhexagonalregion.+rgshiftedcenter2

    c6recenter = recenter(c6,secondcenter)
    c3recenter = recenter(c3,secondcenter)

    siteAregion1 = intersect(intersect(secondhexagonalregion,secondrgshiftedhexagonalregion1 ),secondrgshiftedhexagonalregion2)
    siteAregionfock1 = quantize(siteAregion1,noofflavourpermode)
    siteAregion2 = (c3recenter)*siteAregion1
    siteAregionfock2 = quantize(siteAregion2,noofflavourpermode)
    siteAregion3 = (c3recenter)*siteAregion2
    siteAregionfock3 = quantize(siteAregion3,noofflavourpermode)

    siteBregion1 = (c6recenter)*siteAregion1
    siteBregionfock1 = quantize(siteBregion1,noofflavourpermode)
    siteBregion2 = (c3recenter)*siteBregion1
    siteBregionfock2 = quantize(siteBregion2,noofflavourpermode)
    siteBregion3 = (c3recenter)*siteBregion2
    siteBregionfock3 = quantize(siteBregion3,noofflavourpermode)

    nooffrozenmodes = length(secondhexagonalregionfock|>orderedmodes)-noofcouriermodes
    nooffilledmodes = div(nooffrozenmodes,2)
    noofemptymodes = nooffilledmodes
    @info("no of filled modes = ",nooffilledmodes)
    @info("no of courier modes = ",noofcouriermodes)
    @info("no of empty modes = ",noofemptymodes)

    localcorrelations = regioncorrelations(firstgmeracorrelations,secondhexagonalregionfock)
    localspectrum = localcorrelations|>eigspech
    display(localspectrum|>visualize)
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofcouriermodes, noofemptymodes])
    localcourierisometry = localstates[2]|>FockMap
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[3]|>FockMap
    courierproj = localcourierisometry*localcourierisometry'
    nooflocalcourierseeds = div(noofcouriermodes,6)
    reduandancy = length((siteAregionfock1|>orderedmodes))-nooflocalcourierseeds

    localAcorrelations1 = columns(rows(courierproj,siteAregionfock1),siteAregionfock1)
    localAstates1 = getregionstates(localcorrelations=localAcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
    localAseeds1 = localAstates1[2]

    localAcorrelations2 = columns(rows(courierproj,siteAregionfock2),siteAregionfock2)
    localAstates2 = getregionstates(localcorrelations=localAcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
    localAseeds2 = localAstates2[2]


    localAcorrelations3 = columns(rows(courierproj,siteAregionfock3),siteAregionfock3)
    localAstates3 = getregionstates(localcorrelations=localAcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
    localAseeds3 = localAstates3[2]

    allAseedsstate = localAseeds1+localAseeds2+localAseeds3

    localBcorrelations1 = columns(rows(courierproj,siteBregionfock1),siteBregionfock1)
    localBstates1 = getregionstates(localcorrelations=localBcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
    localBseeds1 = localBstates1[2]

    localBcorrelations2 = columns(rows(courierproj,siteBregionfock2),siteBregionfock2)
    localBstates2 = getregionstates(localcorrelations=localBcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
    localBseeds2 = localBstates2[2]

    localBcorrelations3 = columns(rows(courierproj,siteBregionfock3),siteBregionfock3)
    localBstates3 = getregionstates(localcorrelations=localBcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
    localBseeds3 = localBstates3[2]

    allBseedsstate = localBseeds1+localBseeds2+localBseeds3

    allseedsstates = (allAseedsstate+allBseedsstate)
    courierseeds = allseedsstates|>FockMap
    localcourierisometry

    wanniercourier = localwannierization(localcourierisometry, courierseeds)
    display((wanniercourier'*localcorrelations*wanniercourier)|>eigspec|>visualize)

    localwannierresults =  Dict(:wanniercourier => wanniercourier, :localfilledisometry => localfilledisometry,
                                :localemptyisometry => localemptyisometry, :courierseeds => courierseeds)

    wannierinfos =  Dict(secondhexagonalregionfock=>localwannierresults)

    shiftedsecondcenter1 = [0,1] ∈ firstgmeraspace
    shiftedsecondcenter2 = [-1,1] ∈ firstgmeraspace

    secondshiftedhexagonalcenterlist = [shiftedsecondcenter1, shiftedsecondcenter2]

    for hexagonalcenter in secondshiftedhexagonalcenterlist
        shiftedsecondhexagonalregion = secondhexagonalregion.+hexagonalcenter
        shiftedsecondhexagonalregionfock = quantize(shiftedsecondhexagonalregion,noofflavourpermode)
        c6recenter = recenter(c6,secondcenter+hexagonalcenter)
        c3recenter = recenter(c3,secondcenter+hexagonalcenter)
    
        shiftedsiteAregion1 = siteAregion1.+hexagonalcenter
        shiftedsiteAregionfock1 = quantize(shiftedsiteAregion1,noofflavourpermode)
        shiftedsiteAregion2 = (c3recenter)*shiftedsiteAregion1
        shiftedsiteAregionfock2 = quantize(shiftedsiteAregion2,noofflavourpermode)
        shiftedsiteAregion3 = (c3recenter)*shiftedsiteAregion2
        shiftedsiteAregionfock3 = quantize(shiftedsiteAregion3,noofflavourpermode)
    
        shiftedsiteBregion1 = (c6recenter)*shiftedsiteAregion1
        shiftedsiteBregionfock1 = quantize(shiftedsiteBregion1,noofflavourpermode)
        shiftedsiteBregion2 = (c3recenter)*shiftedsiteBregion1
        shiftedsiteBregionfock2 = quantize(shiftedsiteBregion2,noofflavourpermode)
        shiftedsiteBregion3 = (c3recenter)*shiftedsiteBregion2
        shiftedsiteBregionfock3 = quantize(shiftedsiteBregion3,noofflavourpermode)
    
        shiftedlocalcorrelations = regioncorrelations(firstgmeracorrelations,shiftedsecondhexagonalregionfock)
        shiftedlocalstates = getregionstates(localcorrelations=shiftedlocalcorrelations, grouping=[nooffilledmodes, noofcouriermodes, noofemptymodes])
        shiftedlocalcourierisometry = shiftedlocalstates[2]|>FockMap
        shiftedlocalfilledisometry = shiftedlocalstates[1]|>FockMap
        shiftedlocalemptyisometry = shiftedlocalstates[3]|>FockMap
        shiftedcourierproj = shiftedlocalcourierisometry*shiftedlocalcourierisometry'
    
        shiftedlocalAcorrelations1 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock1),shiftedsiteAregionfock1)
        shiftedlocalAstates1 = getregionstates(localcorrelations=shiftedlocalAcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalAseeds1 = shiftedlocalAstates1[2]
    
        shiftedlocalAcorrelations2 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock2),shiftedsiteAregionfock2)
        shiftedlocalAstates2 = getregionstates(localcorrelations=shiftedlocalAcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalAseeds2 = shiftedlocalAstates2[2]
    
        shiftedlocalAcorrelations3 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock3),shiftedsiteAregionfock3)
        shiftedlocalAstates3 = getregionstates(localcorrelations=shiftedlocalAcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalAseeds3 = shiftedlocalAstates3[2]
    
        shiftedallAseedsstate = shiftedlocalAseeds1+shiftedlocalAseeds2+shiftedlocalAseeds3
    
        shiftedlocalBcorrelations1 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock1),shiftedsiteBregionfock1)
        shiftedlocalBstates1 = getregionstates(localcorrelations=shiftedlocalBcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalBseeds1 = shiftedlocalBstates1[2]
    
        shiftedlocalBcorrelations2 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock2),shiftedsiteBregionfock2)
        shiftedlocalBstates2 = getregionstates(localcorrelations=shiftedlocalBcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalBseeds2 = shiftedlocalBstates2[2]
    
        shiftedlocalBcorrelations3 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock3),shiftedsiteBregionfock3)
        shiftedlocalBstates3 = getregionstates(localcorrelations=shiftedlocalBcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalBseeds3 = shiftedlocalBstates3[2]
    
        shiftedallBseedsstate = shiftedlocalBseeds1+shiftedlocalBseeds2+shiftedlocalBseeds3
    
        shiftedallseedsstates = (shiftedallAseedsstate+shiftedallBseedsstate)
        shiftedcourierseeds = shiftedallseedsstates|>FockMap
    
        shiftedwanniercourier = localwannierization(shiftedlocalcourierisometry, shiftedcourierseeds)
    
        shiftedlocalwannierresults = Dict(:wanniercourier => shiftedwanniercourier, :localfilledisometry => shiftedlocalfilledisometry,
                                        :localemptyisometry => shiftedlocalemptyisometry, :courierseeds => shiftedcourierseeds)
        wannierinfos[shiftedsecondhexagonalregionfock] =  shiftedlocalwannierresults
    end

    ref = [quantize(secondhexagonalregion.+hexagonalcenter,noofflavourpermode) for hexagonalcenter in secondshiftedhexagonalcenterlist]
    secondhexagonalregionfocklist = [secondhexagonalregionfock,ref...]

    extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in secondhexagonalregionfocklist)
    extendedwannierizedfilled =  sum(wannierinfos[regionfock][:localfilledisometry] for regionfock in secondhexagonalregionfocklist)
    extendedwannierizedempty =  sum(wannierinfos[regionfock][:localemptyisometry] for regionfock in secondhexagonalregionfocklist) 

    origin = [0, 0] ∈ firstgmeraspace
    refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

    wanniercourierisometry = globalwannierfunction(firstgmeracorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
    globalisometryfilled = globalwannierfunction(firstgmeracorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
    globalisometryempty = globalwannierfunction(firstgmeracorrelations,extendedwannierizedempty[:,refunictcellfockempty])

    couriercorrelations = wanniercourierisometry' * firstgmeracorrelations * wanniercourierisometry
    # display(couriercorrelations|>eigspech|>visualize)

    filledcorrelations = globalisometryfilled' * firstgmeracorrelations * globalisometryfilled
    # display(filledcorrelations|>eigspech|>visualize)

    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

    secondgmeracorrelations = purifiedcouriercorrelations
    return secondgmeracorrelations,globalisometryfilled,globalisometryempty,wanniercourierisometry
end
export secondgmerastepbycount

function secondgmerastepbythreshold(firstgmeracorrelations,threshold::Float64,noofflavourpermode::Number)
    firstgmeracrystalfock = firstgmeracorrelations|>getoutspace
    firstgmeracrystal::Crystal = firstgmeracrystalfock|>getcrystal
    firstgmeraspace::RealSpace = firstgmeracrystal|>getspace

    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    secondcenter = [2/3,-1/3] ∈ firstgmeraspace
    secondhexagonalregion = gethexagonalregion(rot=refrot,crystal=firstgmeracrystal, center=secondcenter, metricspace=firstgmeraspace)
    secondhexagonalregionfock = quantize(secondhexagonalregion,noofflavourpermode)
    rgshiftedcenter1 = [2/3,-1/3] ∈ firstgmeraspace
    secondrgshiftedhexagonalregion1 = secondhexagonalregion.+rgshiftedcenter1
    rgshiftedcenter2 = [1/3,1/3] ∈ firstgmeraspace
    secondrgshiftedhexagonalregion2 = secondhexagonalregion.+rgshiftedcenter2

    c6recenter = recenter(c6,secondcenter)
    c3recenter = recenter(c3,secondcenter)

    siteAregion1 = intersect(intersect(secondhexagonalregion,secondrgshiftedhexagonalregion1 ),secondrgshiftedhexagonalregion2)
    siteAregionfock1 = quantize(siteAregion1,noofflavourpermode)
    siteAregion2 = (c3recenter)*siteAregion1
    siteAregionfock2 = quantize(siteAregion2,noofflavourpermode)
    siteAregion3 = (c3recenter)*siteAregion2
    siteAregionfock3 = quantize(siteAregion3,noofflavourpermode)

    siteBregion1 = (c6recenter)*siteAregion1
    siteBregionfock1 = quantize(siteBregion1,noofflavourpermode)
    siteBregion2 = (c3recenter)*siteBregion1
    siteBregionfock2 = quantize(siteBregion2,noofflavourpermode)
    siteBregion3 = (c3recenter)*siteBregion2
    siteBregionfock3 = quantize(siteBregion3,noofflavourpermode)

    localcorrelations = regioncorrelations(firstgmeracorrelations,secondhexagonalregionfock)
    localspectrum = localcorrelations|>eigspec
    display(localspectrum|>visualize)

    reffilledmodes = Subset([m for (m, v) in localspectrum|>geteigenvalues if norm(v) < threshold])
    refemptymodes = Subset([m for (m, v) in localspectrum|>geteigenvalues if norm(v) > 1-threshold])

    if length(reffilledmodes)%3==0
        nooffilledmodes = length(reffilledmodes)
    else
        nooffilledmodes = length(reffilledmodes)-(length(reffilledmodes)%3)
    end
    # if length(refemptymodes)%3==0
    #     noofemptymodes = length(refemptymodes)
    # else
    #     noofemptymodes = length(refemptymodes)-(length(refemptymodes)%3)
    # end    
    noofemptymodes = nooffilledmodes
    nooffrozenmodes = nooffilledmodes + noofemptymodes
    noofcouriermodes = length(secondhexagonalregionfock|>orderedmodes)-nooffrozenmodes
    @info("no of filled modes = ",nooffilledmodes)
    @info("no of courier modes = ",noofcouriermodes)
    @info("no of empty modes = ",noofemptymodes)

    rankedmodeswifeval = sort([(m,norm(v)) for (m, v) in localspectrum|>geteigenvalues],by=x->x[2])
    pickedfilledeval = rankedmodeswifeval[nooffilledmodes][2]
    pickedemptyeval = rankedmodeswifeval[end-noofemptymodes+1][2]

    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofcouriermodes, noofemptymodes])
    localcourierisometry = localstates[2]|>FockMap
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[3]|>FockMap
    courierproj = localcourierisometry*localcourierisometry'
    nooflocalcourierseeds = div(noofcouriermodes,6)
    reduandancy = length((siteAregionfock1|>orderedmodes))-nooflocalcourierseeds

    localAcorrelations1 = columns(rows(courierproj,siteAregionfock1),siteAregionfock1)
    localAstates1 = getregionstates(localcorrelations=localAcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
    localAseeds1 = localAstates1[2]

    localAcorrelations2 = columns(rows(courierproj,siteAregionfock2),siteAregionfock2)
    localAstates2 = getregionstates(localcorrelations=localAcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
    localAseeds2 = localAstates2[2]


    localAcorrelations3 = columns(rows(courierproj,siteAregionfock3),siteAregionfock3)
    localAstates3 = getregionstates(localcorrelations=localAcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
    localAseeds3 = localAstates3[2]

    allAseedsstate = localAseeds1+localAseeds2+localAseeds3

    localBcorrelations1 = columns(rows(courierproj,siteBregionfock1),siteBregionfock1)
    localBstates1 = getregionstates(localcorrelations=localBcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
    localBseeds1 = localBstates1[2]

    localBcorrelations2 = columns(rows(courierproj,siteBregionfock2),siteBregionfock2)
    localBstates2 = getregionstates(localcorrelations=localBcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
    localBseeds2 = localBstates2[2]

    localBcorrelations3 = columns(rows(courierproj,siteBregionfock3),siteBregionfock3)
    localBstates3 = getregionstates(localcorrelations=localBcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
    localBseeds3 = localBstates3[2]

    allBseedsstate = localBseeds1+localBseeds2+localBseeds3

    allseedsstates = (allAseedsstate+allBseedsstate)
    courierseeds = allseedsstates|>FockMap
    localcourierisometry

    wanniercourier = localwannierization(localcourierisometry, courierseeds)
    display((wanniercourier'*localcorrelations*wanniercourier)|>eigspec|>visualize)

    localwannierresults =  Dict(:wanniercourier => wanniercourier, :localfilledisometry => localfilledisometry,
                                :localemptyisometry => localemptyisometry, :courierseeds => courierseeds)

    wannierinfos =  Dict(secondhexagonalregionfock=>localwannierresults)

    shiftedsecondcenter1 = [0,1] ∈ firstgmeraspace
    shiftedsecondcenter2 = [-1,1] ∈ firstgmeraspace

    secondshiftedhexagonalcenterlist = [shiftedsecondcenter1, shiftedsecondcenter2]

    for hexagonalcenter in secondshiftedhexagonalcenterlist
        shiftedsecondhexagonalregion = secondhexagonalregion.+hexagonalcenter
        shiftedsecondhexagonalregionfock = quantize(shiftedsecondhexagonalregion,noofflavourpermode)
        c6recenter = recenter(c6,secondcenter+hexagonalcenter)
        c3recenter = recenter(c3,secondcenter+hexagonalcenter)
    
        shiftedsiteAregion1 = siteAregion1.+hexagonalcenter
        shiftedsiteAregionfock1 = quantize(shiftedsiteAregion1,noofflavourpermode)
        shiftedsiteAregion2 = (c3recenter)*shiftedsiteAregion1
        shiftedsiteAregionfock2 = quantize(shiftedsiteAregion2,noofflavourpermode)
        shiftedsiteAregion3 = (c3recenter)*shiftedsiteAregion2
        shiftedsiteAregionfock3 = quantize(shiftedsiteAregion3,noofflavourpermode)
    
        shiftedsiteBregion1 = (c6recenter)*shiftedsiteAregion1
        shiftedsiteBregionfock1 = quantize(shiftedsiteBregion1,noofflavourpermode)
        shiftedsiteBregion2 = (c3recenter)*shiftedsiteBregion1
        shiftedsiteBregionfock2 = quantize(shiftedsiteBregion2,noofflavourpermode)
        shiftedsiteBregion3 = (c3recenter)*shiftedsiteBregion2
        shiftedsiteBregionfock3 = quantize(shiftedsiteBregion3,noofflavourpermode)
    
        shiftedlocalcorrelations = regioncorrelations(firstgmeracorrelations,shiftedsecondhexagonalregionfock)
        shiftedlocalstates = getregionstates(localcorrelations=shiftedlocalcorrelations, grouping=[nooffilledmodes, noofcouriermodes, noofemptymodes])
        shiftedlocalcourierisometry = shiftedlocalstates[2]|>FockMap
        shiftedlocalfilledisometry = shiftedlocalstates[1]|>FockMap
        shiftedlocalemptyisometry = shiftedlocalstates[3]|>FockMap
        shiftedcourierproj = shiftedlocalcourierisometry*shiftedlocalcourierisometry'
    
        shiftedlocalAcorrelations1 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock1),shiftedsiteAregionfock1)
        shiftedlocalAstates1 = getregionstates(localcorrelations=shiftedlocalAcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalAseeds1 = shiftedlocalAstates1[2]
    
        shiftedlocalAcorrelations2 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock2),shiftedsiteAregionfock2)
        shiftedlocalAstates2 = getregionstates(localcorrelations=shiftedlocalAcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalAseeds2 = shiftedlocalAstates2[2]
    
        shiftedlocalAcorrelations3 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock3),shiftedsiteAregionfock3)
        shiftedlocalAstates3 = getregionstates(localcorrelations=shiftedlocalAcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalAseeds3 = shiftedlocalAstates3[2]
    
        shiftedallAseedsstate = shiftedlocalAseeds1+shiftedlocalAseeds2+shiftedlocalAseeds3
    
        shiftedlocalBcorrelations1 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock1),shiftedsiteBregionfock1)
        shiftedlocalBstates1 = getregionstates(localcorrelations=shiftedlocalBcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalBseeds1 = shiftedlocalBstates1[2]
    
        shiftedlocalBcorrelations2 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock2),shiftedsiteBregionfock2)
        shiftedlocalBstates2 = getregionstates(localcorrelations=shiftedlocalBcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalBseeds2 = shiftedlocalBstates2[2]
    
        shiftedlocalBcorrelations3 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock3),shiftedsiteBregionfock3)
        shiftedlocalBstates3 = getregionstates(localcorrelations=shiftedlocalBcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
        shiftedlocalBseeds3 = shiftedlocalBstates3[2]
    
        shiftedallBseedsstate = shiftedlocalBseeds1+shiftedlocalBseeds2+shiftedlocalBseeds3
    
        shiftedallseedsstates = (shiftedallAseedsstate+shiftedallBseedsstate)
        shiftedcourierseeds = shiftedallseedsstates|>FockMap
    
        shiftedwanniercourier = localwannierization(shiftedlocalcourierisometry, shiftedcourierseeds)
    
        shiftedlocalwannierresults = Dict(:wanniercourier => shiftedwanniercourier, :localfilledisometry => shiftedlocalfilledisometry,
                                        :localemptyisometry => shiftedlocalemptyisometry, :courierseeds => shiftedcourierseeds)
        wannierinfos[shiftedsecondhexagonalregionfock] =  shiftedlocalwannierresults
    end

    ref = [quantize(secondhexagonalregion.+hexagonalcenter,noofflavourpermode) for hexagonalcenter in secondshiftedhexagonalcenterlist]
    secondhexagonalregionfocklist = [secondhexagonalregionfock,ref...]

    extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in secondhexagonalregionfocklist)
    extendedwannierizedfilled =  sum(wannierinfos[regionfock][:localfilledisometry] for regionfock in secondhexagonalregionfocklist)
    extendedwannierizedempty =  sum(wannierinfos[regionfock][:localemptyisometry] for regionfock in secondhexagonalregionfocklist) 

    origin = [0, 0] ∈ firstgmeraspace
    refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

    wanniercourierisometry = globalwannierfunction(firstgmeracorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
    globalfilledisometry = globalwannierfunction(firstgmeracorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
    globalemptyisometry = globalwannierfunction(firstgmeracorrelations,extendedwannierizedempty[:,refunictcellfockempty])

    @info "Computing local courier states..."
    leftrestrict = fourier(wanniercourierisometry|>getoutspace,secondhexagonalregionfock) / (firstgmeracrystal|>vol|>sqrt)
    rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
    wanniercourierstate = leftrestrict' * wanniercourierisometry * rightrestrict

    couriercorrelations = wanniercourierisometry' * firstgmeracorrelations * wanniercourierisometry
    # display(couriercorrelations|>eigspech|>visualize)

    filledcorrelations = globalfilledisometry' * firstgmeracorrelations * globalfilledisometry
    # display(filledcorrelations|>eigspech|>visualize)

    emptycorrelations = globalemptyisometry' * firstgmeracorrelations * globalemptyisometry

    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

    return Dict(
        :couriercorrelations=>purifiedcouriercorrelations,
        :filledcorrelations=>filledcorrelations,
        :emptycorrelations=>emptycorrelations,
        :wanniercourierisometry=>wanniercourierisometry,
        :globalfilledisometry=>globalfilledisometry,
        :globalemptyisometry=>globalemptyisometry,
        :pickedfilledeval=>pickedfilledeval,
        :pickedemptyeval=>pickedemptyeval,
        :nooflocalcouriermodes=>noofcouriermodes,
        :wanniercourierstates=>wanniercourierstate|>RegionState,
        :rawcouriercorrelations=>couriercorrelations)
end
export secondgmerastepbythreshold

function thirdgmerastepbycount(secondgmeracorrelations,noofcouriermodes::Number,noofflavourpermode::Number)
    secondgmeracrystalfock = secondgmeracorrelations|>getoutspace
    secondgmeracrystal::Crystal = secondgmeracrystalfock|>getcrystal
    secondgmeraspace::RealSpace = secondgmeracrystal|>getspace

    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    thirdcenter = [1/3,1/3] ∈ secondgmeraspace
    thirdhexagonalregion = gethexagonalregion(rot = refrot,crystal=secondgmeracrystal, center=thirdcenter, metricspace=secondgmeraspace)
    thirdhexagonalregionfock = quantize(thirdhexagonalregion,noofflavourpermode)
    rgshiftedcenter1 = [2/3,-1/3] ∈ secondgmeraspace
    thirdrgshiftedhexagonalregion1 = thirdhexagonalregion.+rgshiftedcenter1
    rgshiftedcenter2 = [1/3,1/3] ∈ secondgmeraspace
    thirdrgshiftedhexagonalregion2 = thirdhexagonalregion.+rgshiftedcenter2

    c6recenter = recenter(c6,thirdcenter)
    c3recenter = recenter(c3,thirdcenter)

    siteAregion1 = intersect(intersect(thirdhexagonalregion,thirdrgshiftedhexagonalregion1),thirdrgshiftedhexagonalregion2)
    siteAregionfock1 = quantize(siteAregion1,noofflavourpermode)
    siteAregion2 = (c3recenter)*siteAregion1
    siteAregionfock2 = quantize(siteAregion2,noofflavourpermode)
    siteAregion3 = (c3recenter)*siteAregion2
    siteAregionfock3 = quantize(siteAregion3,noofflavourpermode)

    siteBregion1 = (c6recenter)*siteAregion1
    siteBregionfock1 = quantize(siteBregion1,noofflavourpermode)
    siteBregion2 = (c3recenter)*siteBregion1
    siteBregionfock2 = quantize(siteBregion2,noofflavourpermode)
    siteBregion3 = (c3recenter)*siteBregion2
    siteBregionfock3 = quantize(siteBregion3,noofflavourpermode)

    nooffrozenmodes = length(thirdhexagonalregionfock|>orderedmodes)-noofcouriermodes
    nooffilledmodes = div(nooffrozenmodes,2)
    noofemptymodes = nooffilledmodes
    @info("no of filled modes = ",nooffilledmodes)
    @info("no of courier modes = ",noofcouriermodes)
    @info("no of empty modes = ",noofemptymodes)

    localcorrelations = regioncorrelations(secondgmeracorrelations,thirdhexagonalregionfock)
    localspectrum = localcorrelations|>eigspech
    display(localspectrum|>visualize)

    if noofcouriermodes ==0 
        localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofemptymodes])
        localfilledisometry = localstates[1]|>FockMap
        localemptyisometry = localstates[2]|>FockMap

        localwannierresults =  Dict(:localfilledisometry => localfilledisometry,
                                    :localemptyisometry => localemptyisometry)

        wannierinfos =  Dict(thirdhexagonalregionfock=>localwannierresults)
        globalisometryfilled = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localfilledisometry ])
        globalisometryempty = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localemptyisometry])
        filledcorrelations = globalisometryfilled' * secondgmeracorrelations * globalisometryfilled
        display(filledcorrelations|>eigspec|>visualize)
        return globalisometryfilled,globalisometryempty
    else
        localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofcouriermodes, noofemptymodes])
        localcourierisometry = localstates[2]|>FockMap
        localfilledisometry = localstates[1]|>FockMap
        localemptyisometry = localstates[3]|>FockMap
        courierproj = localcourierisometry*localcourierisometry'
        nooflocalcourierseeds = div(noofcouriermodes,6)
        reduandancy = length((siteAregionfock1|>orderedmodes))-nooflocalcourierseeds

        localAcorrelations1 = columns(rows(courierproj,siteAregionfock1),siteAregionfock1)
        localAstates1 = getregionstates(localcorrelations=localAcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
        localAseeds1 = localAstates1[2]

        localAcorrelations2 = columns(rows(courierproj,siteAregionfock2),siteAregionfock2)
        localAstates2 = getregionstates(localcorrelations=localAcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
        localAseeds2 = localAstates2[2]


        localAcorrelations3 = columns(rows(courierproj,siteAregionfock3),siteAregionfock3)
        localAstates3 = getregionstates(localcorrelations=localAcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
        localAseeds3 = localAstates3[2]

        allAseedsstate = localAseeds1+localAseeds2+localAseeds3

        localBcorrelations1 = columns(rows(courierproj,siteBregionfock1),siteBregionfock1)
        localBstates1 = getregionstates(localcorrelations=localBcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
        localBseeds1 = localBstates1[2]

        localBcorrelations2 = columns(rows(courierproj,siteBregionfock2),siteBregionfock2)
        localBstates2 = getregionstates(localcorrelations=localBcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
        localBseeds2 = localBstates2[2]

        localBcorrelations3 = columns(rows(courierproj,siteBregionfock3),siteBregionfock3)
        localBstates3 = getregionstates(localcorrelations=localBcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
        localBseeds3 = localBstates3[2]

        allBseedsstate = localBseeds1+localBseeds2+localBseeds3

        allseedsstates = (allAseedsstate+allBseedsstate)
        courierseeds = allseedsstates|>FockMap
        localcourierisometry

        wanniercourier = localwannierization(localcourierisometry, courierseeds)
        display((wanniercourier'*localcorrelations*wanniercourier)|>eigspec|>visualize)

        localwannierresults =  Dict(:wanniercourier => wanniercourier,:localfilledisometry => localfilledisometry,
                                    :localemptyisometry => localemptyisometry, :courierseeds => courierseeds)

        wannierinfos =  Dict(thirdhexagonalregionfock=>localwannierresults)

        wanniercourierisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:wanniercourier])
        globalisometryfilled = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localfilledisometry ])
        globalisometryempty = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localemptyisometry])

        couriercorrelations = wanniercourierisometry' * secondgmeracorrelations * wanniercourierisometry

        couriercorrelationspectrum = couriercorrelations |> crystalspectrum
        purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
        purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

        display(couriercorrelations|>eigspech|>visualize)

        thirdgmeracorrelations = purifiedcouriercorrelations
        return thirdgmeracorrelations,globalisometryfilled,globalisometryempty,wanniercourierisometry
    end
end
export thirdgmerastepbycount

function thirdgmerastepbythreshold(secondgmeracorrelations,threshold::Float64,noofflavourpermode::Number)
    secondgmeracrystalfock = secondgmeracorrelations|>getoutspace
    secondgmeracrystal::Crystal = secondgmeracrystalfock|>getcrystal
    secondgmeraspace::RealSpace = secondgmeracrystal|>getspace

    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    thirdcenter = [1/3,1/3] ∈ secondgmeraspace
    thirdhexagonalregion = gethexagonalregion(rot = refrot,crystal=secondgmeracrystal, center=thirdcenter, metricspace=secondgmeraspace)
    thirdhexagonalregionfock = quantize(thirdhexagonalregion,noofflavourpermode)
    rgshiftedcenter1 = [2/3,-1/3] ∈ secondgmeraspace
    thirdrgshiftedhexagonalregion1 = thirdhexagonalregion.+rgshiftedcenter1
    rgshiftedcenter2 = [1/3,1/3] ∈ secondgmeraspace
    thirdrgshiftedhexagonalregion2 = thirdhexagonalregion.+rgshiftedcenter2

    c6recenter = recenter(c6,thirdcenter)
    c3recenter = recenter(c3,thirdcenter)

    siteAregion1 = intersect(intersect(thirdhexagonalregion,thirdrgshiftedhexagonalregion1),thirdrgshiftedhexagonalregion2)
    siteAregionfock1 = quantize(siteAregion1,noofflavourpermode)
    siteAregion2 = (c3recenter)*siteAregion1
    siteAregionfock2 = quantize(siteAregion2,noofflavourpermode)
    siteAregion3 = (c3recenter)*siteAregion2
    siteAregionfock3 = quantize(siteAregion3,noofflavourpermode)

    siteBregion1 = (c6recenter)*siteAregion1
    siteBregionfock1 = quantize(siteBregion1,noofflavourpermode)
    siteBregion2 = (c3recenter)*siteBregion1
    siteBregionfock2 = quantize(siteBregion2,noofflavourpermode)
    siteBregion3 = (c3recenter)*siteBregion2
    siteBregionfock3 = quantize(siteBregion3,noofflavourpermode)

    localcorrelations = regioncorrelations(secondgmeracorrelations,thirdhexagonalregionfock)
    localspectrum = localcorrelations|>eigspech
    display(localspectrum|>visualize)

    reffilledmodes = Subset([m for (m, v) in localspectrum|>geteigenvalues if norm(v) < threshold])
    refemptymodes = Subset([m for (m, v) in localspectrum|>geteigenvalues if norm(v) > 1-threshold])

    if length(reffilledmodes)%3==0
        nooffilledmodes = length(reffilledmodes)
    else
        nooffilledmodes = length(reffilledmodes)-(length(reffilledmodes)%3)
    end
    # if length(refemptymodes)%3==0
    #     noofemptymodes = length(refemptymodes)
    # else
    #     noofemptymodes = length(refemptymodes)-(length(refemptymodes)%3)
    # end    
    noofemptymodes = nooffilledmodes
    nooffrozenmodes = nooffilledmodes + noofemptymodes
    noofcouriermodes = length(thirdhexagonalregionfock|>orderedmodes)-nooffrozenmodes
    @info("no of filled modes = ",nooffilledmodes)
    @info("no of courier modes = ",noofcouriermodes)
    @info("no of empty modes = ",noofemptymodes)

    rankedmodeswifeval = sort([(m,norm(v)) for (m, v) in localspectrum|>geteigenvalues],by=x->x[2])
    pickedfilledeval = rankedmodeswifeval[nooffilledmodes][2]
    pickedemptyeval = rankedmodeswifeval[end-noofemptymodes+1][2]


    if noofcouriermodes ==0 
        localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofemptymodes])
        localfilledisometry = localstates[1]|>FockMap
        localemptyisometry = localstates[2]|>FockMap

        localwannierresults =  Dict(:localfilledisometry => localfilledisometry,
                                    :localemptyisometry => localemptyisometry)

        wannierinfos =  Dict(thirdhexagonalregionfock=>localwannierresults)
        globalfilledisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localfilledisometry ])
        globalemptyisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localemptyisometry])
        filledcorrelations = globalfilledisometry' * secondgmeracorrelations * globalfilledisometry
        emptycorrelations = globalemptyisometry' * secondgmeracorrelations * globalemptyisometry

        display(filledcorrelations|>eigspec|>visualize)
        return Dict(
            :filledcorrelations=>filledcorrelations,
            :emptycorrelations=>emptycorrelations,
            :globalfilledisometry=>globalfilledisometry,
            :globalemptyisometry=>globalemptyisometry,
            :pickedfilledeval=>pickedfilledeval,
            :pickedemptyeval=>pickedemptyeval)
    else
        localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofcouriermodes, noofemptymodes])
        localcourierisometry = localstates[2]|>FockMap
        localfilledisometry = localstates[1]|>FockMap
        localemptyisometry = localstates[3]|>FockMap
        courierproj = localcourierisometry*localcourierisometry'
        nooflocalcourierseeds = div(noofcouriermodes,6)
        reduandancy = length((siteAregionfock1|>orderedmodes))-nooflocalcourierseeds

        localAcorrelations1 = columns(rows(courierproj,siteAregionfock1),siteAregionfock1)
        localAstates1 = getregionstates(localcorrelations=localAcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
        localAseeds1 = localAstates1[2]

        localAcorrelations2 = columns(rows(courierproj,siteAregionfock2),siteAregionfock2)
        localAstates2 = getregionstates(localcorrelations=localAcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
        localAseeds2 = localAstates2[2]


        localAcorrelations3 = columns(rows(courierproj,siteAregionfock3),siteAregionfock3)
        localAstates3 = getregionstates(localcorrelations=localAcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
        localAseeds3 = localAstates3[2]

        allAseedsstate = localAseeds1+localAseeds2+localAseeds3

        localBcorrelations1 = columns(rows(courierproj,siteBregionfock1),siteBregionfock1)
        localBstates1 = getregionstates(localcorrelations=localBcorrelations1, grouping=[reduandancy, nooflocalcourierseeds])
        localBseeds1 = localBstates1[2]

        localBcorrelations2 = columns(rows(courierproj,siteBregionfock2),siteBregionfock2)
        localBstates2 = getregionstates(localcorrelations=localBcorrelations2, grouping=[reduandancy, nooflocalcourierseeds])
        localBseeds2 = localBstates2[2]

        localBcorrelations3 = columns(rows(courierproj,siteBregionfock3),siteBregionfock3)
        localBstates3 = getregionstates(localcorrelations=localBcorrelations3, grouping=[reduandancy, nooflocalcourierseeds])
        localBseeds3 = localBstates3[2]

        allBseedsstate = localBseeds1+localBseeds2+localBseeds3

        allseedsstates = (allAseedsstate+allBseedsstate)
        courierseeds = allseedsstates|>FockMap
        localcourierisometry

        wanniercourier = localwannierization(localcourierisometry, courierseeds)
        display((wanniercourier'*localcorrelations*wanniercourier)|>eigspec|>visualize)

        localwannierresults =  Dict(:wanniercourier => wanniercourier,:localfilledisometry => localfilledisometry,
                                    :localemptyisometry => localemptyisometry, :courierseeds => courierseeds)

        wannierinfos =  Dict(thirdhexagonalregionfock=>localwannierresults)

        wanniercourierisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:wanniercourier])
        globalfilledisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localfilledisometry ])
        globalemptyisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localemptyisometry])

        @info "Computing local courier states..."
        leftrestrict = fourier(wanniercourierisometry|>getoutspace,thirdhexagonalregionfock) / (secondgmeracrystal|>vol|>sqrt)
        rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
        wanniercourierstate = leftrestrict' * wanniercourierisometry * rightrestrict

        couriercorrelations = wanniercourierisometry' * secondgmeracorrelations * wanniercourierisometry

        filledcorrelations = globalfilledisometry' * secondgmeracorrelations * globalfilledisometry
        emptycorrelations = globalemptyisometry' * secondgmeracorrelations * globalemptyisometry

        couriercorrelationspectrum = couriercorrelations |> crystalspectrum
        purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
        purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

        display(couriercorrelations|>eigspech|>visualize)

        return Dict(
            :couriercorrelations=>purifiedcouriercorrelations,
            :filledcorrelations=>filledcorrelations,
            :emptycorrelations=>emptycorrelations,
            :wanniercourierisometry=>wanniercourierisometry,
            :globalfilledisometry=>globalfilledisometry,
            :globalemptyisometry=>globalemptyisometry,
            :pickedfilledeval=>pickedfilledeval,
            :pickedemptyeval=>pickedemptyeval,
            :nooflocalcouriermodes=>noofcouriermodes,
            :wanniercourierstates=>wanniercourierstate|>RegionState,
            :rawcouriercorrelations=>couriercorrelations)
    end
end
export thirdgmerastepbythreshold