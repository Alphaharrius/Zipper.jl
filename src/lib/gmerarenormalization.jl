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
    reciprocalhashcalibration(crystal.sizes)
    
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

# nearestneighbor = [
#         (m1, m0) => tₙ,
#         (m2, m5|> setattr(:r => [-1, 1] ∈ triangular)) => -tₙ,
#         (m1, m2) => tₙ,
#         (m3, m2) => tₙ,
#         (m4, m1|> setattr(:r => [1, 0] ∈ triangular)) => -tₙ,
#         (m3, m4) => tₙ,
#         (m5, m4) => tₙ,
#         (m0, m3|> setattr(:r => [0, -1] ∈ triangular)) => -tₙ,
#         (m5, m0) => tₙ]

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
    H = CrystalFockMap(energyspectrum) |> CrystalDenseMap

    return bonds,correlations,H
end
export generatesystem

function focktraceL1norm(fockmap,systemsize)
    fockmap = fockmap|>CrystalFockMap
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end
export focktraceL1norm

function focktraceL2norm(fockmap,systemsize)
    return (real(sqrt(tr(fockmap*fockmap'|>rep)/systemsize)))
end
export focktraceL2norm

function blocking(correlations,H,scaling)
    crystalfock = correlations|>getoutspace
    scale = Scale([scaling 0; 0 scaling], crystalfock|>getcrystal|>getspace)
    @info("Performing rgblocking...",scale)
    @info("Generating rgblocking transformation...")
    blocker = @time scale * crystalfock
    @info("Performing rgblocking on correlations...")
    blockedcorrelations = @time blocker * correlations * blocker'
    blockedH = @time blocker * H * blocker'
    return blockedcorrelations,blockedH,blocker
end
export blocking

localisometries(
    correlations::FockMap, regionfock::RegionFock;
    selectionstrategy::Function = modeselectionbythreshold(1e-3))::Dict{Symbol, FockMap} = (
    regioncorrelations(correlations, regionfock) |> selectionstrategy)
export localisometries

function sortgroupdictwifdist(dict::Dict{Mode, Number},rev::Bool)
    sortedtupledata = sort([(abs(value),key) for (key,value) in dict],rev=rev,by=first)
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

function localwannierization(localbasis::FockMap, localseeds::FockMap, svdorthothreshold::Number = 1e-1)
    U, Σ, Vt = svd(localbasis'*localseeds)
    minsvdvalue::Number = minimum(v for (_, v) in Σ)
    @info("min svdvalue", minsvdvalue)
    precarioussvdvalues::Vector = []
    svdvalues = [v for (_, v) in Σ]
    if minsvdvalue < svdorthothreshold
        push!(precarioussvdvalues, minsvdvalue)
    end
    if (precarioussvdvalues |> length) > 0
        @warn "Precarious wannier projection with minimum svdvalue of $(precarioussvdvalues |> minimum)"
        minsvd = minimum(precarioussvdvalues)
    else
        @info("all svd values greater than threshold")
        minsvd = minsvdvalue
    end
    unitary::FockMap = U * Vt
    wannierizedbasis = localbasis*unitary

    return wannierizedbasis,minsvd
end
export localwannierization


function globalwannierfunction(correlations::FockMap, localisometry::FockMap)
    wanniercrystalisos = gmeracrystalisometries(localisometry=localisometry, crystalfock=correlations|>getoutspace,addinspacemomentuminfo=true)

    wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b))+(mode |> getattr(:r)) for mode in localisometry|>getinspace |> orderedmodes)
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

function checkintersectbyeculidean(region::Subset{Offset},blockedregion::Subset{Offset})
    regiondict = Dict(euclidean(pt)=>pt for pt in region)
    blockedregiondict = Dict(euclidean(pt)=>pt for pt in blockedregion)
    euclideanregion = Subset(key for (key,val) in regiondict)
    euclideanblockedregion = Subset(key for (key,val) in blockedregiondict)
    intersection = intersect(euclideanregion,euclideanblockedregion)
    return Subset(blockedregiondict[pt] for pt in intersection)
end
export checkintersectbyeculidean

# function eecontrour()

function gmerafirststepbycount(blockedcorrelations,blockedH,noofcouriermodes::Number,noofflavourpermode::Number)
    blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal
    blockedspace::RealSpace = blockedcrystal|>getspace
    
    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    firstcenter = [0,0] ∈ blockedspace
    firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=firstcenter, metricspace=blockedspace)
    firsthexagonalregionfock = quantize(firsthexagonalregion,noofflavourpermode)
    rgshiftedcenter1 = [2/3,-1/3] ∈ blockedspace
    firstrgshiftedhexagonalregion1 = firsthexagonalregion.+rgshiftedcenter1
    rgshiftedcenter2 = [1/3,1/3] ∈ blockedspace
    firstrgshiftedhexagonalregion2 = firsthexagonalregion.+rgshiftedcenter2

    siteAregion1 = intersect(intersect(firsthexagonalregion,firstrgshiftedhexagonalregion1 ),firstrgshiftedhexagonalregion2)
    siteAregionfock1 = quantize(siteAregion1,noofflavourpermode)
    siteAregion2 = c3*siteAregion1
    siteAregionfock2 = quantize(siteAregion2,noofflavourpermode)
    siteAregion3 = (c3)*siteAregion2
    siteAregionfock3 = quantize(siteAregion3,noofflavourpermode)

    siteBregion1 = (c6)*siteAregion1
    siteBregionfock1 = quantize(siteBregion1,noofflavourpermode)
    siteBregion2 = (c3)*siteBregion1
    siteBregionfock2 = quantize(siteBregion2,noofflavourpermode)
    siteBregion3 = (c3)*siteBregion2
    siteBregionfock3 = quantize(siteBregion3,noofflavourpermode)

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

    fiosave(localcourierisometry, name="localcourierisometry")
    fiosave(localfilledisometry, name="localfilledisometry")
    fiosave(localemptyisometry, name="localemptyisometry")

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

    localwanniercourier,minsvd = localwannierization(localcourierisometry, courierseeds)
    fiosave(localwanniercourier, name="localwanniercourier")
    # display((wanniercourier'*localcorrelations*wanniercourier)|>eigspec|>visualize)

    localwannierresults =  Dict(:localwanniercourier => localwanniercourier,:localfilledisometry => localfilledisometry,
                            :localemptyisometry => localemptyisometry, :courierseeds => courierseeds)

    wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

    shiftedfirstcenter1 = [1,0] ∈ blockedspace
    shiftedfirstcenter2 = [0,1] ∈ blockedspace
    shiftedfirstcenter3 = [1,1] ∈ blockedspace
    
    firstshiftedhexagonalcenters = [shiftedfirstcenter1, shiftedfirstcenter2, shiftedfirstcenter3]
    
    for hexagonalcenter in firstshiftedhexagonalcenters
        shiftedfirsthexagonalregion = firsthexagonalregion.+hexagonalcenter
        shiftedfirsthexagonalregionfock = quantize(shiftedfirsthexagonalregion,noofflavourpermode)
        c6recenter = recenter(c6,hexagonalcenter)
        c3recenter = recenter(c3,hexagonalcenter)
        
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

        shiftedlocalwanniercourier,minsvd = localwannierization(shiftedlocalcourierisometry, shiftedcourierseeds)

        shiftedlocalwannierresults =  Dict(:localwanniercourier => shiftedlocalwanniercourier, :localfilledisometry => shiftedlocalfilledisometry,
                                        :localemptyisometry => shiftedlocalemptyisometry, :courierseeds => shiftedcourierseeds)
        wannierinfos[shiftedfirsthexagonalregionfock] =  shiftedlocalwannierresults
    end

    ref = [quantize(firsthexagonalregion.+hexagonalcenter,noofflavourpermode) for hexagonalcenter in firstshiftedhexagonalcenters]
    firsthexagonalregionfocklist = [firsthexagonalregionfock,ref...]

    extendedwannierizedcourier =  sum(wannierinfos[regionfock][:localwanniercourier] for regionfock in firsthexagonalregionfocklist)
    extendedwannierizedfilled =  sum(wannierinfos[regionfock][:localfilledisometry] for regionfock in firsthexagonalregionfocklist)
    extendedwannierizedempty =  sum(wannierinfos[regionfock][:localemptyisometry] for regionfock in firsthexagonalregionfocklist)
 
                
    origin = [0, 0] ∈ blockedspace
    refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

    wanniercourierisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
    globalfilledisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
    globalemptyisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedempty[:,refunictcellfockempty])

    fiosave(wanniercourierisometry, name="wanniercourierisometry")
    fiosave(globalfilledisometry, name="globalfilledisometry")
    fiosave(globalfilledisometry, name="globalemptyisometry")

    @info "Computing local courier states..."
    leftrestrict = fourier(wanniercourierisometry|>getoutspace, firsthexagonalregionfock) / (blockedcrystal|>vol|>sqrt)
    rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
    wanniercourierstate = leftrestrict' * wanniercourierisometry * rightrestrict


    couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
    courierH = wanniercourierisometry' * blockedH * wanniercourierisometry

    filledcorrelations = globalfilledisometry' * blockedcorrelations * globalfilledisometry
    filledH = globalfilledisometry' * blockedH * globalfilledisometry

    emptycorrelations = globalemptyisometry' * blockedcorrelations * globalemptyisometry
    emptyH = globalemptyisometry' * blockedH * globalemptyisometry

    couriercorrelationspectrum = couriercorrelations |> crystalspectrum

    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap |> CrystalDenseMap

    return Dict(
        :localcorrelations=>localcorrelations,
        :couriercorrelations=>purifiedcouriercorrelations,
        :filledcorrelations=>filledcorrelations,
        :emptycorrelations=>emptycorrelations,
        :wanniercourierisometry=>wanniercourierisometry,
        :globalfilledisometry=>globalfilledisometry,
        :globalemptyisometry=>globalemptyisometry,
        :wanniercourierstates=>wanniercourierstate|>RegionState,
        :minsvdforcourier=>minsvd,
        :rawcouriercorrelations=>couriercorrelations,
        :courierH=>courierH,
        :filledH=>filledH,
        :emptyH=>emptyH)
end
export gmerafirststepbycount

function gmerasecondstepbycount(firstgmeracorrelations,firstgmeraH,noofcouriermodes::Number,noofflavourpermode::Number)
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

    localwanniercourier,minsvd = localwannierization(localcourierisometry, courierseeds)
    # display((localwanniercourier'*localcorrelations*localwanniercourier)|>eigspec|>visualize)

    localwannierresults =  Dict(:localwanniercourier => localwanniercourier, :localfilledisometry => localfilledisometry,
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
    
        shiftedlocalwanniercourier,minsvd = localwannierization(shiftedlocalcourierisometry, shiftedcourierseeds)
    
        shiftedlocalwannierresults = Dict(:localwanniercourier => shiftedlocalwanniercourier, :localfilledisometry => shiftedlocalfilledisometry,
                                        :localemptyisometry => shiftedlocalemptyisometry, :courierseeds => shiftedcourierseeds)
        wannierinfos[shiftedsecondhexagonalregionfock] =  shiftedlocalwannierresults
    end

    ref = [quantize(secondhexagonalregion.+hexagonalcenter,noofflavourpermode) for hexagonalcenter in secondshiftedhexagonalcenterlist]
    secondhexagonalregionfocklist = [secondhexagonalregionfock,ref...]

    extendedwannierizedcourier =  sum(wannierinfos[regionfock][:localwanniercourier] for regionfock in secondhexagonalregionfocklist)
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
    courierH = wanniercourierisometry' * firstgmeraH * wanniercourierisometry
    # display(couriercorrelations|>eigspech|>visualize)

    filledcorrelations = globalfilledisometry' * firstgmeracorrelations * globalfilledisometry
    filledH = globalfilledisometry' * firstgmeraH * globalfilledisometry
    # display(filledcorrelations|>eigspech|>visualize)

    emptycorrelations = globalemptyisometry' * firstgmeracorrelations * globalemptyisometry
    emptyH = globalemptyisometry' * firstgmeraH * globalemptyisometry

    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap |> CrystalDenseMap

    return Dict(
        :localcorrelations=>localcorrelations,
        :couriercorrelations=>purifiedcouriercorrelations,
        :filledcorrelations=>filledcorrelations,
        :emptycorrelations=>emptycorrelations,
        :wanniercourierisometry=>wanniercourierisometry,
        :globalfilledisometry=>globalfilledisometry,
        :globalemptyisometry=>globalemptyisometry,
        :wanniercourierstates=>wanniercourierstate|>RegionState,
        :minsvdforcourier=>minsvd,
        :rawcouriercorrelations=>couriercorrelations,
        :courierH=>courierH,
        :filledH=>filledH,
        :emptyH=>emptyH)
end
export gmerasecondstepbycount

function gmerathirdstepbycount(secondgmeracorrelations,secondgmeraH,noofcouriermodes::Number,noofflavourpermode::Number)
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
        globalfilledisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localfilledisometry ])
        globalemptyisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localemptyisometry])
        filledcorrelations = globalfilledisometry' * secondgmeracorrelations * globalfilledisometry
        emptycorrelations = globalemptyisometry' * secondgmeracorrelations * globalemptyisometry
        filledH = globalfilledisometry' * secondgmeraH * globalfilledisometry
        emptyH = globalemptyisometry' * secondgmeraH * globalemptyisometry
        # display(filledcorrelations|>eigspec|>visualize)
        return return Dict(
            :localcorrelations=>localcorrelations,
            :filledcorrelations=>filledcorrelations,
            :emptycorrelations=>emptycorrelations,
            :globalfilledisometry=>globalfilledisometry,
            :globalemptyisometry=>globalemptyisometry,
            filledH=>filledH,
            emptyH=>emptyH)
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

        localwanniercourier,minsvd = localwannierization(localcourierisometry, courierseeds)
        # display((wanniercourier'*localcorrelations*wanniercourier)|>eigspec|>visualize)

        localwannierresults =  Dict(:localwanniercourier => localwanniercourier,:localfilledisometry => localfilledisometry,
                                    :localemptyisometry => localemptyisometry, :courierseeds => courierseeds)

        wannierinfos =  Dict(thirdhexagonalregionfock=>localwannierresults)

        wanniercourierisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localwanniercourier])
        globalfilledisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localfilledisometry ])
        globalemptyisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localemptyisometry])

        @info "Computing local courier states..."
        leftrestrict = fourier(wanniercourierisometry|>getoutspace,thirdhexagonalregionfock) / (secondgmeracrystal|>vol|>sqrt)
        rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
        wanniercourierstate = leftrestrict' * wanniercourierisometry * rightrestrict

        couriercorrelations = wanniercourierisometry' * secondgmeracorrelations * wanniercourierisometry
        filledcorrelations = globalfilledisometry' * secondgmeracorrelations * globalfilledisometry
        emptycorrelations = globalemptyisometry' * secondgmeracorrelations * globalemptyisometry

        courierH = wanniercourierisometry' * secondgmeraH * wanniercourierisometry
        filledH = globalfilledisometry' * secondgmeraH * globalfilledisometry
        emptyH = globalemptyisometry' * secondgmeraH * globalemptyisometry

        couriercorrelationspectrum = couriercorrelations |> crystalspectrum
        purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
        purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap |> CrystalDenseMap

        # display(couriercorrelations|>eigspech|>visualize)

        return return Dict(
            :localcorrelations=>localcorrelations,
            :couriercorrelations=>purifiedcouriercorrelations,
            :filledcorrelations=>filledcorrelations,
            :emptycorrelations=>emptycorrelations,
            :wanniercourierisometry=>wanniercourierisometry,
            :globalfilledisometry=>globalfilledisometry,
            :globalemptyisometry=>globalemptyisometry,
            :wanniercourierstates=>wanniercourierstate|>RegionState,
            :minsvdforcourier=>minsvd,
            :rawcouriercorrelations=>couriercorrelations,
            :courierH=>courierH,
            :filledH=>filledH,
            :emptyH=>emptyH)
    end
end
export gmerathirdstepbycount

function gmerafinalstep(correlations,H,noofflavourpermode::Number)
    crystalfock = correlations|>getoutspace
    crystal::Crystal = crystalfock|>getcrystal
    space::RealSpace = crystal|>getspace
    if (crystal|>size)[1]==1
        firstcenter = [0,0] ∈ space
        refrot = inv([2/3 -1/3; -1/3 2/3]')
        hexagonalregion = gethexagonalregion(rot = refrot,crystal=crystal, center=firstcenter, metricspace=space)
        hexagonalregionfock = quantize(hexagonalregion,noofflavourpermode)

        nooffrozenmodes = length(hexagonalregionfock|>orderedmodes)
        nooffilledmodes = div(nooffrozenmodes,2)
        noofemptymodes = nooffilledmodes
        @info("no of filled modes = ",nooffilledmodes)
        @info("no of empty modes = ",noofemptymodes)

        localcorrelations = regioncorrelations(correlations,hexagonalregionfock)
        localspectrum = localcorrelations|>eigspec
        display(localspectrum|>visualize)
        localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofemptymodes])
        localfilledisometry = localstates[1]|>FockMap
        localemptyisometry = localstates[2]|>FockMap

        localwannierresults =  Dict(:localfilledisometry => localfilledisometry,
                                    :localemptyisometry => localemptyisometry)

        wannierinfos =  Dict(hexagonalregionfock=>localwannierresults)
        globalfilledisometry = globalwannierfunction(correlations,wannierinfos[hexagonalregionfock][:localfilledisometry ])
        globalemptyisometry = globalwannierfunction(correlations,wannierinfos[hexagonalregionfock][:localemptyisometry])

        filledcorrelations = globalfilledisometry' * correlations * globalfilledisometry
        emptycorrelations = globalemptyisometry' * correlations * globalemptyisometry
        filledH = globalfilledisometry' * H * globalfilledisometry
        emptyH = globalemptyisometry' * H * globalemptyisometry

        return Dict(
        :localcorrelations=>localcorrelations,
        :filledcorrelations=>filledcorrelations,
        :emptycorrelations=>emptycorrelations,
        :globalfilledisometry=>globalfilledisometry,
        :globalemptyisometry=>globalemptyisometry,
        :filledH=>filledH,
        :emptyH=>emptyH)
    else
        @warn("haven't reached the last step of RG")
        return
    end
end
export gmerafinalstep

function firstgmerastep(correlations,H,scaling)
    blockedcorrelations,blockedH,blocker = blocking(correlations,H,scaling)
    fiosave(blockedcorrelations, name="blockedcorrelations")
    fiosave(blockedH, name="blockedH")
    fiosave(blocker, name="blocker")
    noofmodesinlocalreg = (6*(scaling^2))|>Int
    noofdistillablemodes = ((1/4)*noofmodesinlocalreg)|>Int
    noofcouriermodesinfirststep = noofmodesinlocalreg - noofdistillablemodes
    println(noofcouriermodesinfirststep)
    noofcouriermodesinsecondstep = noofcouriermodesinfirststep - noofdistillablemodes
    noofcouriermodesinthirdstep = noofcouriermodesinsecondstep - noofdistillablemodes

    gmera1firststepdata = gmerafirststepbycount(blockedcorrelations,blockedH,noofcouriermodesinfirststep,1)
    gmera1secondstepdata = gmerasecondstepbycount(gmera1firststepdata[:couriercorrelations],gmera1firststepdata[:courierH],noofcouriermodesinsecondstep,(noofcouriermodesinfirststep/6)|>Int)
    gmera1thirdstepdata = gmerathirdstepbycount(gmera1secondstepdata[:couriercorrelations],gmera1secondstepdata[:courierH],noofcouriermodesinthirdstep,(noofcouriermodesinsecondstep/6)|>Int)

    noofflavourpermodeforlaterrg = (noofcouriermodesinthirdstep/6)|>Int

    gmera1firstapproximation = gmera1firststepdata[:globalemptyisometry]*gmera1firststepdata[:globalemptyisometry]'
    gmera1secondapproximation = (gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:globalemptyisometry])*(gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:globalemptyisometry])'
    gmera1thirdapproximation = (gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:wanniercourierisometry]*gmera1thirdstepdata[:globalemptyisometry])*(gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:wanniercourierisometry]*gmera1thirdstepdata[:globalemptyisometry])'

    couriercomposemapgmera1 = (gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:wanniercourierisometry]*gmera1thirdstepdata[:wanniercourierisometry])

    gmera1approximatecorrelation = gmera1firstapproximation + gmera1secondapproximation + gmera1thirdapproximation
    fiosave(gmera1approximatecorrelation, name="gmera1approximatecorrelation")
    fiosave(gmera1firststepdata[:localcorrelations], name="gmera1firstlocalcorrelations")
    fiosave(gmera1firststepdata[:couriercorrelations], name="gmera1firstcouriercorrelations")
    fiosave(gmera1firststepdata[:filledcorrelations], name="gmera1firstfilledcorrelations")
    fiosave(gmera1firststepdata[:emptycorrelations], name="gmera1firstemptycorrelations")
    fiosave(gmera1firststepdata[:courierH], name="gmera1firstcourierH")
    fiosave(gmera1firststepdata[:filledH], name="gmera1firstfilledH")
    fiosave(gmera1firststepdata[:emptyH], name="gmera1firstemptyH")
    fiosave(gmera1firststepdata[:rawcouriercorrelations], name="gmera1firstrawcouriercorrelations")
    fiosave(gmera1firststepdata[:wanniercourierisometry], name="gmera1firstwanniercourierisometry")
    fiosave(gmera1firststepdata[:globalemptyisometry], name="gmera1firstglobalemptyisometry")
    fiosave(gmera1firststepdata[:globalfilledisometry], name="gmera1firstglobalfilledisometry")
    fiosave(gmera1firststepdata[:wanniercourierstates], name="gmera1firstwanniercourierstates")
    fiosave(gmera1firststepdata[:minsvdforcourier], name="gmera1firstminsvdforcourier")

    fiosave(gmera1secondstepdata[:localcorrelations], name="gmera1secondlocalcorrelations")
    fiosave(gmera1secondstepdata[:couriercorrelations], name="gmera1secondcouriercorrelations")
    fiosave(gmera1secondstepdata[:filledcorrelations], name="gmera1secondfilledcorrelations")
    fiosave(gmera1secondstepdata[:emptycorrelations], name="gmera1secondemptycorrelations")
    fiosave(gmera1secondstepdata[:courierH], name="gmera1secondcourierH")
    fiosave(gmera1secondstepdata[:filledH], name="gmera1secondfilledH")
    fiosave(gmera1secondstepdata[:emptyH], name="gmera1secondemptyH")
    fiosave(gmera1secondstepdata[:rawcouriercorrelations], name="gmera1secondrawcouriercorrelations")
    fiosave(gmera1secondstepdata[:wanniercourierisometry], name="gmera1secondwanniercourierisometry")
    fiosave(gmera1secondstepdata[:globalemptyisometry], name="gmera1secondglobalemptyisometry")
    fiosave(gmera1secondstepdata[:globalfilledisometry], name="gmera1secondglobalfilledisometry")
    fiosave(gmera1secondstepdata[:wanniercourierstates], name="gmera1secondwanniercourierstates")
    fiosave(gmera1secondstepdata[:minsvdforcourier], name="gmera1secondminsvdforcourier")

    fiosave(gmera1thirdstepdata[:localcorrelations], name="gmera1thirdlocalcorrelations")
    fiosave(gmera1thirdstepdata[:couriercorrelations], name="gmera1thirdcouriercorrelations")
    fiosave(gmera1thirdstepdata[:filledcorrelations], name="gmera1thirdfilledcorrelations")
    fiosave(gmera1thirdstepdata[:emptycorrelations], name="gmera1thirdemptycorrelations")
    fiosave(gmera1thirdstepdata[:courierH], name="gmera1thirdcourierH")
    fiosave(gmera1thirdstepdata[:filledH], name="gmera1thirdfilledH")
    fiosave(gmera1thirdstepdata[:emptyH], name="gmera1thirdemptyH")
    fiosave(gmera1thirdstepdata[:rawcouriercorrelations], name="gmera1thirdrawcouriercorrelations")
    fiosave(gmera1thirdstepdata[:wanniercourierisometry], name="gmera1thirdwanniercourierisometry")
    fiosave(gmera1thirdstepdata[:globalemptyisometry], name="gmera1thirdglobalemptyisometry")
    fiosave(gmera1thirdstepdata[:globalfilledisometry], name="gmera1thirdglobalfilledisometry")
    fiosave(gmera1thirdstepdata[:wanniercourierstates], name="gmera1thirdwanniercourierstates")
    fiosave(gmera1thirdstepdata[:minsvdforcourier], name="gmera1thirdminsvdforcourier")

    rg1correlations = gmera1thirdstepdata[:couriercorrelations]

    rg1H = gmera1thirdstepdata[:courierH]   
    return  rg1H,rg1correlations,couriercomposemapgmera1,gmera1approximatecorrelation,blockedcorrelations,noofflavourpermodeforlaterrg
end
export firstgmerastep

function intermediategmerastep(rgcorrelations,rgH,couriercomposemap,gmerprevsumapproximatecorrelation,rgstep,noofflavourpermode)
    rgblockedcorrelations,rgblockedH,rgblocker = blocking(rgcorrelations,rgH,2)
    prevrgstep = rgstep-1
    fiosave(rgblockedcorrelations, name="rg$prevrgstep"*"blockedcorrelations")
    fiosave(rgblockedH, name="rg$prevrgstep"*"blockedH")
    fiosave(rgblocker, name="rg$prevrgstep"*"blocker")
    couriercomposemap = couriercomposemap*rgblocker'

    gmerafirststepdata = gmerafirststepbycount(rgblockedcorrelations,rgblockedH,18*noofflavourpermode,noofflavourpermode)
    gmerasecondstepdata = gmerasecondstepbycount(gmerafirststepdata[:couriercorrelations],gmerafirststepdata[:courierH],12*noofflavourpermode,3*noofflavourpermode)
    gmerathirdstepdata = gmerathirdstepbycount(gmerasecondstepdata[:couriercorrelations],gmerasecondstepdata[:courierH],6*noofflavourpermode,2*noofflavourpermode)

    
    gmerafirstapproximation = (couriercomposemap*gmerafirststepdata[:globalemptyisometry])*(couriercomposemap*gmerafirststepdata[:globalemptyisometry])'
    gmerasecondapproximation = (couriercomposemap*gmerafirststepdata[:wanniercourierisometry]*gmerasecondstepdata[:globalemptyisometry])*(couriercomposemap*gmerafirststepdata[:wanniercourierisometry]*gmerasecondstepdata[:globalemptyisometry])'
    gmerathirdapproximation = (couriercomposemap*gmerafirststepdata[:wanniercourierisometry]*gmerasecondstepdata[:wanniercourierisometry]*gmerathirdstepdata[:globalemptyisometry])*(couriercomposemap*gmerafirststepdata[:wanniercourierisometry]*gmerasecondstepdata[:wanniercourierisometry]*gmerathirdstepdata[:globalemptyisometry])'
    
    couriercomposemapgmera = couriercomposemap*(gmerafirststepdata[:wanniercourierisometry]*gmerasecondstepdata[:wanniercourierisometry]*gmerathirdstepdata[:wanniercourierisometry])
    rgsteps = prod([string(i) for i in rgstep])
    fiosave(couriercomposemapgmera, name="gmera$rgsteps"*"couriercomposemapgmera")

    gmeraapproximatecorrelation = gmerafirstapproximation + gmerasecondapproximation  + gmerathirdapproximation
    fiosave(gmeraapproximatecorrelation, name="gmera$rgstep"*"approximatecorrelation")

    fiosave(gmerafirststepdata[:localcorrelations], name="gmera$rgstep"*"firstlocalcorrelations")
    fiosave(gmerafirststepdata[:couriercorrelations], name="gmera$rgstep"*"firstcouriercorrelations")
    fiosave(gmerafirststepdata[:filledcorrelations], name="gmera$rgstep"*"firstfilledcorrelations")
    fiosave(gmerafirststepdata[:emptycorrelations], name="gmera$rgstep"*"firstemptycorrelations")
    fiosave(gmerafirststepdata[:courierH], name="gmera$rgstep"*"firstcourierH")
    fiosave(gmerafirststepdata[:filledH], name="gmera$rgstep"*"firstfilledH")
    fiosave(gmerafirststepdata[:emptyH], name="gmera$rgstep"*"firstemptyH")
    fiosave(gmerafirststepdata[:rawcouriercorrelations], name="gmera$rgstep"*"firstrawcouriercorrelations")
    fiosave(gmerafirststepdata[:wanniercourierisometry], name="gmera$rgstep"*"firstwanniercourierisometry")
    fiosave(gmerafirststepdata[:globalemptyisometry], name="gmera$rgstep"*"firstglobalemptyisometry")
    fiosave(gmerafirststepdata[:globalfilledisometry], name="gmera$rgstep"*"firstglobalfilledisometry")
    fiosave(gmerafirststepdata[:wanniercourierstates], name="gmera$rgstep"*"firstwanniercourierstates")
    fiosave(gmerafirststepdata[:minsvdforcourier], name="gmera$rgstep"*"firstminsvdforcourier")

    fiosave(gmerasecondstepdata[:localcorrelations], name="gmera$rgstep"*"secondlocalcorrelations")
    fiosave(gmerasecondstepdata[:couriercorrelations], name="gmera$rgstep"*"secondcouriercorrelations")
    fiosave(gmerasecondstepdata[:filledcorrelations], name="gmera$rgstep"*"secondfilledcorrelations")
    fiosave(gmerasecondstepdata[:emptycorrelations], name="gmera$rgstep"*"secondemptycorrelations")
    fiosave(gmerasecondstepdata[:courierH], name="gmera$rgstep"*"secondcourierH")
    fiosave(gmerasecondstepdata[:filledH], name="gmera$rgstep"*"secondfilledH")
    fiosave(gmerasecondstepdata[:emptyH], name="gmera$rgstep"*"secondemptyH")
    fiosave(gmerasecondstepdata[:rawcouriercorrelations], name="gmera$rgstep"*"secondrawcouriercorrelations")
    fiosave(gmerasecondstepdata[:wanniercourierisometry], name="gmera$rgstep"*"secondwanniercourierisometry")
    fiosave(gmerasecondstepdata[:globalemptyisometry], name="gmera$rgstep"*"secondglobalemptyisometry")
    fiosave(gmerasecondstepdata[:globalfilledisometry], name="gmera$rgstep"*"secondglobalfilledisometry")
    fiosave(gmerasecondstepdata[:wanniercourierstates], name="gmera$rgstep"*"secondwanniercourierstates")
    fiosave(gmerasecondstepdata[:minsvdforcourier], name="gmera$rgstep"*"secondminsvdforcourier")

    fiosave(gmerathirdstepdata[:localcorrelations], name="gmera$rgstep"*"thirdlocalcorrelations")
    fiosave(gmerathirdstepdata[:couriercorrelations], name="gmera$rgstep"*"thirdcouriercorrelations")
    fiosave(gmerathirdstepdata[:filledcorrelations], name="gmera$rgstep"*"thirdfilledcorrelations")
    fiosave(gmerathirdstepdata[:emptycorrelations], name="gmera$rgstep"*"thirdemptycorrelations")
    fiosave(gmerathirdstepdata[:courierH], name="gmera$rgstep"*"thirdcourierH")
    fiosave(gmerathirdstepdata[:filledH], name="gmera$rgstep"*"thirdfilledH")
    fiosave(gmerathirdstepdata[:emptyH], name="gmera$rgstep"*"thirdemptyH")
    fiosave(gmerathirdstepdata[:rawcouriercorrelations], name="gmera$rgstep"*"thirdrawcouriercorrelations")
    fiosave(gmerathirdstepdata[:wanniercourierisometry], name="gmera$rgstep"*"thirdwanniercourierisometry")
    fiosave(gmerathirdstepdata[:globalemptyisometry], name="gmera$rgstep"*"thirdglobalemptyisometry")
    fiosave(gmerathirdstepdata[:globalfilledisometry], name="gmera$rgstep"*"thirdglobalfilledisometry")
    fiosave(gmerathirdstepdata[:wanniercourierstates], name="gmera$rgstep"*"thirdwanniercourierstates")
    fiosave(gmerathirdstepdata[:minsvdforcourier], name="gmera$rgstep"*"thirdminsvdforcourier")

    rgcorrelations = gmerathirdstepdata[:couriercorrelations]
    rgH = gmerathirdstepdata[:courierH]   
    gmerasumapproximatecorrelationsofar = gmerprevsumapproximatecorrelation+gmeraapproximatecorrelation
    fiosave(gmerasumapproximatecorrelationsofar, name="gmera$rgsteps"*"approximatecorrelationsofar")
    return  rgH,rgcorrelations,couriercomposemapgmera,gmerasumapproximatecorrelationsofar
end
export intermediategmerastep

function finalgmerastep(rgcorrelations,rgH,couriercomposemap,gmerprevsumapproximatecorrelation,blockedcorrelations,rgstep,noofflavourpermode,systemsize)
    rgblockedcorrelations,rgblockedH,rgblocker = blocking(rgcorrelations,rgH,2)
    prevrgstep = rgstep-1
    fiosave(rgblockedcorrelations, name="rg$prevrgstep"*"blockedcorrelations")
    fiosave(rgblockedH, name="rg$prevrgstep"*"blockedH")
    fiosave(rgblocker, name="rg$prevrgstep"*"blocker")
    couriercomposemap = couriercomposemap*rgblocker'

    gmerafinalstepdata = gmerafinalstep(rgblockedcorrelations,rgblockedH,noofflavourpermode)
    fiosave(gmerafinalstepdata[:localcorrelations], name="gmera$rgstep"*"finallocalcorrelations")
    fiosave(gmerafinalstepdata[:filledcorrelations], name="gmera$rgstep"*"finalfilledcorrelations")
    fiosave(gmerafinalstepdata[:emptycorrelations], name="gmera$rgstep"*"finalemptycorrelations")
    fiosave(gmerafinalstepdata[:filledH], name="gmera$rgstep"*"finalfilledH")
    fiosave(gmerafinalstepdata[:emptyH], name="gmera$rgstep"*"finalemptyH")
    fiosave(gmerafinalstepdata[:globalemptyisometry], name="gmera$rgstep"*"finalglobalemptyisometry")
    fiosave(gmerafinalstepdata[:globalfilledisometry], name="gmera$rgstep"*"finalglobalfilledisometry")

    gmerafinalapproximation = (couriercomposemap*gmerafinalstepdata[:globalemptyisometry])*(couriercomposemap*gmerafinalstepdata[:globalemptyisometry])'
    fiosave(gmerafinalapproximation, name="gmera$rgstep"*"finalapproximation")
    gmeraallapproximatecorrelation = gmerprevsumapproximatecorrelation+gmerafinalapproximation
    fiosave(gmeraallapproximatecorrelation, name="gmeraallapproximatecorrelation")

    gmeraalldiff = blockedcorrelations - gmeraallapproximatecorrelation
    fiosave(gmeraalldiff, name="gmeraalldiff")

    @info("calculating L1 norm of the difference between the blocked correlations and the approximate correlation")
    gmeraalldiffL1norm = focktraceL1norm(gmeraalldiff,systemsize^2*6)
    fiosave(gmeraalldiffL1norm, name="gmeraalldiffL1norm")

    # @info("calculating L2 norm of the difference between the blocked correlations and the approximate correlation")
    # gmeraalldiffL2norm = focktraceL2norm(gmeraalldiff,systemsize^2*6)
    # fiosave(gmeraalldiffL2norm, name="gmeraalldiffL2norm")
end
export finalgmerastep