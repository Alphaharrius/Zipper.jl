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

    return correlations,H
end
export generatesystem

function ee(eval::Number)
    if 0<round(eval,digits=10)<1
        return (-eval*log(eval)-(1-eval)*log((1-eval)))
    else 
        return 0 
    end
end
export ee

function entanglementcontour(localcorrelations::SparseFockMap)
    localspectrum = localcorrelations|>eigspec
    emodewifevals = localspectrum|>geteigenvalues
    evectors = localspectrum|>geteigenvectors
    rsfock = localcorrelations|>getoutspace
    result::Dict{Mode, Number} = Dict()
    for rsmode in rsfock
        save = []
        for (emode,eval) in emodewifevals
            append!(save,norm(evectors[rsmode,emode])^2*ee(norm(eval)))
        end
        result[rsmode] = sum(save)
    end
    return result
end
export entanglementcontour

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

function crystalisometry(spectrum::CrystalSpectrum)::FockMap
    blocks::Dict = paralleltasks(
        name="crystalprojector",
        tasks=(()->((k, k)=>u) for (k, u) in spectrum|>geteigenvectors),
        count=spectrum|>getcrystal|>vol)|>parallel|>Dict
    return crystalfockmap(spectrum|>getcrystal, spectrum|>getcrystal, blocks)
end
export crystalisometry

function sortgroupdictwifvalue(dict::Dict{Mode, Number},rev::Bool,roundingdigit::Int64=8)
    sortedtupledata = sort([(value,key) for (key,value) in dict],rev=rev,by=first)
    refvalue = round(sortedtupledata[1][1],digits=roundingdigit)
    result = []
    subresult = Subset(sortedtupledata[1][2])
    for pair in sortedtupledata
        if round(pair[1],digits=roundingdigit)==refvalue
            subresult = subresult+Subset(pair[2])
        else
            append!(result,tuple([refvalue,subresult]))
            refvalue = round(pair[1],digits=roundingdigit)
            subresult = Subset(pair[2])
        end
    end
    append!(result,tuple([refvalue,subresult]))
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

function offsetofmodes(modes::Subset)
    return Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in modes)
end
export offsetofmodes

function basisoffsetofmodes(modes::Subset)
    return Subset((mode|>getattr(:b)) for mode in modes)
end
export basisoffsetofmodes

# function gmerafirststep(blockedcorrelations,blockedH)
#     blockedcrystalfock = blockedcorrelations|>getoutspace
#     blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal
#     blockedspace::RealSpace = blockedcrystal|>getspace
    
#     refrot = inv([2/3 -1/3; -1/3 2/3]')
#     c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
#     c3 = c6^2

#     firstcenter = [0,0] ∈ blockedspace
#     firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=firstcenter, metricspace=blockedspace)
#     firsthexagonalregionfock = quantize(firsthexagonalregion,1)
#     firsthexagonalmodes = Subset(mode for mode in firsthexagonalregionfock)

#     # @info "Computing local correlations..."
#     localrestrict = fourier(blockedcrystalfock, firsthexagonalregionfock) / (blockedcrystal|>vol|>sqrt)
#     localcorrelations = localrestrict'*blockedcorrelations*localrestrict
#     potentialfilledseeds = Subset([mode for mode in firsthexagonalregionfock if (localcorrelations[mode,mode]|>rep)[1,1]|>abs<0.5])
#     potentialemptyseeds = Subset([mode for mode in firsthexagonalregionfock if (localcorrelations[mode,mode]|>rep)[1,1]|>abs>0.5])

#     localspectrum = localcorrelations|>eigspec
#     entanglementcontourinfo = entanglementcontour(localcorrelations)
#     groupedandsortedmodeswifecontour = sortgroupdictwifvalue(entanglementcontourinfo,false)
#     nooffrozenmodes = div(firsthexagonalregionfock|>length,4)*1
#     frozenseeds = groupedandsortedmodeswifecontour[1][2]
#     for (entanglementval,modes) in groupedandsortedmodeswifecontour[2:end]
#         if (frozenseeds|>length)<nooffrozenmodes
#             @info "include more seeds"
#             frozenseeds = frozenseeds+modes
#         elseif (frozenseeds|>length)==nooffrozenmodes
#             @info "just enough frozenseeds"
#             break
#         else 
#             @error "too many frozenseeds, sth wrong"
#         end
#     end

#     filledseeds = intersect(frozenseeds,potentialfilledseeds)
#     emptyseeds = intersect(frozenseeds,potentialemptyseeds)
#     courierseeds = firsthexagonalmodes - frozenseeds
#     # courierseeds|>offsetofmodes|>visualize
#     # frozenseeds|>offsetofmodes|>visualize
#     filledseedsfock = quantize(filledseeds|>offsetofmodes,1)
#     emptyseedsfock = quantize(emptyseeds|>offsetofmodes,1)
#     courierseedsfock = quantize(courierseeds|>offsetofmodes,1)

#     noooffilledmodes = length(filledseeds)
#     nooofemptymodes = length(emptyseeds)
#     noofcouriermodes = length(courierseeds)
    
#     localstates = getregionstates(localcorrelations=localcorrelations, grouping=[noooffilledmodes,noofcouriermodes,nooofemptymodes])
#     localcourieriso = localstates[2]|>FockMap
#     localfillediso = localstates[1]|>FockMap
#     localemptyiso = localstates[3]|>FockMap

#     iden = idmap(localcorrelations|>getoutspace, localcorrelations|>getinspace)
#     localfilledseedsiso = iden[:,filledseedsfock]
#     localemptyseedsiso = iden[:,emptyseedsfock]
#     localcourierseedsiso = iden[:,courierseedsfock]


#     localwannierfillediso,minsvdempty = localwannierization(localfillediso, localfilledseedsiso)
#     localwannieremptyiso,minsvdfilled = localwannierization(localemptyiso, localemptyseedsiso)
#     localwanniercourieriso,minsvdcourier = localwannierization(localcourieriso, localcourierseedsiso)
#     localdisentangler = localwanniercourieriso+localwannierfillediso+localwannieremptyiso
#     redefinelocaldisentangler = FockMap(localdisentangler|>getoutspace,localdisentangler|>getinspace|>RegionFock,localdisentangler|>rep)

#     newcrystalfock = getcrystalfock(redefinelocaldisentangler|>getinspace|>unitcellfock, blockedcrystal)
#     newlocalrestrict = fourier(newcrystalfock, redefinelocaldisentangler|>getinspace) * (blockedcrystal|>vol|>sqrt)

#     globaldisentangler = broadcast(*,(localrestrict*redefinelocaldisentangler), newlocalrestrict')
#     firstgmeracorrelations = globaldisentangler'*blockedcorrelations*globaldisentangler
#     firstgmeraH = globaldisentangler'*blockedH*globaldisentangler


#     return Dict(
#         :localcorrelations=>localcorrelations,
#         :globaldisentangler=>globaldisentangler,
#         :firstgmeracorrelations=>firstgmeracorrelations,
#         :firstgmeraH=>firstgmeraH,
#         :minsvdcourier=>minsvdcourier,
#         :minsvdempty=>minsvdempty,
#         :minsvdfilled=>minsvdfilled)
# end
# export gmerafirststep

# function gmerasecondstep(firstgmeracorrelations,firstgmeraH)
#     firstgmeracrystalfock = firstgmeracorrelations|>getoutspace
#     firstgmeracrystal::Crystal = firstgmeracrystalfock|>getcrystal
#     firstgmeraspace::RealSpace = firstgmeracrystal|>getspace

#     refrot = inv([2/3 -1/3; -1/3 2/3]')
#     c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
#     c3 = c6^2

#     secondcenter = [2/3,-1/3] ∈ firstgmeraspace
#     secondhexagonalregion = gethexagonalregion(rot=refrot,crystal=firstgmeracrystal, center=secondcenter, metricspace=firstgmeraspace)
#     secondhexagonalregionfock = quantize(secondhexagonalregion,1)
#     secondhexagonalmodes = Subset(mode for mode in secondhexagonalregionfock)

#     c6recenter = recenter(c6,secondcenter)
#     c3recenter = recenter(c3,secondcenter)

#     # @info "Computing local correlations..."
#     localrestrict = fourier(firstgmeracrystalfock, secondhexagonalregionfock) / (firstgmeracrystal|>vol|>sqrt)
#     localcorrelations = localrestrict'*firstgmeracorrelations*localrestrict
#     potentialfilledseeds = Subset([mode for mode in secondhexagonalregionfock if (localcorrelations[mode,mode]|>rep)[1,1]|>abs<0.5])
#     potentialemptyseeds = Subset([mode for mode in secondhexagonalregionfock if (localcorrelations[mode,mode]|>rep)[1,1]|>abs>0.5])

#     localspectrum = localcorrelations|>eigspech
#     localspectrum|>visualize|>display
#     entanglementcontourinfo = entanglementcontour(localcorrelations)
#     groupedandsortedmodeswifecontour = sortgroupdictwifvalue(entanglementcontourinfo,false)
#     nooffrozenmodes = div(secondhexagonalregionfock|>length,4)*2
#     frozenseeds = groupedandsortedmodeswifecontour[1][2]
#     for (entanglementval,modes) in groupedandsortedmodeswifecontour[2:end]
#         if (frozenseeds|>length)<nooffrozenmodes
#             @info "include more seeds"
#             frozenseeds = frozenseeds+modes
#         elseif (frozenseeds|>length)==nooffrozenmodes
#             @info "just enough frozenseeds"
#             break
#         else 
#             @error "too many frozenseeds, sth wrong"
#         end
#     end
#     filledseeds = intersect(frozenseeds,potentialfilledseeds)
#     emptyseeds = intersect(frozenseeds,potentialemptyseeds)
#     courierseeds = secondhexagonalmodes - frozenseeds
#     filledseedsfock = quantize(filledseeds|>offsetofmodes,1)
#     emptyseedsfock = quantize(emptyseeds|>offsetofmodes,1)
#     courierseedsfock = quantize(courierseeds|>offsetofmodes,1)

#     noooffilledmodes = length(filledseeds)
#     nooofemptymodes = length(emptyseeds)
#     noofcouriermodes = length(courierseeds)
    
#     localstates = getregionstates(localcorrelations=localcorrelations, grouping=[noooffilledmodes,noofcouriermodes,nooofemptymodes])
#     localcourieriso = localstates[2]|>FockMap
#     localfillediso = localstates[1]|>FockMap
#     localemptyiso = localstates[3]|>FockMap

#     iden = idmap(localcorrelations|>getoutspace, localcorrelations|>getinspace)
#     localfilledseedsiso = iden[:,filledseedsfock]
#     localemptyseedsiso = iden[:,emptyseedsfock]
#     localcourierseedsiso = iden[:,courierseedsfock]


#     localwannierfillediso,minsvdempty = localwannierization(localfillediso, localfilledseedsiso)
#     localwannieremptyiso,minsvdfilled = localwannierization(localemptyiso, localemptyseedsiso)
#     localwanniercourieriso,minsvdcourier = localwannierization(localcourieriso, localcourierseedsiso)
#     localdisentangler = localwanniercourieriso+localwannierfillediso+localwannieremptyiso
#     redefinelocaldisentangler = FockMap(localdisentangler|>getoutspace,localdisentangler|>getinspace|>RegionFock,localdisentangler|>rep)

#     newcrystalfock = getcrystalfock(redefinelocaldisentangler|>getinspace|>unitcellfock, firstgmeracrystal)
#     newlocalrestrict = fourier(newcrystalfock, redefinelocaldisentangler|>getinspace) * (firstgmeracrystal|>vol|>sqrt)

#     globaldisentangler = broadcast(*,(localrestrict*redefinelocaldisentangler), newlocalrestrict')
#     secondgmeracorrelations = globaldisentangler'*firstgmeracorrelations*globaldisentangler
#     secondgmeraH = globaldisentangler'*firstgmeraH*globaldisentangler

#     return Dict(
#         :localcorrelations=>localcorrelations,
#         :globaldisentangler=>globaldisentangler,
#         :secondgmeracorrelations=>secondgmeracorrelations,
#         :secondgmeraH=>secondgmeraH,
#         :minsvdcourier=>minsvdcourier,
#         :minsvdempty=>minsvdempty,
#         :minsvdfilled=>minsvdfilled)
#     end
# export gmerasecondstep

# function gmerathirdstep(secondgmeracorrelations,secondgmeraH)
#     secondgmeracrystalfock = secondgmeracorrelations|>getoutspace
#     secondgmeracrystal::Crystal = secondgmeracrystalfock|>getcrystal
#     secondgmeraspace::RealSpace = secondgmeracrystal|>getspace

#     refrot = inv([2/3 -1/3; -1/3 2/3]')
#     c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
#     c3 = c6^2

#     thirdcenter = [1/3,1/3] ∈ secondgmeraspace
#     thirdhexagonalregion = gethexagonalregion(rot = refrot,crystal=secondgmeracrystal, center=thirdcenter, metricspace=secondgmeraspace)
#     thirdhexagonalregionfock = quantize(thirdhexagonalregion,1)
#     thirdhexagonalmodes = Subset(mode for mode in thirdhexagonalregionfock)

#     c6recenter = recenter(c6,thirdcenter)
#     c3recenter = recenter(c3,thirdcenter)

#     # @info "Computing local correlations..."
#     localrestrict = fourier(secondgmeracrystalfock, thirdhexagonalregionfock) / (secondgmeracrystal|>vol|>sqrt)
#     localcorrelations = localrestrict'*secondgmeracorrelations*localrestrict
#     potentialfilledseeds = Subset([mode for mode in thirdhexagonalregionfock if (localcorrelations[mode,mode]|>rep)[1,1]|>abs<0.5])
#     potentialemptyseeds = Subset([mode for mode in thirdhexagonalregionfock if (localcorrelations[mode,mode]|>rep)[1,1]|>abs>0.5])

#     localspectrum = localcorrelations|>eigspech
#     localspectrum|>visualize|>display
#     entanglementcontourinfo = entanglementcontour(localcorrelations)
#     groupedandsortedmodeswifecontour = sortgroupdictwifvalue(entanglementcontourinfo,false)
#     nooffrozenmodes = div(thirdhexagonalregionfock|>length,4)*3
#     frozenseeds = groupedandsortedmodeswifecontour[1][2]
#     for (entanglementval,modes) in groupedandsortedmodeswifecontour[2:end]
#         if (frozenseeds|>length)<nooffrozenmodes
#             @info "include more seeds"
#             frozenseeds = frozenseeds+modes
#         elseif (frozenseeds|>length)==nooffrozenmodes
#             @info "just enough frozenseeds"
#             break
#         else 
#             @error "too many frozenseeds, sth wrong"
#         end
#     end
#     filledseeds = intersect(frozenseeds,potentialfilledseeds)
#     emptyseeds = intersect(frozenseeds,potentialemptyseeds)
#     courierseeds = thirdhexagonalmodes - frozenseeds
#     filledseedsfock = quantize(filledseeds|>offsetofmodes,1)
#     emptyseedsfock = quantize(emptyseeds|>offsetofmodes,1)
#     courierseedsfock = quantize(courierseeds|>offsetofmodes,1)

#     noooffilledmodes = length(filledseeds)
#     nooofemptymodes = length(emptyseeds)
#     noofcouriermodes = length(courierseeds)

#     localstates = getregionstates(localcorrelations=localcorrelations, grouping=[noooffilledmodes,noofcouriermodes,nooofemptymodes])
#     localcourierisometry = localstates[2]|>FockMap
#     localfilledisometry = localstates[1]|>FockMap
#     localemptyisometry = localstates[3]|>FockMap

#     iden = idmap(localcorrelations|>getoutspace, localcorrelations|>getinspace)
#     localfilledseedsisometry = iden[:,filledseedsfock]
#     localemptyseedsisometry = iden[:,emptyseedsfock]
#     localcourierseedsisometry = iden[:,courierseedsfock]

#     localwannierfilled,minsvdempty = localwannierization(localfilledisometry, localfilledseedsisometry)
#     localwannierempty,minsvdfilled = localwannierization(localemptyisometry, localemptyseedsisometry)
#     localwanniercourier,minsvdcourier = localwannierization(localcourierisometry, localcourierseedsisometry)
#     localdisentangler = localwanniercourier+localwannierfilled+localwannierempty
#     redefinelocaldisentangler = FockMap(localdisentangler|>getoutspace,localdisentangler|>getinspace|>RegionFock,localdisentangler|>rep)

#     newcrystalfock = getcrystalfock(redefinelocaldisentangler|>getinspace|>unitcellfock, secondgmeracrystal)

#     newlocalrestrict = fourier(newcrystalfock, redefinelocaldisentangler|>getinspace) * (secondgmeracrystal|>vol|>sqrt)
#     globaldisentangler = broadcast(*,(localrestrict*redefinelocaldisentangler), newlocalrestrict')
#     thirdgmeracorrelations = globaldisentangler'*secondgmeracorrelations*globaldisentangler
#     thirdgmeraH = globaldisentangler'*secondgmeraH*globaldisentangler


#     return return Dict(
#         :localcorrelations=>localcorrelations,
#         :globaldisentangler=>globaldisentangler,
#         :thirdgmeracorrelations=>thirdgmeracorrelations,
#         :thirdgmeraH=>thirdgmeraH,
#         :minsvdcourier=>minsvdcourier,
#         :minsvdempty=>minsvdempty,
#         :minsvdfilled=>minsvdfilled)
#     end
# export gmerathirdstep

function gmeraredefinestep(blockedcorrelations,blockedH,noofflavourpermode::Number)
    blockedcrystalfock = blockedcorrelations|>getoutspace
    blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal
    blockedspace::RealSpace = blockedcrystal|>getspace

    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    origin = [0,0] ∈ blockedspace
    hexagonalregion = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=origin, metricspace=blockedspace)
    hexagonalregionfock = quantize(hexagonalregion,noofflavourpermode)
    rgshiftedcenter1 = [2/3,-1/3] ∈ blockedspace
    rgshiftedhexagonalregion1 = hexagonalregion.+rgshiftedcenter1
    rgshiftedcenter2 = [1/3,1/3] ∈ blockedspace
    rgshiftedhexagonalregion2 = hexagonalregion.+rgshiftedcenter2

    siteAregion1 = intersect(intersect(hexagonalregion,rgshiftedhexagonalregion1 ),rgshiftedhexagonalregion2)
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

    # @info "Computing local correlations..."
    localrestrict = fourier(blockedcrystalfock, hexagonalregionfock) / (blockedcrystal|>vol|>sqrt)
    localcorrelations = localrestrict'*blockedcorrelations*localrestrict
    localspectrum = localcorrelations|>eigspec
    display(localspectrum|>visualize)

    totalnumberofmodes = hexagonalregionfock|>length

    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[totalnumberofmodes])
    localunitary = localstates[1]|>FockMap

    numberofmodespercorner = div(totalnumberofmodes,6)

    localA1restrict = fourier(blockedcrystalfock, siteAregionfock1) / (blockedcrystal|>vol|>sqrt)
    localAcorrelations1 = localA1restrict'*blockedcorrelations*localA1restrict
    localAstates1 = getregionstates(localcorrelations=localAcorrelations1, grouping=[numberofmodespercorner])
    localAseeds1 = localAstates1[1]

    localA2restrict = fourier(blockedcrystalfock, siteAregionfock2) / (blockedcrystal|>vol|>sqrt)
    localAcorrelations2 = localA2restrict'*blockedcorrelations*localA2restrict
    localAstates2 = getregionstates(localcorrelations=localAcorrelations2, grouping=[numberofmodespercorner])
    localAseeds2 = localAstates2[1]


    localA3restrict = fourier(blockedcrystalfock, siteAregionfock3) / (blockedcrystal|>vol|>sqrt)
    localAcorrelations3 = localA3restrict'*blockedcorrelations*localA3restrict
    localAstates3 = getregionstates(localcorrelations=localAcorrelations3, grouping=[numberofmodespercorner])
    localAseeds3 = localAstates3[1]

    allAseedsstate = localAseeds1+localAseeds2+localAseeds3

    localB1restrict = fourier(blockedcrystalfock, siteBregionfock1) / (blockedcrystal|>vol|>sqrt)
    localBcorrelations1 = localB1restrict'*blockedcorrelations*localB1restrict
    localBstates1 = getregionstates(localcorrelations=localBcorrelations1, grouping=[numberofmodespercorner])
    localBseeds1 = localBstates1[1]

    localB2restrict = fourier(blockedcrystalfock, siteBregionfock2) / (blockedcrystal|>vol|>sqrt)
    localBcorrelations2 = localB2restrict'*blockedcorrelations*localB2restrict
    localBstates2 = getregionstates(localcorrelations=localBcorrelations2, grouping=[numberofmodespercorner])
    localBseeds2 = localBstates2[1]

    localB3restrict = fourier(blockedcrystalfock, siteBregionfock3) / (blockedcrystal|>vol|>sqrt)
    localBcorrelations3 = localB3restrict'*blockedcorrelations*localB3restrict
    localBstates3 = getregionstates(localcorrelations=localBcorrelations3, grouping=[numberofmodespercorner])
    localBseeds3 = localBstates3[1]

    allBseedsstate = localBseeds1+localBseeds2+localBseeds3

    allseedsstates = (allAseedsstate+allBseedsstate)
    allseeds = allseedsstates|>FockMap

    localwannierunitary,minsvd = localwannierization(localunitary, allseeds)

    redunitcellfock = localwannierunitary|>getinspace|>unitcellfock
    redunitcellmodes = Subset(m for m in redunitcellfock)
    newblockedcrystalfock = getcrystalfock(redunitcellfock, Crystal(redunitcellmodes|>basisoffsetofmodes,blockedcrystal.sizes))

    newlocalrestrict = fourier(newblockedcrystalfock, localwannierunitary|>getinspace) * (blockedcrystal|>vol|>sqrt)

    globalunitary = broadcast(*,(localrestrict*localwannierunitary), newlocalrestrict')

    rotcorrelations = globalunitary'*blockedcorrelations*globalunitary
    rotH = globalunitary'*blockedH*globalunitary

    return Dict(
        :localcorrelations=>localcorrelations,
        :globalunitary=>globalunitary,
        :rotcorrelations=>rotcorrelations,
        :rotH=>rotH,
        :minsvdredef=>minsvd)
end
export gmeraredefinestep

function gmerathirdstep(rotcorrelations,rotH,noofcouriermodes::Number,noofflavourpermode::Number)
    rotcrystalfock = rotcorrelations|>getoutspace
    rotcrystal = rotcorrelations|>getoutspace|>getcrystal
    rotspace::RealSpace = rotcrystal|>getspace
        
    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    firstcenter = [1/3,1/3] ∈ rotspace
    firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=rotcrystal, center=firstcenter, metricspace=rotspace)
    firsthexagonalregionfock = quantize(firsthexagonalregion,noofflavourpermode)
    rgshiftedcenter1 = [2/3,-1/3] ∈ rotspace
    firstrgshiftedhexagonalregion1 = firsthexagonalregion.+rgshiftedcenter1
    rgshiftedcenter2 = [1/3,1/3] ∈ rotspace
    firstrgshiftedhexagonalregion2 = firsthexagonalregion.+rgshiftedcenter2

    c6recenter = recenter(c6,firstcenter)
    c3recenter = recenter(c3,firstcenter)

    siteAregion1 = intersect(intersect(firsthexagonalregion,firstrgshiftedhexagonalregion1 ),firstrgshiftedhexagonalregion2)
    siteAregionfock1 = quantize(siteAregion1,noofflavourpermode)
    siteAregion2 = c3recenter*siteAregion1
    siteAregionfock2 = quantize(siteAregion2,noofflavourpermode)
    siteAregion3 = (c3recenter)*siteAregion2
    siteAregionfock3 = quantize(siteAregion3,noofflavourpermode)

    siteBregion1 = (c6recenter)*siteAregion1
    siteBregionfock1 = quantize(siteBregion1,noofflavourpermode)
    siteBregion2 = (c3recenter)*siteBregion1
    siteBregionfock2 = quantize(siteBregion2,noofflavourpermode)
    siteBregion3 = (c3recenter)*siteBregion2
    siteBregionfock3 = quantize(siteBregion3,noofflavourpermode)

    nooffrozenmodes = length(firsthexagonalregionfock|>orderedmodes)-noofcouriermodes
    nooffilledmodes = div(nooffrozenmodes,2)
    noofemptymodes = nooffilledmodes
    @info("no of filled modes = ",nooffilledmodes)
    @info("no of courier modes = ",noofcouriermodes)
    @info("no of empty modes = ",noofemptymodes)

    # @info "Computing local correlations..."
    localrestrict = fourier(rotcrystalfock, firsthexagonalregionfock) / (rotcrystal|>vol|>sqrt)
    localcorrelations = localrestrict'*rotcorrelations*localrestrict

    localspectrum = localcorrelations|>eigspec
    localspectrum|>visualize|>display

    if noofcouriermodes != 0 
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

        localwanniercourierisometry,minsvd = localwannierization(localcourierisometry, courierseeds)
        courierunitcellfock = localwanniercourierisometry|>getinspace|>unitcellfock
        courierunitcellmodes = Subset(m for m in courierunitcellfock)

        filledunitcellfock = localfilledisometry|>getinspace|>unitcellfock
        filledunitcellmodes = Subset(m for m in filledunitcellfock)
        emptyunitcellfock = localemptyisometry|>getinspace|>unitcellfock
        emptyunitcellmodes = Subset(m for m in emptyunitcellfock)

        newcouriercrystalfock = getcrystalfock(localwanniercourierisometry|>getinspace|>unitcellfock, Crystal(courierunitcellmodes|>basisoffsetofmodes,rotcrystal.sizes))
        newlocalcourierrestrict = fourier(newcouriercrystalfock, localwanniercourierisometry|>getinspace) * (rotcrystal|>vol|>sqrt)

        newfilledcrystalfock = getcrystalfock(localfilledisometry|>getinspace|>unitcellfock, Crystal(filledunitcellmodes|>basisoffsetofmodes,rotcrystal.sizes))
        newlocalfilledrestrict = fourier(newfilledcrystalfock, localfilledisometry|>getinspace) * (rotcrystal|>vol|>sqrt)
        newemptycrystalfock = getcrystalfock(localemptyisometry|>getinspace|>unitcellfock, Crystal(emptyunitcellmodes|>basisoffsetofmodes,rotcrystal.sizes))
        newlocalemptyrestrict = fourier(newemptycrystalfock, localemptyisometry|>getinspace) * (rotcrystal|>vol|>sqrt)

        globalcourierisometry = broadcast(*,(localrestrict*localwanniercourierisometry), newlocalcourierrestrict')
        globalfilledisometry = broadcast(*,(localrestrict*localfilledisometry), newlocalfilledrestrict')
        globalemptyisometry = broadcast(*,(localrestrict*localemptyisometry), newlocalemptyrestrict')

        rawcouriercorrelations = globalcourierisometry'*rotcorrelations*globalcourierisometry
        couriercorrelationspectrum = rawcouriercorrelations|>crystalspectrum
        couriercorrelations = roundingpurification(couriercorrelationspectrum)|>crystalfockmap
        courierH = globalcourierisometry'*rotH*globalcourierisometry

        filledcorrelations = globalfilledisometry'*rotcorrelations*globalfilledisometry
        filledH = globalfilledisometry'*rotH*globalfilledisometry
        emptycorrelations = globalemptyisometry'*rotcorrelations*globalemptyisometry
        emptyH = globalemptyisometry'*rotH*globalemptyisometry

        return Dict(
            :localcorrelations=>localcorrelations,
            :globalcourierisometry=>globalcourierisometry,
            :globalfilledisometry=>globalfilledisometry,
            :globalemptyisometry=>globalemptyisometry,
            :couriercorrelations=>couriercorrelations,
            :courierH=>courierH,
            :filledcorrelations=>filledcorrelations,
            :filledH=>filledH,
            :emptycorrelations=>emptycorrelations,
            :emptyH=>emptyH,
            :rawcouriercorrelations=>rawcouriercorrelations,
            :minsvdcourier=>minsvd)
    else
        localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofemptymodes])
        localfilledisometry = localstates[1]|>FockMap
        localemptyisometry = localstates[2]|>FockMap

        nooflocalcourierseeds = div(noofcouriermodes,6)
        reduandancy = length((siteAregionfock1|>orderedmodes))-nooflocalcourierseeds


        filledunitcellfock = localfilledisometry|>getinspace|>unitcellfock
        filledunitcellmodes = Subset(m for m in filledunitcellfock)
        emptyunitcellfock = localemptyisometry|>getinspace|>unitcellfock
        emptyunitcellmodes = Subset(m for m in emptyunitcellfock)

        newfilledcrystalfock = getcrystalfock(localfilledisometry|>getinspace|>unitcellfock, Crystal(filledunitcellmodes|>basisoffsetofmodes,rotcrystal.sizes))
        newlocalfilledrestrict = fourier(newfilledcrystalfock, localfilledisometry|>getinspace) * (rotcrystal|>vol|>sqrt)
        newemptycrystalfock = getcrystalfock(localemptyisometry|>getinspace|>unitcellfock, Crystal(emptyunitcellmodes|>basisoffsetofmodes,rotcrystal.sizes))
        newlocalemptyrestrict = fourier(newemptycrystalfock, localemptyisometry|>getinspace) * (rotcrystal|>vol|>sqrt)

        globalfilledisometry = broadcast(*,(localrestrict*localfilledisometry), newlocalfilledrestrict')
        globalemptyisometry = broadcast(*,(localrestrict*localemptyisometry), newlocalemptyrestrict')

        filledcorrelations = globalfilledisometry'*rotcorrelations*globalfilledisometry
        filledH = globalfilledisometry'*rotH*globalfilledisometry
        emptycorrelations = globalemptyisometry'*rotcorrelations*globalemptyisometry
        emptyH = globalemptyisometry'*rotH*globalemptyisometry

        return Dict(
            :localcorrelations=>localcorrelations,
            :globalfilledisometry=>globalfilledisometry,
            :globalemptyisometry=>globalemptyisometry,
            :filledcorrelations=>filledcorrelations,
            :filledH=>filledH,
            :emptycorrelations=>emptycorrelations,
            :emptyH=>emptyH)
    end
end
export gmerathirdstep

function gmerasecondstep(firstgmeracorrelations,firstgmeraH,noofcouriermodes::Number,noofflavourpermode::Number)
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

    # @info "Computing local correlations..."
    localrestrict = fourier(firstgmeracrystalfock, secondhexagonalregionfock) / (firstgmeracrystal|>vol|>sqrt)
    localcorrelations = localrestrict'*firstgmeracorrelations*localrestrict

    localspectrum = localcorrelations|>eigspec
    localspectrum|>visualize|>display

    if noofcouriermodes != 0 
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

        localwanniercourierisometry,minsvd = localwannierization(localcourierisometry, courierseeds)
        courierunitcellfock = localwanniercourierisometry|>getinspace|>unitcellfock
        courierunitcellmodes = Subset(m for m in courierunitcellfock)

        filledunitcellfock = localfilledisometry|>getinspace|>unitcellfock
        filledunitcellmodes = Subset(m for m in filledunitcellfock)
        emptyunitcellfock = localemptyisometry|>getinspace|>unitcellfock
        emptyunitcellmodes = Subset(m for m in emptyunitcellfock)

        newcouriercrystalfock = getcrystalfock(localwanniercourierisometry|>getinspace|>unitcellfock, Crystal(courierunitcellmodes|>basisoffsetofmodes,firstgmeracrystal.sizes))
        newlocalcourierrestrict = fourier(newcouriercrystalfock, localwanniercourierisometry|>getinspace) * (firstgmeracrystal|>vol|>sqrt)

        newfilledcrystalfock = getcrystalfock(localfilledisometry|>getinspace|>unitcellfock, Crystal(filledunitcellmodes|>basisoffsetofmodes,firstgmeracrystal.sizes))
        newlocalfilledrestrict = fourier(newfilledcrystalfock, localfilledisometry|>getinspace) * (firstgmeracrystal|>vol|>sqrt)
        newemptycrystalfock = getcrystalfock(localemptyisometry|>getinspace|>unitcellfock, Crystal(emptyunitcellmodes|>basisoffsetofmodes,firstgmeracrystal.sizes))
        newlocalemptyrestrict = fourier(newemptycrystalfock, localemptyisometry|>getinspace) * (firstgmeracrystal|>vol|>sqrt)

        globalcourierisometry = broadcast(*,(localrestrict*localwanniercourierisometry), newlocalcourierrestrict')
        globalfilledisometry = broadcast(*,(localrestrict*localfilledisometry), newlocalfilledrestrict')
        globalemptyisometry = broadcast(*,(localrestrict*localemptyisometry), newlocalemptyrestrict')

        rawcouriercorrelations = globalcourierisometry'*firstgmeracorrelations*globalcourierisometry
        couriercorrelationspectrum = rawcouriercorrelations|>crystalspectrum
        couriercorrelations = roundingpurification(couriercorrelationspectrum)|>crystalfockmap
        courierH = globalcourierisometry'*firstgmeraH*globalcourierisometry

        filledcorrelations = globalfilledisometry'*firstgmeracorrelations*globalfilledisometry
        filledH = globalfilledisometry'*firstgmeraH*globalfilledisometry
        emptycorrelations = globalemptyisometry'*firstgmeracorrelations*globalemptyisometry
        emptyH = globalemptyisometry'*firstgmeraH*globalemptyisometry

        return Dict(
            :localcorrelations=>localcorrelations,
            :globalcourierisometry=>globalcourierisometry,
            :globalfilledisometry=>globalfilledisometry,
            :globalemptyisometry=>globalemptyisometry,
            :couriercorrelations=>couriercorrelations,
            :courierH=>courierH,
            :filledcorrelations=>filledcorrelations,
            :filledH=>filledH,
            :emptycorrelations=>emptycorrelations,
            :emptyH=>emptyH,
            :rawcouriercorrelations=>rawcouriercorrelations,
            :minsvdcourier=>minsvd)
    else
        localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofemptymodes])
        localfilledisometry = localstates[1]|>FockMap
        localemptyisometry = localstates[2]|>FockMap

        filledunitcellfock = localfilledisometry|>getinspace|>unitcellfock
        filledunitcellmodes = Subset(m for m in filledunitcellfock)
        emptyunitcellfock = localemptyisometry|>getinspace|>unitcellfock
        emptyunitcellmodes = Subset(m for m in emptyunitcellfock)

        newfilledcrystalfock = getcrystalfock(localfilledisometry|>getinspace|>unitcellfock, Crystal(filledunitcellmodes|>basisoffsetofmodes,firstgmeracrystal.sizes))
        newlocalfilledrestrict = fourier(newfilledcrystalfock, localfilledisometry|>getinspace) * (firstgmeracrystal|>vol|>sqrt)
        newemptycrystalfock = getcrystalfock(localemptyisometry|>getinspace|>unitcellfock, Crystal(emptyunitcellmodes|>basisoffsetofmodes,firstgmeracrystal.sizes))
        newlocalemptyrestrict = fourier(newemptycrystalfock, localemptyisometry|>getinspace) * (firstgmeracrystal|>vol|>sqrt)

        globalfilledisometry = broadcast(*,(localrestrict*localfilledisometry), newlocalfilledrestrict')
        globalemptyisometry = broadcast(*,(localrestrict*localemptyisometry), newlocalemptyrestrict')

        filledcorrelations = globalfilledisometry'*firstgmeracorrelations*globalfilledisometry
        filledH = globalfilledisometry'*firstgmeraH*globalfilledisometry
        emptycorrelations = globalemptyisometry'*firstgmeracorrelations*globalemptyisometry
        emptyH = globalemptyisometry'*firstgmeraH*globalemptyisometry

        return Dict(
            :localcorrelations=>localcorrelations,
            :globalfilledisometry=>globalfilledisometry,
            :globalemptyisometry=>globalemptyisometry,
            :filledcorrelations=>filledcorrelations,
            :filledH=>filledH,
            :emptycorrelations=>emptycorrelations,
            :emptyH=>emptyH)
    end 
end
export gmerasecondstep

function gmerafirststep(secondgmeracorrelations,secondgmeraH,noofcouriermodes::Number,noofflavourpermode::Number)
    secondgmeracrystalfock = secondgmeracorrelations|>getoutspace
    secondgmeracrystal::Crystal = secondgmeracrystalfock|>getcrystal
    secondgmeraspace::RealSpace = secondgmeracrystal|>getspace

    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    thirdcenter = [0,0] ∈ secondgmeraspace
    thirdhexagonalregion = gethexagonalregion(rot = refrot,crystal=secondgmeracrystal, center=thirdcenter, metricspace=secondgmeraspace)
    thirdhexagonalregionfock = quantize(thirdhexagonalregion,noofflavourpermode)
    rgshiftedcenter1 = [2/3,-1/3] ∈ secondgmeraspace
    thirdrgshiftedhexagonalregion1 = thirdhexagonalregion.+rgshiftedcenter1
    rgshiftedcenter2 = [1/3,1/3] ∈ secondgmeraspace
    thirdrgshiftedhexagonalregion2 = thirdhexagonalregion.+rgshiftedcenter2

    siteAregion1 = intersect(intersect(thirdhexagonalregion,thirdrgshiftedhexagonalregion1),thirdrgshiftedhexagonalregion2)
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

    nooffrozenmodes = length(thirdhexagonalregionfock|>orderedmodes)-noofcouriermodes
    nooffilledmodes = div(nooffrozenmodes,2)
    noofemptymodes = nooffilledmodes
    @info("no of filled modes = ",nooffilledmodes)
    @info("no of courier modes = ",noofcouriermodes)
    @info("no of empty modes = ",noofemptymodes)

    # @info "Computing local correlations..."
    localrestrict = fourier(secondgmeracrystalfock, thirdhexagonalregionfock) / (secondgmeracrystal|>vol|>sqrt)
    localcorrelations = localrestrict'*secondgmeracorrelations*localrestrict

    localspectrum = localcorrelations|>eigspec
    localspectrum|>visualize|>display

    if noofcouriermodes != 0 
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

        localwanniercourierisometry,minsvd = localwannierization(localcourierisometry, courierseeds)
        courierunitcellfock = localwanniercourierisometry|>getinspace|>unitcellfock
        courierunitcellmodes = Subset(m for m in courierunitcellfock)

        filledunitcellfock = localfilledisometry|>getinspace|>unitcellfock
        filledunitcellmodes = Subset(m for m in filledunitcellfock)
        emptyunitcellfock = localemptyisometry|>getinspace|>unitcellfock
        emptyunitcellmodes = Subset(m for m in emptyunitcellfock)

        newcouriercrystalfock = getcrystalfock(localwanniercourierisometry|>getinspace|>unitcellfock, Crystal(courierunitcellmodes|>basisoffsetofmodes,secondgmeracrystal.sizes))
        newlocalcourierrestrict = fourier(newcouriercrystalfock, localwanniercourierisometry|>getinspace) * (secondgmeracrystal|>vol|>sqrt)

        newfilledcrystalfock = getcrystalfock(localfilledisometry|>getinspace|>unitcellfock, Crystal(filledunitcellmodes|>basisoffsetofmodes,secondgmeracrystal.sizes))
        newlocalfilledrestrict = fourier(newfilledcrystalfock, localfilledisometry|>getinspace) * (secondgmeracrystal|>vol|>sqrt)
        newemptycrystalfock = getcrystalfock(localemptyisometry|>getinspace|>unitcellfock, Crystal(emptyunitcellmodes|>basisoffsetofmodes,secondgmeracrystal.sizes))
        newlocalemptyrestrict = fourier(newemptycrystalfock, localemptyisometry|>getinspace) * (secondgmeracrystal|>vol|>sqrt)

        globalcourierisometry = broadcast(*,(localrestrict*localwanniercourierisometry), newlocalcourierrestrict')
        globalfilledisometry = broadcast(*,(localrestrict*localfilledisometry), newlocalfilledrestrict')
        globalemptyisometry = broadcast(*,(localrestrict*localemptyisometry), newlocalemptyrestrict')

        rawcouriercorrelations = globalcourierisometry'*secondgmeracorrelations*globalcourierisometry
        couriercorrelationspectrum = rawcouriercorrelations|>crystalspectrum
        couriercorrelations = roundingpurification(couriercorrelationspectrum)|>crystalfockmap
        courierH = globalcourierisometry'*secondgmeraH*globalcourierisometry

        filledcorrelations = globalfilledisometry'*secondgmeracorrelations*globalfilledisometry
        filledH = globalfilledisometry'*secondgmeraH*globalfilledisometry
        emptycorrelations = globalemptyisometry'*secondgmeracorrelations*globalemptyisometry
        emptyH = globalemptyisometry'*secondgmeraH*globalemptyisometry


        return Dict(
            :localcorrelations=>localcorrelations,
            :globalcourierisometry=>globalcourierisometry,
            :globalfilledisometry=>globalfilledisometry,
            :globalemptyisometry=>globalemptyisometry,
            :couriercorrelations=>couriercorrelations,
            :courierH=>courierH,
            :filledcorrelations=>filledcorrelations,
            :filledH=>filledH,
            :emptycorrelations=>emptycorrelations,
            :emptyH=>emptyH,
            :rawcouriercorrelations=>rawcouriercorrelations,
            :minsvdcourier=>minsvd)
    else
        localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofemptymodes])
        localfilledisometry = localstates[1]|>FockMap
        localemptyisometry = localstates[2]|>FockMap


        filledunitcellfock = localfilledisometry|>getinspace|>unitcellfock
        filledunitcellmodes = Subset(m for m in filledunitcellfock)
        emptyunitcellfock = localemptyisometry|>getinspace|>unitcellfock
        emptyunitcellmodes = Subset(m for m in emptyunitcellfock)

        newfilledcrystalfock = getcrystalfock(localfilledisometry|>getinspace|>unitcellfock, Crystal(filledunitcellmodes|>basisoffsetofmodes,secondgmeracrystal.sizes))
        newlocalfilledrestrict = fourier(newfilledcrystalfock, localfilledisometry|>getinspace) * (secondgmeracrystal|>vol|>sqrt)
        newemptycrystalfock = getcrystalfock(localemptyisometry|>getinspace|>unitcellfock, Crystal(emptyunitcellmodes|>basisoffsetofmodes,secondgmeracrystal.sizes))
        newlocalemptyrestrict = fourier(newemptycrystalfock, localemptyisometry|>getinspace) * (secondgmeracrystal|>vol|>sqrt)

        globalfilledisometry = broadcast(*,(localrestrict*localfilledisometry), newlocalfilledrestrict')
        globalemptyisometry = broadcast(*,(localrestrict*localemptyisometry), newlocalemptyrestrict')

        filledcorrelations = globalfilledisometry'*secondgmeracorrelations*globalfilledisometry
        filledH = globalfilledisometry'*secondgmeraH*globalfilledisometry
        emptycorrelations = globalemptyisometry'*secondgmeracorrelations*globalemptyisometry
        emptyH = globalemptyisometry'*secondgmeraH*globalemptyisometry


        return Dict(
            :localcorrelations=>localcorrelations,
            :globalfilledisometry=>globalfilledisometry,
            :globalemptyisometry=>globalemptyisometry,
            :filledcorrelations=>filledcorrelations,
            :filledH=>filledH,
            :emptycorrelations=>emptycorrelations,
            :emptyH=>emptyH)
    end
end
export gmerafirststep

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

        # @info "Computing local correlations..."
        localrestrict = fourier(crystalfock, hexagonalregionfock) / (crystal|>vol|>sqrt)
        localcorrelations = localrestrict'*correlations*localrestrict
        localspectrum = localcorrelations|>eigspech
        localspectrum|>visualize|>display

        localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofemptymodes])
        localfilledisometry = localstates[1]|>FockMap
        localemptyisometry = localstates[2]|>FockMap

        filledunitcellfock = localfilledisometry|>getinspace|>unitcellfock
        filledunitcellmodes = Subset(m for m in filledunitcellfock)
        emptyunitcellfock = localemptyisometry|>getinspace|>unitcellfock
        emptyunitcellmodes = Subset(m for m in emptyunitcellfock)

        
        newfilledcrystalfock = getcrystalfock(localfilledisometry|>getinspace|>unitcellfock, Crystal(filledunitcellmodes|>basisoffsetofmodes,crystal.sizes))
        newemptycrystalfock = getcrystalfock(localemptyisometry|>getinspace|>unitcellfock, Crystal(emptyunitcellmodes|>basisoffsetofmodes,crystal.sizes))

        newlocalfilledrestrict = fourier(newfilledcrystalfock, localfilledisometry|>getinspace) * (crystal|>vol|>sqrt)
        newlocalemptyrestrict = fourier(newemptycrystalfock, localemptyisometry|>getinspace) * (crystal|>vol|>sqrt)

        globalfilledisometry = broadcast(*,(localrestrict*localfilledisometry), newlocalfilledrestrict')
        globalemptyisometry = broadcast(*,(localrestrict*localemptyisometry), newlocalemptyrestrict')

        filledcorrelations = globalfilledisometry'*correlations*globalfilledisometry
        emptycorrelations = globalemptyisometry'*correlations*globalemptyisometry

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

function firstgmerastep(correlations,H,scaling,systemsize)
    blockedcorrelations,blockedH,blocker = blocking(correlations,H,scaling)
    fiosave(blockedcorrelations, name="blockedcorrelations")
    fiosave(blockedH, name="blockedH")
    fiosave(blocker, name="blocker")

    noofmodesinlocalreg = (6*(scaling^2))|>Int
    noofdistillablemodes = ((1/4)*noofmodesinlocalreg)|>Int
    noofcouriermodesinfirststep = noofmodesinlocalreg - noofdistillablemodes
    # println(noofcouriermodesinfirststep)
    noofcouriermodesinsecondstep = noofcouriermodesinfirststep - noofdistillablemodes
    noofcouriermodesinthirdstep = noofcouriermodesinsecondstep - noofdistillablemodes

    # gmera1redefinedata = gmeraredefinestep(blockedcorrelations,blockedH,1)
    gmera1firststepdata = gmerafirststep(blockedcorrelations,blockedH,noofcouriermodesinfirststep,1)
    gmera1firststepterminateddata = gmerafirststep(blockedcorrelations,blockedH,0,1)

    gmera1secondstepdata = gmerasecondstep(gmera1firststepdata[:couriercorrelations],gmera1firststepdata[:courierH],noofcouriermodesinsecondstep,div(noofcouriermodesinfirststep,6))
    gmera1secondstepterminateddata = gmerasecondstep(gmera1firststepdata[:couriercorrelations],gmera1firststepdata[:courierH],0,div(noofcouriermodesinfirststep,6))

    gmera1thirdstepdata = gmerathirdstep(gmera1secondstepdata[:couriercorrelations],gmera1secondstepdata[:courierH],noofcouriermodesinthirdstep,div(noofcouriermodesinsecondstep,6))
    gmera1thirdsteterminateddata = gmerathirdstep(gmera1secondstepdata[:couriercorrelations],gmera1secondstepdata[:courierH],0,div(noofcouriermodesinsecondstep,6))

    noofflavourpermodeforlaterrg = (noofcouriermodesinthirdstep/6)|>Int

    gmera1firstapproximation = (gmera1firststepdata[:globalemptyisometry])*(gmera1firststepdata[:globalemptyisometry])'
    gmera1firstterminatedapproximation = (gmera1firststepterminateddata[:globalemptyisometry])*(gmera1firststepterminateddata[:globalemptyisometry])'
    fiosave(gmera1firstterminatedapproximation, name="gmera1firstterminatedapproximation")

    gmera1firstalldiff = blockedcorrelations - gmera1firstterminatedapproximation
    fiosave(gmera1firstalldiff, name="gmera1firstalldiff")

    @info("calculating L1 norm of the difference between the blocked correlations and the terminated approximate correlation at gmera step1")
    gmera1firstalldiffL1norm = focktraceL1norm(gmera1firstalldiff,systemsize^2*6)
    fiosave(gmera1firstalldiffL1norm, name="gmera1firstalldiffL1norm")

    gmera1secondapproximation = (gmera1firststepdata[:globalcourierisometry]*gmera1secondstepdata[:globalemptyisometry])*(gmera1firststepdata[:globalcourierisometry]*gmera1secondstepdata[:globalemptyisometry])'
    gmera1secondterminatedapproximation = gmera1firstapproximation + (gmera1firststepdata[:globalcourierisometry]*gmera1secondstepterminateddata[:globalemptyisometry])*(gmera1firststepdata[:globalcourierisometry]*gmera1secondstepterminateddata[:globalemptyisometry])'
    fiosave(gmera1secondterminatedapproximation, name="gmera1secondterminatedapproximation")

    gmera1secondalldiff = blockedcorrelations - gmera1secondterminatedapproximation 
    fiosave(gmera1secondalldiff, name="gmera1secondalldiff")

    @info("calculating L1 norm of the difference between the blocked correlations and the terminated approximate correlation at gmera step2")
    gmera1secondalldiffL1norm = focktraceL1norm(gmera1secondalldiff,systemsize^2*6)
    fiosave(gmera1secondalldiffL1norm, name="gmera1secondalldiffL1norm")

    gmera1thirdapproximation = (gmera1firststepdata[:globalcourierisometry]*gmera1secondstepdata[:globalcourierisometry]*gmera1thirdstepdata[:globalemptyisometry])*(gmera1firststepdata[:globalcourierisometry]*gmera1secondstepdata[:globalcourierisometry]*gmera1thirdstepdata[:globalemptyisometry])'
    gmera1thirdterminatedapproximation = gmera1firstapproximation + gmera1secondapproximation + (gmera1firststepdata[:globalcourierisometry]*gmera1secondstepdata[:globalcourierisometry]*gmera1thirdsteterminateddata[:globalemptyisometry])*(gmera1firststepdata[:globalcourierisometry]*gmera1secondstepdata[:globalcourierisometry]*gmera1thirdsteterminateddata[:globalemptyisometry])'
    fiosave(gmera1thirdterminatedapproximation, name="gmera1thirdterminatedapproximation")

    gmera1thirdalldiff = blockedcorrelations - gmera1thirdterminatedapproximation
    fiosave(gmera1thirdalldiff, name="gmera1thirdalldiff")

    @info("calculating L1 norm of the difference between the blocked correlations and the terminated approximate correlation at gmera step3")
    gmera1thirdalldiffL1norm = focktraceL1norm(gmera1thirdalldiff,systemsize^2*6)
    fiosave(gmera1thirdalldiffL1norm, name="gmera1thirdalldiffL1norm")


    couriercomposemapgmera1 = (gmera1firststepdata[:globalcourierisometry]*gmera1secondstepdata[:globalcourierisometry]*gmera1thirdstepdata[:globalcourierisometry])
    fiosave(couriercomposemapgmera1, name="gmera1couriercomposemapgmera")

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
    fiosave(gmera1firststepdata[:globalcourierisometry], name="gmera1firstglobalcourierisometry")
    fiosave(gmera1firststepdata[:globalemptyisometry], name="gmera1firstglobalemptyisometry")
    fiosave(gmera1firststepdata[:globalfilledisometry], name="gmera1firstglobalfilledisometry")
    fiosave(gmera1firststepdata[:minsvdcourier], name="gmera1firstminsvdforcourier")

    fiosave(gmera1secondstepdata[:localcorrelations], name="gmera1secondlocalcorrelations")
    fiosave(gmera1secondstepdata[:couriercorrelations], name="gmera1secondcouriercorrelations")
    fiosave(gmera1secondstepdata[:filledcorrelations], name="gmera1secondfilledcorrelations")
    fiosave(gmera1secondstepdata[:emptycorrelations], name="gmera1secondemptycorrelations")
    fiosave(gmera1secondstepdata[:courierH], name="gmera1secondcourierH")
    fiosave(gmera1secondstepdata[:filledH], name="gmera1secondfilledH")
    fiosave(gmera1secondstepdata[:emptyH], name="gmera1secondemptyH")
    fiosave(gmera1secondstepdata[:rawcouriercorrelations], name="gmera1secondrawcouriercorrelations")
    fiosave(gmera1secondstepdata[:globalcourierisometry], name="gmera1secondglobalcourierisometry")
    fiosave(gmera1secondstepdata[:globalemptyisometry], name="gmera1secondglobalemptyisometry")
    fiosave(gmera1secondstepdata[:globalfilledisometry], name="gmera1secondglobalfilledisometry")
    fiosave(gmera1secondstepdata[:minsvdcourier], name="gmera1secondminsvdforcourier")

    fiosave(gmera1thirdstepdata[:localcorrelations], name="gmera1thirdlocalcorrelations")
    fiosave(gmera1thirdstepdata[:couriercorrelations], name="gmera1thirdcouriercorrelations")
    fiosave(gmera1thirdstepdata[:filledcorrelations], name="gmera1thirdfilledcorrelations")
    fiosave(gmera1thirdstepdata[:emptycorrelations], name="gmera1thirdemptycorrelations")
    fiosave(gmera1thirdstepdata[:courierH], name="gmera1thirdcourierH")
    fiosave(gmera1thirdstepdata[:filledH], name="gmera1thirdfilledH")
    fiosave(gmera1thirdstepdata[:emptyH], name="gmera1thirdemptyH")
    fiosave(gmera1thirdstepdata[:rawcouriercorrelations], name="gmera1thirdrawcouriercorrelations")
    fiosave(gmera1thirdstepdata[:globalcourierisometry], name="gmera1thirdglobalcourierisometry")
    fiosave(gmera1thirdstepdata[:globalemptyisometry], name="gmera1thirdglobalemptyisometry")
    fiosave(gmera1thirdstepdata[:globalfilledisometry], name="gmera1thirdglobalfilledisometry")
    fiosave(gmera1thirdstepdata[:minsvdcourier], name="gmera1thirdminsvdforcourier")

    rg1correlations = gmera1thirdstepdata[:couriercorrelations]
    rg1H = gmera1thirdstepdata[:courierH]   

    return  rg1H,rg1correlations,couriercomposemapgmera1,gmera1approximatecorrelation,blockedcorrelations,noofflavourpermodeforlaterrg
end
export firstgmerastep

function intermediategmerastep(rgcorrelations,rgH,couriercomposemap,gmerprevsumapproximatecorrelation,rgstep,noofflavourpermode,blockedcorrelations,systemsize)
    rgblockedcorrelations,rgblockedH,rgblocker = blocking(rgcorrelations,rgH,2)
    prevrgstep = rgstep-1
    fiosave(rgblockedcorrelations, name="rg$prevrgstep"*"blockedcorrelations")
    fiosave(rgblockedH, name="rg$prevrgstep"*"blockedH")
    fiosave(rgblocker, name="rg$prevrgstep"*"blocker")
    couriercomposemap = couriercomposemap*rgblocker'

    # gmeraredefinedata = gmeraredefinestep(rgblockedcorrelations,rgblockedH,noofflavourpermode)
    gmerafirststepdata = gmerafirststep(rgblockedcorrelations,rgblockedH,18*noofflavourpermode,noofflavourpermode)
    gmerafirststepterminateddata = gmerafirststep(rgblockedcorrelations,rgblockedH,0,noofflavourpermode)

    gmerasecondstepdata = gmerasecondstep(gmerafirststepdata[:couriercorrelations],gmerafirststepdata[:courierH],12*noofflavourpermode,3*noofflavourpermode)
    gmerasecondstepterminateddata = gmerasecondstep(gmerafirststepdata[:couriercorrelations],gmerafirststepdata[:courierH],0,3*noofflavourpermode)

    gmerathirdstepdata = gmerathirdstep(gmerasecondstepdata[:couriercorrelations],gmerasecondstepdata[:courierH],6*noofflavourpermode,2*noofflavourpermode)
    gmerathirdstepterminateddata = gmerathirdstep(gmerasecondstepdata[:couriercorrelations],gmerasecondstepdata[:courierH],0,2*noofflavourpermode)

    gmerafirstapproximation = (couriercomposemap*gmerafirststepdata[:globalemptyisometry])*(couriercomposemap*gmerafirststepdata[:globalemptyisometry])'
    gmerafirstterminatedapproximation = gmerprevsumapproximatecorrelation+(couriercomposemap*gmerafirststepterminateddata[:globalemptyisometry])*(couriercomposemap*gmerafirststepterminateddata[:globalemptyisometry])'
    fiosave(gmerafirstterminatedapproximation, name="gmera$rgstep"*"firstterminatedapproximation")

    gmerafirstalldiff = blockedcorrelations - gmerafirstterminatedapproximation
    fiosave(gmerafirstalldiff, name="gmera$rgstep"*"firstalldiff")

    @info("calculating L1 norm of the difference between the blocked correlations and the terminated approximate correlation at gmera step1")
    gmerafirstalldiffL1norm = focktraceL1norm(gmerafirstalldiff,systemsize^2*6)
    fiosave(gmerafirstalldiffL1norm, name="gmera$rgstep"*"firstalldiffL1norm")

    gmerasecondapproximation = (couriercomposemap*gmerafirststepdata[:globalcourierisometry]*gmerasecondstepdata[:globalemptyisometry])*(couriercomposemap*gmerafirststepdata[:globalcourierisometry]*gmerasecondstepdata[:globalemptyisometry])'
    gmerasecondterminatedapproximation = gmerprevsumapproximatecorrelation+gmerasecondapproximation+(couriercomposemap*gmerafirststepdata[:globalcourierisometry]*gmerasecondstepterminateddata[:globalemptyisometry])*(couriercomposemap*gmerafirststepdata[:globalcourierisometry]*gmerasecondstepterminateddata[:globalemptyisometry])'
    fiosave(gmerasecondterminatedapproximation, name="gmera$rgstep"*"secondterminatedapproximation")

    gmerasecondalldiff = blockedcorrelations - gmerasecondterminatedapproximation
    fiosave(gmerasecondalldiff, name="gmera$rgstep"*"secondalldiff")

    @info("calculating L1 norm of the difference between the blocked correlations and the terminated approximate correlation at gmera step2")
    gmerasecondalldiffL1norm = focktraceL1norm(gmerasecondalldiff,systemsize^2*6)
    fiosave(gmerasecondalldiffL1norm, name="gmera$rgstep"*"secondalldiffL1norm")

    gmerathirdapproximation = (couriercomposemap*gmerafirststepdata[:globalcourierisometry]*gmerasecondstepdata[:globalcourierisometry]*gmerathirdstepdata[:globalemptyisometry])*(couriercomposemap*gmerafirststepdata[:globalcourierisometry]*gmerasecondstepdata[:globalcourierisometry]*gmerathirdstepdata[:globalemptyisometry])'
    gmerathirdterminatedapproximation = gmerprevsumapproximatecorrelation+gmerasecondapproximation+gmerathirdapproximation+(couriercomposemap*gmerafirststepdata[:globalcourierisometry]*gmerasecondstepdata[:globalcourierisometry]*gmerathirdstepterminateddata[:globalemptyisometry])*(couriercomposemap*gmerafirststepdata[:globalcourierisometry]*gmerasecondstepdata[:globalcourierisometry]*gmerathirdstepterminateddata[:globalemptyisometry])'
    fiosave(gmerathirdterminatedapproximation, name="gmera$rgstep"*"thirdterminatedapproximation")

    gmerathirdalldiff = blockedcorrelations - gmerasecondterminatedapproximation
    fiosave(gmerathirdalldiff, name="gmera$rgstep"*"thirdalldiff")

    @info("calculating L1 norm of the difference between the blocked correlations and the terminated approximate correlation at gmera step2")
    gmerathirdalldiffL1norm = focktraceL1norm(gmerathirdalldiff,systemsize^2*6)
    fiosave(gmerathirdalldiffL1norm, name="gmera$rgstep"*"thirdalldiffL1norm")
    
    couriercomposemapgmera = couriercomposemap*(gmerafirststepdata[:globalcourierisometry]*gmerasecondstepdata[:globalcourierisometry]*gmerathirdstepdata[:globalcourierisometry])
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
    fiosave(gmerafirststepdata[:globalcourierisometry], name="gmera$rgstep"*"firstglobalcourierisometry")
    fiosave(gmerafirststepdata[:globalemptyisometry], name="gmera$rgstep"*"firstglobalemptyisometry")
    fiosave(gmerafirststepdata[:globalfilledisometry], name="gmera$rgstep"*"firstglobalfilledisometry")
    fiosave(gmerafirststepdata[:minsvdcourier], name="gmera$rgstep"*"firstminsvdforcourier")

    fiosave(gmerasecondstepdata[:localcorrelations], name="gmera$rgstep"*"secondlocalcorrelations")
    fiosave(gmerasecondstepdata[:couriercorrelations], name="gmera$rgstep"*"secondcouriercorrelations")
    fiosave(gmerasecondstepdata[:filledcorrelations], name="gmera$rgstep"*"secondfilledcorrelations")
    fiosave(gmerasecondstepdata[:emptycorrelations], name="gmera$rgstep"*"secondemptycorrelations")
    fiosave(gmerasecondstepdata[:courierH], name="gmera$rgstep"*"secondcourierH")
    fiosave(gmerasecondstepdata[:filledH], name="gmera$rgstep"*"secondfilledH")
    fiosave(gmerasecondstepdata[:emptyH], name="gmera$rgstep"*"secondemptyH")
    fiosave(gmerasecondstepdata[:rawcouriercorrelations], name="gmera$rgstep"*"secondrawcouriercorrelations")
    fiosave(gmerasecondstepdata[:globalcourierisometry], name="gmera$rgstep"*"secondglobalcourierisometry")
    fiosave(gmerasecondstepdata[:globalemptyisometry], name="gmera$rgstep"*"secondglobalemptyisometry")
    fiosave(gmerasecondstepdata[:globalfilledisometry], name="gmera$rgstep"*"secondglobalfilledisometry")
    fiosave(gmerasecondstepdata[:minsvdcourier], name="gmera$rgstep"*"secondminsvdforcourier")

    fiosave(gmerathirdstepdata[:localcorrelations], name="gmera$rgstep"*"thirdlocalcorrelations")
    fiosave(gmerathirdstepdata[:couriercorrelations], name="gmera$rgstep"*"thirdcouriercorrelations")
    fiosave(gmerathirdstepdata[:filledcorrelations], name="gmera$rgstep"*"thirdfilledcorrelations")
    fiosave(gmerathirdstepdata[:emptycorrelations], name="gmera$rgstep"*"thirdemptycorrelations")
    fiosave(gmerathirdstepdata[:courierH], name="gmera$rgstep"*"thirdcourierH")
    fiosave(gmerathirdstepdata[:filledH], name="gmera$rgstep"*"thirdfilledH")
    fiosave(gmerathirdstepdata[:emptyH], name="gmera$rgstep"*"thirdemptyH")
    fiosave(gmerathirdstepdata[:rawcouriercorrelations], name="gmera$rgstep"*"thirdrawcouriercorrelations")
    fiosave(gmerathirdstepdata[:globalcourierisometry], name="gmera$rgstep"*"thirdglobalcourierisometry")
    fiosave(gmerathirdstepdata[:globalemptyisometry], name="gmera$rgstep"*"thirdglobalemptyisometry")
    fiosave(gmerathirdstepdata[:globalfilledisometry], name="gmera$rgstep"*"thirdglobalfilledisometry")
    fiosave(gmerathirdstepdata[:minsvdcourier], name="gmera$rgstep"*"thirdminsvdforcourier")

    rgcorrelations = gmerathirdstepdata[:couriercorrelations]
    rgH = gmerathirdstepdata[:courierH]   
    gmerasumapproximatecorrelationsofar = gmerprevsumapproximatecorrelation+gmeraapproximatecorrelation
    fiosave(gmerasumapproximatecorrelationsofar, name="gmera$rgsteps"*"approximatecorrelationsofar")
    return  rgH,rgcorrelations,couriercomposemapgmera,gmerasumapproximatecorrelationsofar
end
export intermediategmerastep

function finalgmerastep(rgcorrelations,rgH,couriercomposemap,gmerprevsumapproximatecorrelation,blockedcorrelations,rgstep,systemsize,noofflavourpermode)
    remainsize = (rgcorrelations|>getinspace|>getcrystal|>size)[1]
    if remainsize == 3
        @info("blocking with scaling 3 cause only a system of 3 by 3 regions is left")
        rgblockedcorrelations,rgblockedH,rgblocker = blocking(rgcorrelations,rgH,3)
    elseif remainsize == 5
        @info("blocking with scaling 5 cause only a system of 5 by 5 regions is left")
        rgblockedcorrelations,rgblockedH,rgblocker = blocking(rgcorrelations,rgH,5)
    else
        @info("blocking with scaling 2")
        rgblockedcorrelations,rgblockedH,rgblocker = blocking(rgcorrelations,rgH,2)
    end
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
end
export finalgmerastep


# function firstgmerastepwifredef(correlations,H,scaling)
#     blockedcorrelations,blockedH,blocker = blocking(correlations,H,scaling)
#     fiosave(blockedcorrelations, name="blockedcorrelations")
#     fiosave(blockedH, name="blockedH")
#     fiosave(blocker, name="blocker")

#     noofmodesinlocalreg = (6*(scaling^2))|>Int
#     noofdistillablemodes = ((1/4)*noofmodesinlocalreg)|>Int
#     noofcouriermodesinfirststep = noofmodesinlocalreg - noofdistillablemodes
#     # println(noofcouriermodesinfirststep)
#     noofcouriermodesinsecondstep = noofcouriermodesinfirststep - noofdistillablemodes
#     noofcouriermodesinthirdstep = noofcouriermodesinsecondstep - noofdistillablemodes

#     gmera1redefinedata = gmeraredefinestep(blockedcorrelations,blockedH,1)
#     gmera1firststepdata = gmerafirststep(gmera1redefinedata[:rotcorrelations],gmera1redefinedata[:rotH],noofcouriermodesinfirststep,div(noofmodesinlocalreg,6))
#     gmera1secondstepdata = gmerasecondstep(gmera1firststepdata[:couriercorrelations],gmera1firststepdata[:courierH],noofcouriermodesinsecondstep,div(noofcouriermodesinfirststep,6))
#     gmera1thirdstepdata = gmerathirdstep(gmera1secondstepdata[:couriercorrelations],gmera1secondstepdata[:courierH],noofcouriermodesinthirdstep,div(noofcouriermodesinsecondstep,6))

#     noofflavourpermodeforlaterrg = (noofcouriermodesinthirdstep/6)|>Int

#     gmera1firstapproximation = (gmera1redefinedata[:globalunitary]*gmera1firststepdata[:globalemptyisometry])*(gmera1redefinedata[:globalunitary]*gmera1firststepdata[:globalemptyisometry])'
#     gmera1secondapproximation = (gmera1redefinedata[:globalunitary]*gmera1firststepdata[:globalcourierisometry]*gmera1secondstepdata[:globalemptyisometry])*(gmera1redefinedata[:globalunitary]*gmera1firststepdata[:globalcourierisometry]*gmera1secondstepdata[:globalemptyisometry])'
#     gmera1thirdapproximation = (gmera1redefinedata[:globalunitary]*gmera1firststepdata[:globalcourierisometry]*gmera1secondstepdata[:globalcourierisometry]*gmera1thirdstepdata[:globalemptyisometry])*(gmera1redefinedata[:globalunitary]*gmera1firststepdata[:globalcourierisometry]*gmera1secondstepdata[:globalcourierisometry]*gmera1thirdstepdata[:globalemptyisometry])'

#     couriercomposemapgmera1 = (gmera1redefinedata[:globalunitary]*gmera1firststepdata[:globalcourierisometry]*gmera1secondstepdata[:globalcourierisometry]*gmera1thirdstepdata[:globalcourierisometry])

#     gmera1approximatecorrelation = gmera1firstapproximation + gmera1secondapproximation + gmera1thirdapproximation
#     fiosave(gmera1approximatecorrelation, name="gmera1approximatecorrelation")

#     fiosave(gmera1redefinedata[:localcorrelations], name="gmera1redeflocalcorrelationss")
#     fiosave(gmera1redefinedata[:globalunitary], name="gmera1redefglobalunitary")
#     fiosave(gmera1redefinedata[:rotcorrelations], name="gmera1rotcorrelations")
#     fiosave(gmera1redefinedata[:rotH], name="gmera1rotH")
#     fiosave(gmera1redefinedata[:minsvdredef], name="gmera1minsvdredef")

#     fiosave(gmera1firststepdata[:localcorrelations], name="gmera1firstlocalcorrelations")
#     fiosave(gmera1firststepdata[:couriercorrelations], name="gmera1firstcouriercorrelations")
#     fiosave(gmera1firststepdata[:filledcorrelations], name="gmera1firstfilledcorrelations")
#     fiosave(gmera1firststepdata[:emptycorrelations], name="gmera1firstemptycorrelations")
#     fiosave(gmera1firststepdata[:courierH], name="gmera1firstcourierH")
#     fiosave(gmera1firststepdata[:filledH], name="gmera1firstfilledH")
#     fiosave(gmera1firststepdata[:emptyH], name="gmera1firstemptyH")
#     fiosave(gmera1firststepdata[:rawcouriercorrelations], name="gmera1firstrawcouriercorrelations")
#     fiosave(gmera1firststepdata[:globalcourierisometry], name="gmera1firstglobalcourierisometry")
#     fiosave(gmera1firststepdata[:globalemptyisometry], name="gmera1firstglobalemptyisometry")
#     fiosave(gmera1firststepdata[:globalfilledisometry], name="gmera1firstglobalfilledisometry")
#     fiosave(gmera1firststepdata[:minsvdcourier], name="gmera1firstminsvdforcourier")

#     fiosave(gmera1secondstepdata[:localcorrelations], name="gmera1secondlocalcorrelations")
#     fiosave(gmera1secondstepdata[:couriercorrelations], name="gmera1secondcouriercorrelations")
#     fiosave(gmera1secondstepdata[:filledcorrelations], name="gmera1secondfilledcorrelations")
#     fiosave(gmera1secondstepdata[:emptycorrelations], name="gmera1secondemptycorrelations")
#     fiosave(gmera1secondstepdata[:courierH], name="gmera1secondcourierH")
#     fiosave(gmera1secondstepdata[:filledH], name="gmera1secondfilledH")
#     fiosave(gmera1secondstepdata[:emptyH], name="gmera1secondemptyH")
#     fiosave(gmera1secondstepdata[:rawcouriercorrelations], name="gmera1secondrawcouriercorrelations")
#     fiosave(gmera1secondstepdata[:globalcourierisometry], name="gmera1secondglobalcourierisometry")
#     fiosave(gmera1secondstepdata[:globalemptyisometry], name="gmera1secondglobalemptyisometry")
#     fiosave(gmera1secondstepdata[:globalfilledisometry], name="gmera1secondglobalfilledisometry")
#     fiosave(gmera1secondstepdata[:minsvdcourier], name="gmera1secondminsvdforcourier")

#     fiosave(gmera1thirdstepdata[:localcorrelations], name="gmera1thirdlocalcorrelations")
#     fiosave(gmera1thirdstepdata[:couriercorrelations], name="gmera1thirdcouriercorrelations")
#     fiosave(gmera1thirdstepdata[:filledcorrelations], name="gmera1thirdfilledcorrelations")
#     fiosave(gmera1thirdstepdata[:emptycorrelations], name="gmera1thirdemptycorrelations")
#     fiosave(gmera1thirdstepdata[:courierH], name="gmera1thirdcourierH")
#     fiosave(gmera1thirdstepdata[:filledH], name="gmera1thirdfilledH")
#     fiosave(gmera1thirdstepdata[:emptyH], name="gmera1thirdemptyH")
#     fiosave(gmera1thirdstepdata[:rawcouriercorrelations], name="gmera1thirdrawcouriercorrelations")
#     fiosave(gmera1thirdstepdata[:globalcourierisometry], name="gmera1thirdglobalcourierisometry")
#     fiosave(gmera1thirdstepdata[:globalemptyisometry], name="gmera1thirdglobalemptyisometry")
#     fiosave(gmera1thirdstepdata[:globalfilledisometry], name="gmera1thirdglobalfilledisometry")
#     fiosave(gmera1thirdstepdata[:minsvdcourier], name="gmera1thirdminsvdforcourier")

#     rg1correlations = gmera1thirdstepdata[:couriercorrelations]
#     rg1H = gmera1thirdstepdata[:courierH]   

#     return  rg1H,rg1correlations,couriercomposemapgmera1,gmera1approximatecorrelation,blockedcorrelations,noofflavourpermodeforlaterrg
# end
# export firstgmerastepwifredef

# function intermediategmerastepwifredef(rgcorrelations,rgH,couriercomposemap,gmerprevsumapproximatecorrelation,rgstep,noofflavourpermode)
#     size = (rgcorrelations|>getinspace|>getcrystal|>size)[1]
#     if size == 3
#         @info("blocking with scaling 3 cause only a system of 3 by 3 regions is left")
#         rgblockedcorrelations,rgblockedH,rgblocker = blocking(rgcorrelations,rgH,3)
#     else
#         @info("blocking with scaling 2")
#         rgblockedcorrelations,rgblockedH,rgblocker = blocking(rgcorrelations,rgH,2)
#     end
#     prevrgstep = rgstep-1
#     # fiosave(rgblockedcorrelations, name="rg$prevrgstep"*"blockedcorrelations")
#     # fiosave(rgblockedH, name="rg$prevrgstep"*"blockedH")
#     # fiosave(rgblocker, name="rg$prevrgstep"*"blocker")
#     couriercomposemap = couriercomposemap*rgblocker'

#     gmeraredefinedata = gmeraredefinestep(rgblockedcorrelations,rgblockedH,noofflavourpermode)
#     gmerafirststepdata = gmerafirststep(gmeraredefinedata[:rotcorrelations],gmeraredefinedata[:rotH],18*noofflavourpermode,4*noofflavourpermode)
#     gmerasecondstepdata = gmerasecondstep(gmerafirststepdata[:couriercorrelations],gmerafirststepdata[:courierH],12*noofflavourpermode,3*noofflavourpermode)
#     gmerathirdstepdata = gmerathirdstep(gmerasecondstepdata[:couriercorrelations],gmerasecondstepdata[:courierH],6*noofflavourpermode,2*noofflavourpermode)

#     gmerafirstapproximation = (couriercomposemap*gmeraredefinedata[:globalunitary]*gmerafirststepdata[:globalemptyisometry])*(couriercomposemap*gmeraredefinedata[:globalunitary]*gmerafirststepdata[:globalemptyisometry])'
#     gmerasecondapproximation = (couriercomposemap*gmeraredefinedata[:globalunitary]*gmerafirststepdata[:globalcourierisometry]*gmerasecondstepdata[:globalemptyisometry])*(couriercomposemap*gmeraredefinedata[:globalunitary]*gmerafirststepdata[:globalcourierisometry]*gmerasecondstepdata[:globalemptyisometry])'
#     gmerathirdapproximation = (couriercomposemap*gmeraredefinedata[:globalunitary]*gmerafirststepdata[:globalcourierisometry]*gmerasecondstepdata[:globalcourierisometry]*gmerathirdstepdata[:globalemptyisometry])*(couriercomposemap*gmeraredefinedata[:globalunitary]*gmerafirststepdata[:globalcourierisometry]*gmerasecondstepdata[:globalcourierisometry]*gmerathirdstepdata[:globalemptyisometry])'
    
#     couriercomposemapgmera = couriercomposemap*(gmeraredefinedata[:globalunitary]*gmerafirststepdata[:globalcourierisometry]*gmerasecondstepdata[:globalcourierisometry]*gmerathirdstepdata[:globalcourierisometry])
#     rgsteps = prod([string(i) for i in rgstep])
#     # fiosave(couriercomposemapgmera, name="gmera$rgsteps"*"couriercomposemapgmera")

#     gmeraapproximatecorrelation = gmerafirstapproximation + gmerasecondapproximation  + gmerathirdapproximation
#     # fiosave(gmeraapproximatecorrelation, name="gmera$rgstep"*"approximatecorrelation")

#     fiosave(gmeraredefinedata[:localcorrelations], name="gmera$rgstep"*"redeflocalcorrelationss")
#     fiosave(gmeraredefinedata[:globalunitary], name="gmera$rgstep"*"redefglobalunitary")
#     fiosave(gmeraredefinedata[:rotcorrelations], name="gmera$rgstep"*"rotcorrelations")
#     fiosave(gmeraredefinedata[:rotH], name="gmera$rgstep"*"rotH")
#     fiosave(gmeraredefinedata[:minsvdredef], name="gmera$rgstep"*"minsvdredef")

#     fiosave(gmerafirststepdata[:localcorrelations], name="gmera$rgstep"*"firstlocalcorrelations")
#     fiosave(gmerafirststepdata[:couriercorrelations], name="gmera$rgstep"*"firstcouriercorrelations")
#     fiosave(gmerafirststepdata[:filledcorrelations], name="gmera$rgstep"*"firstfilledcorrelations")
#     fiosave(gmerafirststepdata[:emptycorrelations], name="gmera$rgstep"*"firstemptycorrelations")
#     fiosave(gmerafirststepdata[:courierH], name="gmera$rgstep"*"firstcourierH")
#     fiosave(gmerafirststepdata[:filledH], name="gmera$rgstep"*"firstfilledH")
#     fiosave(gmerafirststepdata[:emptyH], name="gmera$rgstep"*"firstemptyH")
#     fiosave(gmerafirststepdata[:rawcouriercorrelations], name="gmera$rgstep"*"firstrawcouriercorrelations")
#     fiosave(gmerafirststepdata[:globalcourierisometry], name="gmera$rgstep"*"firstglobalcourierisometry")
#     fiosave(gmerafirststepdata[:globalemptyisometry], name="gmera$rgstep"*"firstglobalemptyisometry")
#     fiosave(gmerafirststepdata[:globalfilledisometry], name="gmera$rgstep"*"firstglobalfilledisometry")
#     fiosave(gmerafirststepdata[:minsvdcourier], name="gmera$rgstep"*"firstminsvdforcourier")

#     fiosave(gmerasecondstepdata[:localcorrelations], name="gmera$rgstep"*"secondlocalcorrelations")
#     fiosave(gmerasecondstepdata[:couriercorrelations], name="gmera$rgstep"*"secondcouriercorrelations")
#     fiosave(gmerasecondstepdata[:filledcorrelations], name="gmera$rgstep"*"secondfilledcorrelations")
#     fiosave(gmerasecondstepdata[:emptycorrelations], name="gmera$rgstep"*"secondemptycorrelations")
#     fiosave(gmerasecondstepdata[:courierH], name="gmera$rgstep"*"secondcourierH")
#     fiosave(gmerasecondstepdata[:filledH], name="gmera$rgstep"*"secondfilledH")
#     fiosave(gmerasecondstepdata[:emptyH], name="gmera$rgstep"*"secondemptyH")
#     fiosave(gmerasecondstepdata[:rawcouriercorrelations], name="gmera$rgstep"*"secondrawcouriercorrelations")
#     fiosave(gmerasecondstepdata[:globalcourierisometry], name="gmera$rgstep"*"secondglobalcourierisometry")
#     fiosave(gmerasecondstepdata[:globalemptyisometry], name="gmera$rgstep"*"secondglobalemptyisometry")
#     fiosave(gmerasecondstepdata[:globalfilledisometry], name="gmera$rgstep"*"secondglobalfilledisometry")
#     fiosave(gmerasecondstepdata[:minsvdcourier], name="gmera$rgstep"*"secondminsvdforcourier")

#     fiosave(gmerathirdstepdata[:localcorrelations], name="gmera$rgstep"*"thirdlocalcorrelations")
#     fiosave(gmerathirdstepdata[:couriercorrelations], name="gmera$rgstep"*"thirdcouriercorrelations")
#     fiosave(gmerathirdstepdata[:filledcorrelations], name="gmera$rgstep"*"thirdfilledcorrelations")
#     fiosave(gmerathirdstepdata[:emptycorrelations], name="gmera$rgstep"*"thirdemptycorrelations")
#     fiosave(gmerathirdstepdata[:courierH], name="gmera$rgstep"*"thirdcourierH")
#     fiosave(gmerathirdstepdata[:filledH], name="gmera$rgstep"*"thirdfilledH")
#     fiosave(gmerathirdstepdata[:emptyH], name="gmera$rgstep"*"thirdemptyH")
#     fiosave(gmerathirdstepdata[:rawcouriercorrelations], name="gmera$rgstep"*"thirdrawcouriercorrelations")
#     fiosave(gmerathirdstepdata[:globalcourierisometry], name="gmera$rgstep"*"thirdglobalcourierisometry")
#     fiosave(gmerathirdstepdata[:globalemptyisometry], name="gmera$rgstep"*"thirdglobalemptyisometry")
#     fiosave(gmerathirdstepdata[:globalfilledisometry], name="gmera$rgstep"*"thirdglobalfilledisometry")
#     fiosave(gmerathirdstepdata[:minsvdcourier], name="gmera$rgstep"*"thirdminsvdforcourier")

#     rgcorrelations = gmerathirdstepdata[:couriercorrelations]
#     rgH = gmerathirdstepdata[:courierH]   
#     gmerasumapproximatecorrelationsofar = gmerprevsumapproximatecorrelation+gmeraapproximatecorrelation
#     fiosave(gmerasumapproximatecorrelationsofar, name="gmera$rgsteps"*"approximatecorrelationsofar")
#     return  rgH,rgcorrelations,couriercomposemapgmera,gmerasumapproximatecorrelationsofar
# end
# export intermediategmerastepwifredef

# function finalgmerastepwifredef(rgcorrelations,rgH,couriercomposemap,gmerprevsumapproximatecorrelation,blockedcorrelations,rgstep,systemsize,noofflavourpermode)
#     rgblockedcorrelations,rgblockedH,rgblocker = blocking(rgcorrelations,rgH,2)
#     prevrgstep = rgstep-1
#     fiosave(rgblockedcorrelations, name="rg$prevrgstep"*"blockedcorrelations")
#     fiosave(rgblockedH, name="rg$prevrgstep"*"blockedH")
#     fiosave(rgblocker, name="rg$prevrgstep"*"blocker")
#     couriercomposemap = couriercomposemap*rgblocker'

#     gmerafinalstepdata = gmerafinalstep(rgblockedcorrelations,rgblockedH,noofflavourpermode)
#     fiosave(gmerafinalstepdata[:localcorrelations], name="gmera$rgstep"*"finallocalcorrelations")
#     fiosave(gmerafinalstepdata[:filledcorrelations], name="gmera$rgstep"*"finalfilledcorrelations")
#     fiosave(gmerafinalstepdata[:emptycorrelations], name="gmera$rgstep"*"finalemptycorrelations")
#     fiosave(gmerafinalstepdata[:filledH], name="gmera$rgstep"*"finalfilledH")
#     fiosave(gmerafinalstepdata[:emptyH], name="gmera$rgstep"*"finalemptyH")
#     fiosave(gmerafinalstepdata[:globalemptyisometry], name="gmera$rgstep"*"finalglobalemptyisometry")
#     fiosave(gmerafinalstepdata[:globalfilledisometry], name="gmera$rgstep"*"finalglobalfilledisometry")

#     gmerafinalapproximation = (couriercomposemap*gmerafinalstepdata[:globalemptyisometry])*(couriercomposemap*gmerafinalstepdata[:globalemptyisometry])'
#     fiosave(gmerafinalapproximation, name="gmera$rgstep"*"finalapproximation")
#     gmeraallapproximatecorrelation = gmerprevsumapproximatecorrelation+gmerafinalapproximation
#     fiosave(gmeraallapproximatecorrelation, name="gmeraallapproximatecorrelation")

#     gmeraalldiff = blockedcorrelations - gmeraallapproximatecorrelation
#     fiosave(gmeraalldiff, name="gmeraalldiff")

#     @info("calculating L1 norm of the difference between the blocked correlations and the approximate correlation")
#     gmeraalldiffL1norm = vectorizeL1norm(gmeraalldiff)/systemsize^2*6
#     fiosave(gmeraalldiffL1norm, name="gmeraalldiffL1norm")
# end
# export finalgmerastepwifredef