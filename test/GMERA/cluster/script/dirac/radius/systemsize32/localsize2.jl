using Zipper
using LinearAlgebra,Plots

setmaxthreads(8)

plotlyjs()
systemsize=32
correlations,H = generatesystem(-0.1+0im, 0.1+0im,-1 + 0im,0im,systemsize)

function blocking(correlations,H,scale)
    crystalfock = correlations|>getoutspace
    scale = Scale([scale 0; 0 scale], crystalfock|>getcrystal|>getspace)
    @info("Performing rgblocking...",scale)
    @info("Generating rgblocking transformation...")
    blocker = @time scale * crystalfock
    @info("Performing rgblocking on correlations...")
    blockedcorrelations = @time blocker * correlations * blocker'
    blockedH = @time blocker * H * blocker'
    return blockedcorrelations,blockedH
end

function gmerafirststepwifentanglementcontour(blockedcorrelations,blockedH,nooffrozenmodes)
    blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal
    blockedcrystalfock = blockedcorrelations|>getoutspace
    blockedspace::RealSpace = blockedcrystal|>getspace

    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    firstcenter = [0,0] ∈ blockedspace
    firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=firstcenter, metricspace=blockedspace)
    firsthexagonalregionfock = quantize(firsthexagonalregion,1)

    localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
    localspectrum = localcorrelations|>eigspec

    frozenseedsoffset = rankedandgroupoffsets(firsthexagonalregion,1)
    courierseedsoffset = firsthexagonalregion-frozenseedsoffset
    courierseedsfock = quantize(courierseedsoffset,1)

    display(visualize(courierseedsoffset))

    nooffilledmodes = div(length(frozenseedsoffset),2)
    noofemptymodes = nooffilledmodes
    noofcouriermodes = length(courierseedsoffset)

    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes,noofcouriermodes,noofemptymodes])
    localcourierisometry = localstates[2]|>FockMap
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[3]|>FockMap

    courierseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in courierseedsfock)]
    localwanniercourier,minsvd,traceofsv=localwannierization(localcourierisometry, courierseeds)

    localwannierresults =  Dict(:localwanniercourier => localwanniercourier, :localfilledisometry => localfilledisometry,
                                :localemptyisometry => localemptyisometry, :courierseeds => courierseeds)

    wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

    shiftedfirstcenter1 = [1,0] ∈ blockedspace
    shiftedfirstcenter2 = [0,1] ∈ blockedspace
    shiftedfirstcenter3 = [1,1] ∈ blockedspace

    firstshiftedhexagonalcenters = [shiftedfirstcenter1, shiftedfirstcenter2, shiftedfirstcenter3]

    for hexagonalcenter in firstshiftedhexagonalcenters
        shiftedfirsthexagonalregion = firsthexagonalregion.+hexagonalcenter
        shiftedfirsthexagonalregionfock = quantize(shiftedfirsthexagonalregion,1)
        
        localcorrelations = regioncorrelations(blockedcorrelations,shiftedfirsthexagonalregionfock)
        localspectrum = localcorrelations|>eigspec
        
        frozenseedsoffset = rankedandgroupoffsets(shiftedfirsthexagonalregion,1)
        courierseedsoffset = shiftedfirsthexagonalregion-frozenseedsoffset
        courierseedsfock = quantize(courierseedsoffset,1)

        nooffilledmodes = div(length(frozenseedsoffset),2)
        noofemptymodes = nooffilledmodes
        noofcouriermodes = length(courierseedsoffset)
    
        localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes,noofcouriermodes,noofemptymodes])
        localcourierisometry = localstates[2]|>FockMap
        localfilledisometry = localstates[1]|>FockMap
        localemptyisometry = localstates[3]|>FockMap
    
    
        courierseeds = idmap(shiftedfirsthexagonalregionfock, shiftedfirsthexagonalregionfock)[:,FockSpace(mode for mode in courierseedsfock)]
    
        localwanniercourier,minsvd,traceofsv = localwannierization(localcourierisometry, courierseeds)
    
    
        shiftedlocalwannierresults =  Dict(:localwanniercourier => localwanniercourier, :localfilledisometry => localfilledisometry,
                                            :localemptyisometry => localemptyisometry, :courierseeds => courierseeds)
        wannierinfos[shiftedfirsthexagonalregionfock] =  shiftedlocalwannierresults
    end

    ref = [quantize(firsthexagonalregion.+hexagonalcenter,1) for hexagonalcenter in firstshiftedhexagonalcenters]
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

    @info "Computing local courier states..."
    leftrestrict = fourier(wanniercourierisometry|>getoutspace, firsthexagonalregionfock) / (blockedcrystal|>vol|>sqrt)
    rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
    wanniercourierstate = leftrestrict' * wanniercourierisometry * rightrestrict

    couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
    couriercorrelations|>crystalspectrum|>visualize
    courierH = wanniercourierisometry' * blockedH * wanniercourierisometry

    filledcorrelations = globalfilledisometry' * blockedcorrelations * globalfilledisometry
    filledH = globalfilledisometry' * blockedH * globalfilledisometry

    emptycorrelations = globalemptyisometry' * blockedcorrelations * globalemptyisometry
    emptyH = globalemptyisometry' * blockedH * globalemptyisometry

    couriercorrelationspectrum = couriercorrelations |> crystalspectrum

    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap
    return Dict(
        :localcorrelations=>localcorrelations,
        :couriercorrelations=>purifiedcouriercorrelations,
        :courierH=>courierH,
        :filledH=>filledH,
        :emptyH=>emptyH,
        :filledcorrelations=>filledcorrelations,
        :emptycorrelations=>emptycorrelations,
        :wanniercourierisometry=>wanniercourierisometry,
        :globalfilledisometry=>globalfilledisometry,
        :globalemptyisometry=>globalemptyisometry,
        :wanniercourierstates=>wanniercourierstate|>RegionState,
        :minsvdforcourier=>minsvd,
        :traceofsv=>traceofsv,
        :rawcouriercorrelations=>couriercorrelations)
end

function gmerasecondstepwifentanglementcontour(firstcouriercorrelations,firstcourierH,nooffrozenmodes)
    firstgmeracrystalfock = firstcouriercorrelations|>getoutspace
    firstgmeracrystal::Crystal = firstgmeracrystalfock|>getcrystal
    firstgmeraspace::RealSpace = firstgmeracrystal|>getspace

    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    secondcenter = [2/3,-1/3] ∈ firstgmeraspace
    secondhexagonalregion = gethexagonalregion(rot=refrot,crystal=firstgmeracrystal, center=secondcenter, metricspace=firstgmeraspace)
    secondhexagonalregionfock = quantize(secondhexagonalregion,1)

    localcorrelations = regioncorrelations(firstcouriercorrelations,secondhexagonalregionfock)
    localspectrum = localcorrelations|>eigspech

    
    frozenseedsoffset = rankedandgroupoffsets(secondhexagonalregion,1)
    courierseedsoffset = secondhexagonalregion-frozenseedsoffset
    courierseedsfock = quantize(courierseedsoffset,1)

    display(visualize(courierseedsoffset))

    nooffilledmodes = div(length(frozenseedsoffset),2)
    noofemptymodes = nooffilledmodes
    noofcouriermodes = length(courierseedsoffset)

    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes,noofcouriermodes,noofemptymodes])
    localcourierisometry = localstates[2]|>FockMap
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[3]|>FockMap

    courierseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in courierseedsfock)]
    localwanniercourier,minsvd,traceofsv=localwannierization(localcourierisometry, courierseeds)

    localwannierresults =  Dict(:localwanniercourier => localwanniercourier, :localfilledisometry => localfilledisometry,
                                :localemptyisometry => localemptyisometry, :courierseeds => courierseeds)

    wannierinfos =  Dict(secondhexagonalregionfock=>localwannierresults)

    shiftedsecondcenter1 = [0,1] ∈ firstgmeraspace
    shiftedsecondcenter2 = [-1,1] ∈ firstgmeraspace

    secondshiftedhexagonalcenterlist = [shiftedsecondcenter1, shiftedsecondcenter2]

    for hexagonalcenter in secondshiftedhexagonalcenterlist
        shiftedsecondhexagonalregion = secondhexagonalregion.+hexagonalcenter
        shiftedsecondhexagonalregionfock = quantize(shiftedsecondhexagonalregion,1)
        
        localcorrelations = regioncorrelations(firstcouriercorrelations,shiftedsecondhexagonalregionfock)
        localspectrum = localcorrelations|>eigspech

        frozenseedsoffset = rankedandgroupoffsets(shiftedsecondhexagonalregion,1)
        courierseedsoffset = shiftedsecondhexagonalregion-frozenseedsoffset
        courierseedsfock = quantize(courierseedsoffset,1)

        nooffilledmodes = div(length(frozenseedsoffset),2)
        noofemptymodes = nooffilledmodes
        noofcouriermodes = length(courierseedsoffset)

        localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes,noofcouriermodes,noofemptymodes])
        localcourierisometry = localstates[2]|>FockMap
        localfilledisometry = localstates[1]|>FockMap
        localemptyisometry = localstates[3]|>FockMap
        
        courierseeds = idmap(shiftedsecondhexagonalregionfock, shiftedsecondhexagonalregionfock)[:,FockSpace(mode for mode in courierseedsfock)]
        localwanniercourier,minsvd,traceofsv=localwannierization(localcourierisometry, courierseeds)
    
    
        shiftedlocalwannierresults = Dict(:localwanniercourier => localwanniercourier, :localfilledisometry => localfilledisometry,
                                    :localemptyisometry => localemptyisometry, :courierseeds => courierseeds)
        wannierinfos[shiftedsecondhexagonalregionfock] =  shiftedlocalwannierresults
    end

    ref = [quantize(secondhexagonalregion.+hexagonalcenter,1) for hexagonalcenter in secondshiftedhexagonalcenterlist]
    secondhexagonalregionfocklist = [secondhexagonalregionfock,ref...]

    extendedwannierizedcourier =  sum(wannierinfos[regionfock][:localwanniercourier] for regionfock in secondhexagonalregionfocklist)
    extendedwannierizedfilled =  sum(wannierinfos[regionfock][:localfilledisometry] for regionfock in secondhexagonalregionfocklist)
    extendedwannierizedempty =  sum(wannierinfos[regionfock][:localemptyisometry] for regionfock in secondhexagonalregionfocklist)
    
                
    origin = [0, 0] ∈ firstgmeraspace
    refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

    wanniercourierisometry = globalwannierfunction(firstcouriercorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
    globalfilledisometry = globalwannierfunction(firstcouriercorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
    globalemptyisometry = globalwannierfunction(firstcouriercorrelations,extendedwannierizedempty[:,refunictcellfockempty])

    @info "Computing local courier states..."
    leftrestrict = fourier(wanniercourierisometry|>getoutspace, secondhexagonalregionfock) / (firstgmeracrystal|>vol|>sqrt)
    rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
    wanniercourierstate = leftrestrict' * wanniercourierisometry * rightrestrict

    couriercorrelations = wanniercourierisometry' * firstcouriercorrelations * wanniercourierisometry
    courierH = wanniercourierisometry' * firstcourierH * wanniercourierisometry

    filledcorrelations = globalfilledisometry' * firstcouriercorrelations * globalfilledisometry
    filledH = globalfilledisometry' * firstcourierH * globalfilledisometry

    emptycorrelations = globalemptyisometry' * firstcouriercorrelations * globalemptyisometry
    emptyH = globalemptyisometry' * firstcourierH * globalemptyisometry

    couriercorrelationspectrum = couriercorrelations |> crystalspectrum

    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap
    return Dict(
        :localcorrelations=>localcorrelations,
        :couriercorrelations=>purifiedcouriercorrelations,
        :courierH=>courierH,
        :filledH=>filledH,
        :emptyH=>emptyH,
        :filledcorrelations=>filledcorrelations,
        :emptycorrelations=>emptycorrelations,
        :wanniercourierisometry=>wanniercourierisometry,
        :globalfilledisometry=>globalfilledisometry,
        :globalemptyisometry=>globalemptyisometry,
        :wanniercourierstates=>wanniercourierstate|>RegionState,
        :minsvdforcourier=>minsvd,
        :traceofsv=>traceofsv,
        :rawcouriercorrelations=>couriercorrelations)
end

function gmerathirdstepwifentanglementcontour(secondcouriercorrelations,secondcourierH,nooffrozenmodes)
    secondgmeracrystalfock = secondcouriercorrelations|>getoutspace
    secondgmeracrystal::Crystal = secondgmeracrystalfock|>getcrystal
    secondgmeraspace::RealSpace = secondgmeracrystal|>getspace

    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    thirdcenter = [1/3,1/3] ∈ secondgmeraspace
    thirdhexagonalregion = gethexagonalregion(rot = refrot,crystal=secondgmeracrystal, center=thirdcenter, metricspace=secondgmeraspace)
    thirdhexagonalregionfock = quantize(thirdhexagonalregion,1)
    
    localcorrelations = regioncorrelations(secondcouriercorrelations,thirdhexagonalregionfock)
    localspectrum = localcorrelations|>eigspech

    frozenseedsoffset = rankedandgroupoffsets(thirdhexagonalregion,1)
    courierseedsoffset = thirdhexagonalregion-frozenseedsoffset
    courierseedsfock = quantize(courierseedsoffset,1)

    display(visualize(courierseedsoffset))

    nooffilledmodes = div(length(frozenseedsoffset),2)
    noofemptymodes = nooffilledmodes
    noofcouriermodes = length(courierseedsoffset)

    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes,noofcouriermodes,noofemptymodes])
    localcourierisometry = localstates[2]|>FockMap
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[3]|>FockMap

    courierseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in courierseedsfock)]
    localwanniercourier,minsvd,traceofsv=localwannierization(localcourierisometry, courierseeds)

    localwannierresults =  Dict(:localwanniercourier => localwanniercourier,:localfilledisometry => localfilledisometry,
                                    :localemptyisometry => localemptyisometry, :courierseeds => courierseeds)

    wannierinfos =  Dict(thirdhexagonalregionfock=>localwannierresults)

    wanniercourierisometry = globalwannierfunction(secondcouriercorrelations,wannierinfos[thirdhexagonalregionfock][:localwanniercourier])
    globalfilledisometry = globalwannierfunction(secondcouriercorrelations,wannierinfos[thirdhexagonalregionfock][:localfilledisometry ])
    globalemptyisometry = globalwannierfunction(secondcouriercorrelations,wannierinfos[thirdhexagonalregionfock][:localemptyisometry])

    @info "Computing local courier states..."
    leftrestrict = fourier(wanniercourierisometry|>getoutspace,thirdhexagonalregionfock) / (secondgmeracrystal|>vol|>sqrt)
    rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
    wanniercourierstate = leftrestrict' * wanniercourierisometry * rightrestrict

    couriercorrelations = wanniercourierisometry' * secondcouriercorrelations * wanniercourierisometry
    courierH = wanniercourierisometry' * secondcourierH * wanniercourierisometry

    filledcorrelations = globalfilledisometry' * secondcouriercorrelations * globalfilledisometry
    filledH = globalfilledisometry' * secondcourierH * globalfilledisometry

    emptycorrelations = globalemptyisometry' * secondcouriercorrelations * globalemptyisometry
    emptyH = globalemptyisometry' * secondcourierH * globalemptyisometry

    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

    return Dict(
        :localcorrelations=>localcorrelations,
        :couriercorrelations=>purifiedcouriercorrelations,
        :courierH=>courierH,
        :filledH=>filledH,
        :emptyH=>emptyH,
        :filledcorrelations=>filledcorrelations,
        :emptycorrelations=>emptycorrelations,
        :wanniercourierisometry=>wanniercourierisometry,
        :globalfilledisometry=>globalfilledisometry,
        :globalemptyisometry=>globalemptyisometry,
        :wanniercourierstates=>wanniercourierstate|>RegionState,
        :minsvdforcourier=>minsvd,
        :traceofsv=>traceofsv,
        :rawcouriercorrelations=>couriercorrelations)
end

blockedcorrelations,blockedH = blocking(correlations,H,2)

blockedH|>crystalspectrum|>visualize

firstgmera = gmerafirststepwifentanglementcontour(blockedcorrelations,blockedH,6)
firstcouriercorrelations = firstgmera[:couriercorrelations]
firstcourierH = firstgmera[:courierH]
firstgmera[:localcorrelations]|>eigspec|>visualize

secondgmera = gmerasecondstepwifentanglementcontour(firstcouriercorrelations,firstcourierH,6)
secondcouriercorrelations = secondgmera[:couriercorrelations]
secondcourierH = secondgmera[:courierH]
secondgmera[:localcorrelations]|>eigspec|>visualize

thirdgmera = gmerathirdstepwifentanglementcontour(secondcouriercorrelations,secondcourierH,6)
rg1correlations = thirdgmera[:couriercorrelations]
rg1H = thirdgmera[:courierH]
thirdgmera[:localcorrelations]|>eigspec|>visualize


rg1blockedcorrelations,rg1blockedH = blocking(rg1correlations,rg1H,2)

rg1firstgmera = gmerafirststepwifentanglementcontour(rg1blockedcorrelations,rg1blockedH,6)
rg1firstcouriercorrelations = rg1firstgmera[:couriercorrelations]
rg1firstcourierH = rg1firstgmera[:courierH]
rg1firstgmera[:localcorrelations]|>eigspec|>visualize

rg1secondgmera = gmerasecondstepwifentanglementcontour(rg1firstcouriercorrelations,rg1firstcourierH,6)
rg1secondcouriercorrelations = rg1secondgmera[:couriercorrelations]
rg1secondcourierH = rg1secondgmera[:courierH]
rg1secondgmera[:localcorrelations]|>eigspec|>visualize

rg1thirdgmera = gmerathirdstepwifentanglementcontour(rg1secondcouriercorrelations,rg1secondcourierH,6)
rg2correlations = rg1thirdgmera[:couriercorrelations]
rg2H = rg1thirdgmera[:courierH]
rg1thirdgmera[:localcorrelations]|>eigspec|>visualize

rg2blockedcorrelations,rg2blockedH = blocking(rg2correlations,rg2H,2)

rg2firstgmera = gmerafirststepwifentanglementcontour(rg2blockedcorrelations,rg2blockedH,6)
rg2firstcouriercorrelations = rg2firstgmera[:couriercorrelations]
rg2firstcourierH = rg2firstgmera[:courierH]
rg2firstgmera[:localcorrelations]|>eigspec|>visualize
rg2firstgmera[:emptycorrelations]|>crystalspectrum|>visualize

rg2secondgmera = gmerasecondstepwifentanglementcontour(rg2firstcouriercorrelations,rg2firstcourierH,6)
rg2secondcouriercorrelations = rg2secondgmera[:couriercorrelations]
rg2secondcourierH = rg2secondgmera[:courierH]
rg2secondgmera[:localcorrelations]|>eigspec|>visualize
rg2secondgmera[:filledcorrelations]|>crystalspectrum|>visualize

rg2thirdgmera = gmerathirdstepwifentanglementcontour(rg2secondcouriercorrelations,rg2secondcourierH,6)
rg3correlations = rg2thirdgmera[:couriercorrelations]
rg3H = rg2thirdgmera[:courierH]
rg2thirdgmera[:localcorrelations]|>eigspec|>visualize

rg3blockedcorrelations,rg3blockedH = blocking(rg3correlations,rg3H,2)

rg3firstgmera = gmerafirststepwifentanglementcontour(rg3blockedcorrelations,rg3blockedH,6)
rg3firstcouriercorrelations = rg3firstgmera[:couriercorrelations]
rg3firstcourierH = rg3firstgmera[:courierH]
rg3firstgmera[:localcorrelations]|>eigspec|>visualize

rg3secondgmera = gmerasecondstepwifentanglementcontour(rg3firstcouriercorrelations,rg3firstcourierH,6)
rg3secondcouriercorrelations = rg3secondgmera[:couriercorrelations]
rg3secondcourierH = rg3secondgmera[:courierH]
rg3secondgmera[:localcorrelations]|>eigspec|>visualize

rg3thirdgmera = gmerathirdstepwifentanglementcontour(rg3secondcouriercorrelations,rg3secondcourierH,6)
rg4correlations = rg3thirdgmera[:couriercorrelations]
rg4H = rg3thirdgmera[:courierH]
rg3thirdgmera[:localcorrelations]|>eigspec|>visualize

rg4blockedcorrelations,rg4blockedH = blocking(rg4correlations,rg4H,2)

rg4firstgmera = gmerafirststepwifentanglementcontour(rg4blockedcorrelations,rg4blockedH,6)
rg4firstcouriercorrelations = rg4firstgmera[:couriercorrelations]
rg4firstcourierH = rg4firstgmera[:courierH]
rg4firstgmera[:localcorrelations]|>eigspec|>visualize
rg4firstgmera[:courierH]|>crystalspectrum|>visualize

rg4secondgmera = gmerasecondstepwifentanglementcontour(rg4firstcouriercorrelations,rg4firstcourierH,6)
rg4secondcouriercorrelations = rg4secondgmera[:couriercorrelations]
rg4secondcourierH = rg4secondgmera[:courierH]
rg4secondgmera[:localcorrelations]|>eigspec|>visualize
rg4secondgmera[:courierH]|>crystalspectrum|>visualize

rg4thirdgmera = gmerathirdstepwifentanglementcontour(rg4secondcouriercorrelations,rg4secondcourierH,6)
rg5correlations = rg4thirdgmera[:couriercorrelations]
rg5H = rg4thirdgmera[:courierH]
rg4thirdgmera[:localcorrelations]|>eigspec|>visualize

rg5blockedcorrelations,rg5blockedH = blocking(rg5correlations,rg5H,2)

rg5firstgmera = gmerafirststepwifentanglementcontour(rg5blockedcorrelations,rg5blockedH,6)
rg5firstcouriercorrelations = rg4firstgmera[:couriercorrelations]
rg5firstcourierH = rg4firstgmera[:courierH]
rg5firstgmera[:localcorrelations]|>eigspec|>visualize
rg5firstgmera[:courierH]|>crystalspectrum|>visualize

