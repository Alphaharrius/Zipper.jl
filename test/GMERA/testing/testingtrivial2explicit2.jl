using Zipper
using LinearAlgebra
using Plots

setmaxthreads(8)

function focktraceL1norm(fockmap,systemsize)
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end

systemsize=40
correlations,H = generatesystem( -0+0im, 0+0im,-1 + 0im,0im,systemsize)
crystalfock = correlations|>getoutspace
crystal = correlations|>getoutspace|>getcrystal
space::RealSpace = crystal|>getspace
refrot = inv([2/3 -1/3; -1/3 2/3]')
refcenter = [0,0] ∈ space

intermediatehexagonalregion = gethexagonalregion(rot = refrot,crystal=crystal, center=refcenter, metricspace=space)

scale = Scale([2 0; 0 2], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * crystalfock
@info("Performing rgblocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'

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

innerhexagonalregion = checkintersectbyeculidean(intermediatehexagonalregion,firsthexagonalregion)
innerhexagonalregionfock = quantize(innerhexagonalregion,1)

#siteAregion1
siteAregion1 = intersect(intersect(firsthexagonalregion,firstrgshiftedhexagonalregion1 ),firstrgshiftedhexagonalregion2)
steAregion1fock = quantize(siteAregion1,1)
regioncorrelations(blockedcorrelations,steAregion1fock)|>eigspech|>visualize
siteAregion1closetriangle1 =  intersect(siteAregion1,innerhexagonalregion)

siteAregion1bdytriangle1 = siteAregion1closetriangle1.+(1/2)*rgshiftedcenter1
siteAregion1bdytriangle2 = siteAregion1closetriangle1.+(1/2)*rgshiftedcenter2
siteAregion1bdytriangle1fock = quantize(siteAregion1bdytriangle1,1)
siteAregion1bdytriangle2fock = quantize(siteAregion1bdytriangle2,1)

siteAregion1centertriangle = siteAregion1 - siteAregion1closetriangle1 - siteAregion1bdytriangle1 - siteAregion1bdytriangle2
siteAregion1centertrianglefock = quantize(siteAregion1centertriangle,1)

#siteAregion2
siteAregion2 = c3*siteAregion1
siteAregion2closetriangle1 =  c3*siteAregion1closetriangle1

siteAregion2bdytriangle1 = c3*siteAregion1bdytriangle1
siteAregion2bdytriangle2 = c3*siteAregion1bdytriangle2
siteAregion2bdytriangle1fock = quantize(siteAregion2bdytriangle1,1)
siteAregion2bdytriangle2fock = quantize(siteAregion2bdytriangle2,1)

siteAregion2centertriangle = c3*siteAregion1centertriangle
siteAregion2centertrianglefock = quantize(siteAregion2centertriangle,1)

#siteAregion3
siteAregion3 = c3*siteAregion2
siteAregion3closetriangle1 =  c3*siteAregion2closetriangle1

siteAregion3bdytriangle1 = c3*siteAregion2bdytriangle1
siteAregion3bdytriangle2 = c3*siteAregion2bdytriangle2
siteAregion3bdytriangle1fock = quantize(siteAregion3bdytriangle1,1)
siteAregion3bdytriangle2fock = quantize(siteAregion3bdytriangle2,1)

siteAregion3centertriangle = c3*siteAregion2centertriangle
siteAregion3centertrianglefock = quantize(siteAregion3centertriangle,1)

#siteBregion1
siteBregion1 = c6*siteAregion1
siteBregion1closetriangle1 =  c6*siteAregion1closetriangle1

siteBregion1bdytriangle1 = c6*siteAregion1bdytriangle1
siteBregion1bdytriangle2 = c6*siteAregion1bdytriangle2
siteBregion1bdytriangle1fock = quantize(siteBregion1bdytriangle1,1)
siteBregion1bdytriangle2fock = quantize(siteBregion1bdytriangle2,1)

siteBregion1centertriangle = c6*siteAregion1centertriangle
siteBregion1centertrianglefock = quantize(siteBregion1centertriangle,1)

#siteBregion2
siteBregion2 = c3*siteBregion1
siteBregion2closetriangle1 =  c3*siteBregion1closetriangle1

siteBregion2bdytriangle1 = c3*siteBregion1bdytriangle1
siteBregion2bdytriangle2 = c3*siteBregion1bdytriangle2
siteBregion2bdytriangle1fock = quantize(siteBregion2bdytriangle1,1)
siteBregion2bdytriangle2fock = quantize(siteBregion2bdytriangle2,1)

siteBregion2centertriangle = c3*siteBregion1centertriangle
siteBregion2centertrianglefock = quantize(siteBregion2centertriangle,1)

#siteBregion3
siteBregion3 = (c3)*siteBregion2
siteBregion3closetriangle1 =  c3*siteBregion2closetriangle1

siteBregion3bdytriangle1 = c3*siteBregion2bdytriangle1
siteBregion3bdytriangle2 = c3*siteBregion2bdytriangle2
siteBregion3bdytriangle1fock = quantize(siteBregion3bdytriangle1,1)
siteBregion3bdytriangle2fock = quantize(siteBregion3bdytriangle2,1)

siteBregion3centertriangle = c3*siteBregion2centertriangle
siteBregion3centertrianglefock = quantize(siteBregion3centertriangle,1)

localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
localspectrum = localcorrelations|>eigspec
# localspectrum|>visualize

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3,18,3])
localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[3]|>FockMap

localfilledcorrelations = idmap(localcorrelations|>getoutspace)-(localfilledisometry*localfilledisometry')
localemptycorrelations = idmap(localcorrelations|>getoutspace)-(localemptyisometry*localemptyisometry')
localcouriercorrelations = idmap(localcorrelations|>getoutspace)-(localcourierisometry*localcourierisometry')

localfilledcorrelations[innerhexagonalregionfock,innerhexagonalregionfock]|>eigspech|>visualize
localcouriercorrelations[siteBregion3bdytriangle2fock,siteBregion3bdytriangle2fock]|>eigspech|>visualize
localcouriercorrelations[siteAregion1centertrianglefock,siteAregion1centertrianglefock]|>eigspech|>visualize

innerhexagonalfilledseeds = getregionstates(localcorrelations=localfilledcorrelations[innerhexagonalregionfock,innerhexagonalregionfock],grouping=[3])[1]
centerfilledseedsB1 =  getregionstates(localcorrelations=localfilledcorrelations[siteBregion1centertrianglefock,siteBregion1centertrianglefock],grouping=[1])[1]
centerfilledseedsB2 =  getregionstates(localcorrelations=localfilledcorrelations[siteBregion2centertrianglefock,siteBregion2centertrianglefock],grouping=[1])[1]
centerfilledseedsB3 =  getregionstates(localcorrelations=localfilledcorrelations[siteBregion3centertrianglefock,siteBregion3centertrianglefock],grouping=[1])[1]
filledseeds = innerhexagonalfilledseeds+centerfilledseedsA1+centerfilledseedsA2+centerfilledseedsA3
filledseeds = filledseeds|>FockMap
filledextension = extensionmap(filledseeds|>getoutspace,localfilledisometry|>getoutspace)
extendedfilledseeds = filledextension'*filledseeds

localwannierfilled = localwannierization(localfilledisometry, extendedfilledseeds)

innerhexagonalemptyseeds = getregionstates(localcorrelations=localemptycorrelations[innerhexagonalregionfock,innerhexagonalregionfock],grouping=[3])[1]
centeremptyseedsA1 =  getregionstates(localcorrelations=localemptycorrelations[siteAregion1centertrianglefock,siteAregion1centertrianglefock],grouping=[1])[1]
centeremptyseedsA2 =  getregionstates(localcorrelations=localemptycorrelations[siteAregion2centertrianglefock,siteAregion2centertrianglefock],grouping=[1])[1]
centeremptyseedsA3 =  getregionstates(localcorrelations=localemptycorrelations[siteAregion3centertrianglefock,siteAregion3centertrianglefock],grouping=[1])[1]
emptyseeds = innerhexagonalemptyseeds+centeremptyseedsA1+centeremptyseedsA2+centeremptyseedsA3
emptyseeds  = emptyseeds |>FockMap
emptyextension = extensionmap(emptyseeds|>getoutspace,localemptyisometry|>getoutspace)
extendedemptyseeds = emptyextension'*emptyseeds

localwannierempty = localwannierization(localemptyisometry, extendedemptyseeds)

bdytriangle1courierseedsA1 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion1bdytriangle1fock,siteAregion1bdytriangle1fock]), grouping=[1])[1]
bdytriangle2courierseedsA1 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion1bdytriangle2fock,siteAregion1bdytriangle2fock]), grouping=[1])[1]
bdytriangle1courierseedsA2 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion2bdytriangle1fock,siteAregion2bdytriangle1fock]), grouping=[1])[1]
bdytriangle2courierseedsA2 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion2bdytriangle2fock,siteAregion2bdytriangle2fock]), grouping=[1])[1]
bdytriangle1courierseedsA3 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion3bdytriangle1fock,siteAregion3bdytriangle1fock]), grouping=[1])[1]
bdytriangle2courierseedsA3 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion3bdytriangle2fock,siteAregion3bdytriangle2fock]), grouping=[1])[1]
bdytriangle1courierseedsB1 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion1bdytriangle1fock,siteBregion1bdytriangle1fock]), grouping=[1])[1]
bdytriangle2courierseedsB1 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion1bdytriangle2fock,siteBregion1bdytriangle2fock]), grouping=[1])[1]
bdytriangle1courierseedsB2 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion2bdytriangle1fock,siteBregion2bdytriangle1fock]), grouping=[1])[1]
bdytriangle2courierseedsB2 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion2bdytriangle2fock,siteBregion2bdytriangle2fock]), grouping=[1])[1]
bdytriangle1courierseedsB3 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion3bdytriangle1fock,siteBregion3bdytriangle1fock]), grouping=[1])[1]
bdytriangle2courierseedsB3 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion3bdytriangle2fock,siteBregion3bdytriangle2fock]), grouping=[1])[1]

courierseeds = (bdytriangle1courierseedsA1+bdytriangle2courierseedsA1+bdytriangle1courierseedsA2+bdytriangle2courierseedsA2+bdytriangle1courierseedsA3+bdytriangle2courierseedsA3
                +bdytriangle1courierseedsB1+bdytriangle2courierseedsB1+bdytriangle1courierseedsB2+bdytriangle2courierseedsB2+bdytriangle1courierseedsB3+bdytriangle2courierseedsB3)|>FockMap
courierextension = extensionmap(courierseeds|>getoutspace,localcourierisometry|>getoutspace)
extendedcourierseeds = courierextension'*courierseeds
localwanniercourier = localwannierization(localcourierisometry, extendedcourierseeds)

localwannierresults =  Dict(:localwanniercourier => localwanniercourier, :localwannierfilled => localwannierfilled,
                                :localwannierempty => localwannierempty, :courierseeds => courierseeds)

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

    localcorrelations = regioncorrelations(blockedcorrelations,shiftedfirsthexagonalregionfock)
    localspectrum = localcorrelations|>eigspec
    localspectrum|>visualize

    localsiteA1correlations = regioncorrelations(blockedcorrelations,shiftedsiteAregionfock1)
    localsiteA2correlations = regioncorrelations(blockedcorrelations,shiftedsiteAregionfock2)
    localsiteA3correlations = regioncorrelations(blockedcorrelations,shiftedsiteAregionfock3)
    localsiteB1correlations = regioncorrelations(blockedcorrelations,shiftedsiteBregionfock1)
    localsiteB2correlations = regioncorrelations(blockedcorrelations,shiftedsiteBregionfock2)
    localsiteB3correlations = regioncorrelations(blockedcorrelations,shiftedsiteBregionfock3)

    # localA1states = getregionstates(localcorrelations=localsiteA1correlations, grouping=[1, 2, 1])
    # localA2states = getregionstates(localcorrelations=localsiteA2correlations, grouping=[1, 2, 1])
    # localA3states = getregionstates(localcorrelations=localsiteA3correlations, grouping=[1, 2, 1])
    # localB1states = getregionstates(localcorrelations=localsiteB1correlations, grouping=[1, 2, 1])
    # localB2states = getregionstates(localcorrelations=localsiteB2correlations, grouping=[1, 2, 1])
    # localB3states = getregionstates(localcorrelations=localsiteB3correlations, grouping=[1, 2, 1])

    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[6,12,6])
    localcourierisometry = localstates[2]|>FockMap
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[3]|>FockMap

    # filledseedsA1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteAregionfock1,shiftedsiteAregionfock1]), grouping=[3,1])[2]
    # filledseedsA2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteAregionfock2,shiftedsiteAregionfock2]), grouping=[3,1])[2]
    # filledseedsA3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteAregionfock3,shiftedsiteAregionfock3]), grouping=[3,1])[2]
    # filledseedsB1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteBregionfock1,shiftedsiteBregionfock1]), grouping=[3,1])[2]
    # filledseedsB2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteBregionfock2,shiftedsiteBregionfock2]), grouping=[3,1])[2]
    # filledseedsB3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteBregionfock3,shiftedsiteBregionfock3]), grouping=[3,1])[2]

    # filledseeds = (filledseedsA1+filledseedsA2+filledseedsA3+filledseedsB1+filledseedsB2+filledseedsB3)|>FockMap
    # localwannierfilled = localwannierization(localfilledisometry, filledseeds)

    # emptyseedsA1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteAregionfock1,shiftedsiteAregionfock1]), grouping=[3,1])[2]
    # emptyseedsA2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteAregionfock2,shiftedsiteAregionfock2]), grouping=[3,1])[2]
    # emptyseedsA3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteAregionfock3,shiftedsiteAregionfock3]), grouping=[3,1])[2]
    # emptyseedsB1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteBregionfock1,shiftedsiteBregionfock1]), grouping=[3,1])[2]
    # emptyseedsB2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteBregionfock2,shiftedsiteBregionfock2]), grouping=[3,1])[2]
    # emptyseedsB3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteBregionfock3,shiftedsiteBregionfock3]), grouping=[3,1])[2]

    # emptyseeds = (emptyseedsA1+emptyseedsA2+emptyseedsA3+emptyseedsB1+emptyseedsB2+emptyseedsB3)|>FockMap
    # localwannierempty = localwannierization(localemptyisometry, emptyseeds)

    courierseedsA1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteAregionfock1,shiftedsiteAregionfock1]), grouping=[2,2])[2]
    courierseedsA2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteAregionfock2,shiftedsiteAregionfock2]), grouping=[2,2])[2]
    courierseedsA3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteAregionfock3,shiftedsiteAregionfock3]), grouping=[2,2])[2]
    courierseedsB1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteBregionfock1,shiftedsiteBregionfock1]), grouping=[2,2])[2]
    courierseedsB2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteBregionfock2,shiftedsiteBregionfock2]), grouping=[2,2])[2]
    courierseedsB3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteBregionfock3,shiftedsiteBregionfock3]), grouping=[2,2])[2]

    courierseeds = (courierseedsA1+courierseedsA2+courierseedsA3+courierseedsB1+courierseedsB2+courierseedsB3)|>FockMap
    localwanniercourier = localwannierization(localcourierisometry, courierseeds)


    shiftedlocalwannierresults =  Dict(:localwanniercourier => localwanniercourier, :localwannierfilled => localwannierfilled,
                                        :localwannierempty => localwannierempty, :courierseeds => courierseeds)
    wannierinfos[shiftedfirsthexagonalregionfock] =  shiftedlocalwannierresults
end

ref = [quantize(firsthexagonalregion.+hexagonalcenter,1) for hexagonalcenter in firstshiftedhexagonalcenters]
firsthexagonalregionfocklist = [firsthexagonalregionfock,ref...]

extendedwannierizedcourier =  sum(wannierinfos[regionfock][:localwanniercourier] for regionfock in firsthexagonalregionfocklist)
extendedwannierizedfilled =  sum(wannierinfos[regionfock][:localwannierfilled] for regionfock in firsthexagonalregionfocklist)
extendedwannierizedempty =  sum(wannierinfos[regionfock][:localwannierempty] for regionfock in firsthexagonalregionfocklist)
 
                
origin = [0, 0] ∈ blockedspace
refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

wanniercourierisometry1 = globalwannierfunction(blockedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
wannierfilledisometry1 = globalwannierfunction(blockedcorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
wannieremptyisometry1 = globalwannierfunction(blockedcorrelations,extendedwannierizedempty[:,refunictcellfockempty])

@info "Computing local courier states..."
leftrestrict = fourier(wanniercourierisometry1|>getoutspace, firsthexagonalregionfock) / (blockedcrystal|>vol|>sqrt)
rightrestrict = fourier(wanniercourierisometry1|>getinspace, wanniercourierisometry1|>getinspace|>unitcellfock|>RegionFock)
wanniercourierstate = leftrestrict' * wanniercourierisometry1 * rightrestrict

couriercorrelations = wanniercourierisometry1' * blockedcorrelations * wanniercourierisometry1
# display(couriercorrelations|>eigspech|>visualize)

filledcorrelations = wannierfilledisometry1' * blockedcorrelations * wannierfilledisometry1
# display(filledcorrelations|>eigspech|>visualize)

emptycorrelations = wannieremptyisometry1' * blockedcorrelations * wannieremptyisometry1

couriercorrelationspectrum = couriercorrelations |> crystalspectrum

purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

firstgmeracorrelations = purifiedcouriercorrelations

firstgmeracrystalfock = firstgmeracorrelations|>getoutspace
firstgmeracrystal::Crystal = firstgmeracrystalfock|>getcrystal
firstgmeraspace::RealSpace = firstgmeracrystal|>getspace

noofflavourpermode = 2

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
localspectrum|>visualize

localsiteA1correlations = regioncorrelations(firstgmeracorrelations,siteAregionfock1)
localsiteA2correlations = regioncorrelations(firstgmeracorrelations,siteAregionfock2)
localsiteA3correlations = regioncorrelations(firstgmeracorrelations,siteAregionfock3)
localsiteB1correlations = regioncorrelations(firstgmeracorrelations,siteBregionfock1)
localsiteB2correlations = regioncorrelations(firstgmeracorrelations,siteBregionfock2)
localsiteB3correlations = regioncorrelations(firstgmeracorrelations,siteBregionfock3)
# localsiteB1correlations|>eigspec|>visualize

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3,6,3])
localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[3]|>FockMap

(localcourierisometry*localcourierisometry')[siteAregionfock1,siteAregionfock1]|>eigspech|>visualize

# filledseedsA1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock1,siteAregionfock1]), grouping=[1,1])[2]
# filledseedsA2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock2,siteAregionfock2]), grouping=[1,1])[2]
# filledseedsA3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock3,siteAregionfock3]), grouping=[1,1])[2]
# filledseedsB1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteBregionfock1,siteBregionfock1]), grouping=[1,1])[2]
# filledseedsB2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteBregionfock2,siteBregionfock2]), grouping=[1,1])[2]
# filledseedsB3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteBregionfock3,siteBregionfock3]), grouping=[1,1])[2]

# filledseeds = (filledseedsA1+filledseedsA2+filledseedsA3+filledseedsB1+filledseedsB2+filledseedsB3)|>FockMap
# localwannierfilled = localwannierization(localfilledisometry, filledseeds)

# emptyseedsA1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteAregionfock1,siteAregionfock1]), grouping=[1,1])[2]
# emptyseedsA2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteAregionfock2,siteAregionfock2]), grouping=[1,1])[2]
# emptyseedsA3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteAregionfock3,siteAregionfock3]), grouping=[1,1])[2]
# emptyseedsB1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteBregionfock1,siteBregionfock1]), grouping=[1,1])[2]
# emptyseedsB2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteBregionfock2,siteBregionfock2]), grouping=[1,1])[2]
# emptyseedsB3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteBregionfock3,siteBregionfock3]), grouping=[1,1])[2]

# emptyseeds = (emptyseedsA1+emptyseedsA2+emptyseedsA3+emptyseedsB1+emptyseedsB2+emptyseedsB3)|>FockMap
# localwannierempty = localwannierization(localemptyisometry, emptyseeds)

courierseedsA1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteAregionfock1,siteAregionfock1]), grouping=[1,1])[2]
courierseedsA2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteAregionfock2,siteAregionfock2]), grouping=[1,1])[2]
courierseedsA3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteAregionfock3,siteAregionfock3]), grouping=[1,1])[2]
courierseedsB1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteBregionfock1,siteBregionfock1]), grouping=[1,1])[2]
courierseedsB2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteBregionfock2,siteBregionfock2]), grouping=[1,1])[2]
courierseedsB3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteBregionfock3,siteBregionfock3]), grouping=[1,1])[2]

courierseeds = (courierseedsA1+courierseedsA2+courierseedsA3+courierseedsB1+courierseedsB2+courierseedsB3)|>FockMap
localwanniercourier = localwannierization(localcourierisometry, courierseeds)

localwannierresults =  Dict(:localwanniercourier => localwanniercourier, :localwannierfilled => localfilledisometry,
                                :localwannierempty => localemptyisometry, :courierseeds => courierseeds)

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

    localcorrelations = regioncorrelations(firstgmeracorrelations,shiftedsecondhexagonalregionfock)
    localspectrum = localcorrelations|>eigspec

    localsiteA1correlations = regioncorrelations(firstgmeracorrelations,siteAregionfock1)
    localsiteA2correlations = regioncorrelations(firstgmeracorrelations,siteAregionfock2)
    localsiteA3correlations = regioncorrelations(firstgmeracorrelations,siteAregionfock3)
    localsiteB1correlations = regioncorrelations(firstgmeracorrelations,siteBregionfock1)
    localsiteB2correlations = regioncorrelations(firstgmeracorrelations,siteBregionfock2)
    localsiteB3correlations = regioncorrelations(firstgmeracorrelations,siteBregionfock3)
    localsiteB1correlations|>eigspec|>visualize

    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3,6,3])
    localcourierisometry = localstates[2]|>FockMap
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[3]|>FockMap

    # filledseedsA1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteAregionfock1,shiftedsiteAregionfock1]), grouping=[1,1])[2]
    # filledseedsA2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteAregionfock2,shiftedsiteAregionfock2]), grouping=[1,1])[2]
    # filledseedsA3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteAregionfock3,shiftedsiteAregionfock3]), grouping=[1,1])[2]
    # filledseedsB1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteBregionfock1,shiftedsiteBregionfock1]), grouping=[1,1])[2]
    # filledseedsB2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteBregionfock2,shiftedsiteBregionfock2]), grouping=[1,1])[2]
    # filledseedsB3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteBregionfock3,shiftedsiteBregionfock3]), grouping=[1,1])[2]

    # filledseeds = (filledseedsA1+filledseedsA2+filledseedsA3+filledseedsB1+filledseedsB2+filledseedsB3)|>FockMap
    # localwannierfilled = localwannierization(localfilledisometry, filledseeds)

    # emptyseedsA1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteAregionfock1,shiftedsiteAregionfock1]), grouping=[1,1])[2]
    # emptyseedsA2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteAregionfock2,shiftedsiteAregionfock2]), grouping=[1,1])[2]
    # emptyseedsA3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteAregionfock3,shiftedsiteAregionfock3]), grouping=[1,1])[2]
    # emptyseedsB1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteBregionfock1,shiftedsiteBregionfock1]), grouping=[1,1])[2]
    # emptyseedsB2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteBregionfock2,shiftedsiteBregionfock2]), grouping=[1,1])[2]
    # emptyseedsB3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteBregionfock3,shiftedsiteBregionfock3]), grouping=[1,1])[2]

    # emptyseeds = (emptyseedsA1+emptyseedsA2+emptyseedsA3+emptyseedsB1+emptyseedsB2+emptyseedsB3)|>FockMap
    # localwannierempty = localwannierization(localemptyisometry, emptyseeds)

    courierseedsA1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteAregionfock1,shiftedsiteAregionfock1]), grouping=[1,1])[2]
    courierseedsA2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteAregionfock2,shiftedsiteAregionfock2]), grouping=[1,1])[2]
    courierseedsA3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteAregionfock3,shiftedsiteAregionfock3]), grouping=[1,1])[2]
    courierseedsB1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteBregionfock1,shiftedsiteBregionfock1]), grouping=[1,1])[2]
    courierseedsB2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteBregionfock2,shiftedsiteBregionfock2]), grouping=[1,1])[2]
    courierseedsB3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteBregionfock3,shiftedsiteBregionfock3]), grouping=[1,1])[2]

    courierseeds = (courierseedsA1+courierseedsA2+courierseedsA3+courierseedsB1+courierseedsB2+courierseedsB3)|>FockMap
    localwanniercourier = localwannierization(localcourierisometry, courierseeds)
    

    shiftedlocalwannierresults = Dict(:localwanniercourier => localwanniercourier, :localwannierfilled => localfilledisometry,
                                    :localwannierempty => localemptyisometry, :courierseeds => courierseeds)
    wannierinfos[shiftedsecondhexagonalregionfock] =  shiftedlocalwannierresults
end

ref = [quantize(secondhexagonalregion.+hexagonalcenter,noofflavourpermode) for hexagonalcenter in secondshiftedhexagonalcenterlist]
secondhexagonalregionfocklist = [secondhexagonalregionfock,ref...]

extendedwannierizedcourier =  sum(wannierinfos[regionfock][:localwanniercourier] for regionfock in secondhexagonalregionfocklist)
extendedwannierizedfilled =  sum(wannierinfos[regionfock][:localwannierfilled] for regionfock in secondhexagonalregionfocklist)
extendedwannierizedempty =  sum(wannierinfos[regionfock][:localwannierempty] for regionfock in secondhexagonalregionfocklist) 

origin = [0, 0] ∈ firstgmeraspace
refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

wanniercourierisometry2 = globalwannierfunction(firstgmeracorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
wannierfilledisometry2 = globalwannierfunction(firstgmeracorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
wannieremptyisometry2 = globalwannierfunction(firstgmeracorrelations,extendedwannierizedempty[:,refunictcellfockempty])

@info "Computing local courier states..."
leftrestrict = fourier(wanniercourierisometry2|>getoutspace,secondhexagonalregionfock) / (firstgmeracrystal|>vol|>sqrt)
rightrestrict = fourier(wanniercourierisometry2|>getinspace, wanniercourierisometry2|>getinspace|>unitcellfock|>RegionFock)
wanniercourierstate = leftrestrict' * wanniercourierisometry2 * rightrestrict

couriercorrelations = wanniercourierisometry2' * firstgmeracorrelations * wanniercourierisometry2
# display(couriercorrelations|>eigspech|>visualize)

filledcorrelations = wannierfilledisometry2' * firstgmeracorrelations * wannierfilledisometry2
# display(filledcorrelations|>eigspech|>visualize)

emptycorrelations = wannieremptyisometry2' * firstgmeracorrelations * wannieremptyisometry2

couriercorrelationspectrum = couriercorrelations |> crystalspectrum
purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

secondgmeracorrelations = purifiedcouriercorrelations

secondgmeracrystalfock = secondgmeracorrelations|>getoutspace
secondgmeracrystal::Crystal = secondgmeracrystalfock|>getcrystal
secondgmeraspace::RealSpace = secondgmeracrystal|>getspace

refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2
noofflavourpermode=1

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

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3,3])
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[2]|>FockMap

(localfilledisometry*localfilledisometry')[siteBregionfock1,siteBregionfock1]|>eigspech|>visualize

# filledseedsA1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock1,siteAregionfock1]), grouping=[1])[1]
# filledseedsA2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock2,siteAregionfock2]), grouping=[1])[1]
# filledseedsA3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock3,siteAregionfock3]), grouping=[1])[1]

# filledseeds = (filledseedsA1+filledseedsA2+filledseedsA3)|>FockMap
# restrictfilledfock = filledseeds|>getoutspace|>RegionFock
# localwannierfilled = localwannierization(localfilledisometry[restrictfilledfock,:], filledseeds)
# localwannierfilled'*localcorrelations*localwannierfilled


# emptyseedsA1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteBregionfock1,siteBregionfock1]), grouping=[1])[1]
# emptyseedsA2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteBregionfock2,siteBregionfock2]), grouping=[1])[1]
# emptyseedsA3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteBregionfock3,siteBregionfock3]), grouping=[1])[1]

# emptyseeds = (emptyseedsA1+emptyseedsA2+emptyseedsA3)|>FockMap
# restrictemptyfock = emptyseeds|>getoutspace|>RegionFock
# localwannierempty = localwannierization(localemptyisometry[restrictemptyfock,:], emptyseeds)


localwannierresults =  Dict(:localfilledisometry => localfilledisometry,
                            :localemptyisometry => localemptyisometry)

wannierinfos =  Dict(thirdhexagonalregionfock=>localwannierresults)
wannierfilledisometry3 = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localfilledisometry ])
wannieremptyisometry3 = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localemptyisometry])
filledcorrelations = wannierfilledisometry3' * secondgmeracorrelations * wannierfilledisometry3
emptycorrelations = wannieremptyisometry3' * secondgmeracorrelations * wannieremptyisometry3
visualize(emptycorrelations|>eigspech)

firstapproximation = wannieremptyisometry1*wannieremptyisometry1'
secondapproximation = (wanniercourierisometry1*wannieremptyisometry2)*(wanniercourierisometry1*wannieremptyisometry2)'
thirdapproximation = (wanniercourierisometry1*wanniercourierisometry2*wannieremptyisometry3)*(wanniercourierisometry1*wanniercourierisometry2*wannieremptyisometry3)'

arrpoximatecorrelationtwo = firstapproximation + secondapproximation + thirdapproximation

difftwo = blockedcorrelations-arrpoximatecorrelationtwo
difftwonorm = focktraceL1norm(difftwo,9600)
difftwonorm = (difftwonorm/6144)*9600