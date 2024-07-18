using Zipper
using LinearAlgebra
using Plots

setmaxthreads(8)

systemsize=30
correlations,H = generatesystem(0im, 0im,-1 + 0im,0im,systemsize)
crystalfock = correlations|>getoutspace

scale = Scale([10 0; 0 10], crystalfock|>getcrystal|>getspace)
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
localspectrum|>visualize

localsiteA1correlations = regioncorrelations(blockedcorrelations,siteAregionfock1)
localsiteA2correlations = regioncorrelations(blockedcorrelations,siteAregionfock2)
localsiteA3correlations = regioncorrelations(blockedcorrelations,siteAregionfock3)
localsiteB1correlations = regioncorrelations(blockedcorrelations,siteBregionfock1)
localsiteB2correlations = regioncorrelations(blockedcorrelations,siteBregionfock2)
localsiteB3correlations = regioncorrelations(blockedcorrelations,siteBregionfock3)
localsiteA1correlations|>eigspech|>visualize
visualize(siteAregion1)

# localA1states = getregionstates(localcorrelations=localsiteA1correlations, grouping=[1, 2, 1])
# localA2states = getregionstates(localcorrelations=localsiteA2correlations, grouping=[1, 2, 1])
# localA3states = getregionstates(localcorrelations=localsiteA3correlations, grouping=[1, 2, 1])
# localB1states = getregionstates(localcorrelations=localsiteB1correlations, grouping=[1, 2, 1])
# localB2states = getregionstates(localcorrelations=localsiteB2correlations, grouping=[1, 2, 1])
# localB3states = getregionstates(localcorrelations=localsiteB3correlations, grouping=[1, 2, 1])

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[270,60,270])
localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[3]|>FockMap
# localfrozenisometry = (localstates[1]+localstates[3])|>FockMap

(localfilledisometry*localfilledisometry')[siteAregionfock1,siteAregionfock1]|>eigspech|>visualize
(localcourierisometry*localcourierisometry')[siteAregionfock1,siteAregionfock1]|>eigspech|>visualize
# (localfrozenisometry*localfrozenisometry')[siteAregionfock1,siteAregionfock1]|>eigspech|>visualize

# frozenseedsA1 = getregionstates(localcorrelations=((localfrozenisometry*localfrozenisometry')[siteAregionfock1,siteAregionfock1]), grouping=[12,52])[2]
# frozenseedsA2 = getregionstates(localcorrelations=((localfrozenisometry*localfrozenisometry')[siteAregionfock2,siteAregionfock2]), grouping=[12,52])[2]
# frozenseedsA3 = getregionstates(localcorrelations=((localfrozenisometry*localfrozenisometry')[siteAregionfock3,siteAregionfock3]), grouping=[12,52])[2]
# frozenseedsB1 = getregionstates(localcorrelations=((localfrozenisometry*localfrozenisometry')[siteBregionfock1,siteBregionfock1]), grouping=[12,52])[2]
# frozenseedsB2 = getregionstates(localcorrelations=((localfrozenisometry*localfrozenisometry')[siteBregionfock2,siteBregionfock2]), grouping=[12,52])[2]
# frozenseedsB3 = getregionstates(localcorrelations=((localfrozenisometry*localfrozenisometry')[siteBregionfock3,siteBregionfock3]), grouping=[12,52])[2]

# frozenseeds = (frozenseedsA1+frozenseedsA2+frozenseedsA3+frozenseedsB1+frozenseedsB2+frozenseedsB3)|>FockMap
# localwannierfrozen = localwannierization(localfrozenisometry, frozenseeds)

# scatter(sort([norm(d) for d in diag((localwannierfrozen'*localcorrelations*localwannierfrozen).rep)]))

# filledseedsA1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock1,siteAregionfock1]), grouping=[15,1])[2]
# filledseedsA2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock2,siteAregionfock2]), grouping=[15,1])[2]
# filledseedsA3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock3,siteAregionfock3]), grouping=[15,1])[2]
# filledseedsB1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteBregionfock1,siteBregionfock1]), grouping=[15,1])[2]
# filledseedsB2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteBregionfock2,siteBregionfock2]), grouping=[15,1])[2]
# filledseedsB3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteBregionfock3,siteBregionfock3]), grouping=[15,1])[2]

# filledseeds = (filledseedsA1+filledseedsA2+filledseedsA3+filledseedsB1+filledseedsB2+filledseedsB3)|>FockMap
# localwannierfilled = localwannierization(localfilledisometry, filledseeds)

# emptyseedsA1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteAregionfock1,siteAregionfock1]), grouping=[10,6])[2]
# emptyseedsA2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteAregionfock2,siteAregionfock2]), grouping=[10,6])[2]
# emptyseedsA3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteAregionfock3,siteAregionfock3]), grouping=[10,6])[2]
# emptyseedsB1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteBregionfock1,siteBregionfock1]), grouping=[10,6])[2]
# emptyseedsB2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteBregionfock2,siteBregionfock2]), grouping=[10,6])[2]
# emptyseedsB3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteBregionfock3,siteBregionfock3]), grouping=[10,6])[2]

# emptyseeds = (emptyseedsA1+emptyseedsA2+emptyseedsA3+emptyseedsB1+emptyseedsB2+emptyseedsB3)|>FockMap
# localwannierempty = localwannierization(localemptyisometry, emptyseeds)

courierseedsA1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteAregionfock1,siteAregionfock1]), grouping=[90,10])[2]
courierseedsA2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteAregionfock2,siteAregionfock2]), grouping=[90,10])[2]
courierseedsA3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteAregionfock3,siteAregionfock3]), grouping=[90,10])[2]
courierseedsB1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteBregionfock1,siteBregionfock1]), grouping=[90,10])[2]
courierseedsB2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteBregionfock2,siteBregionfock2]), grouping=[90,10])[2]
courierseedsB3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteBregionfock3,siteBregionfock3]), grouping=[90,10])[2]

courierseeds = (courierseedsA1+courierseedsA2+courierseedsA3+courierseedsB1+courierseedsB2+courierseedsB3)|>FockMap
localwanniercourier = localwannierization(localcourierisometry, courierseeds)

localwannierresults =  Dict(:localwanniercourier => localwanniercourier, :localwannierfilled => localfilledisometry,
                                :localwannierempty => localemptyisometry, :courierseeds => courierseeds)

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

    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[270,60,270])
    localcourierisometry = localstates[2]|>FockMap
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[3]|>FockMap

    # filledseedsA1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteAregionfock1,shiftedsiteAregionfock1]), grouping=[10,6])[2]
    # filledseedsA2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteAregionfock2,shiftedsiteAregionfock2]), grouping=[10,6])[2]
    # filledseedsA3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteAregionfock3,shiftedsiteAregionfock3]), grouping=[10,6])[2]
    # filledseedsB1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteBregionfock1,shiftedsiteBregionfock1]), grouping=[10,6])[2]
    # filledseedsB2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteBregionfock2,shiftedsiteBregionfock2]), grouping=[10,6])[2]
    # filledseedsB3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteBregionfock3,shiftedsiteBregionfock3]), grouping=[10,6])[2]

    # filledseeds = (filledseedsA1+filledseedsA2+filledseedsA3+filledseedsB1+filledseedsB2+filledseedsB3)|>FockMap
    # localwannierfilled = localwannierization(localfilledisometry, filledseeds)

    # emptyseedsA1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteAregionfock1,shiftedsiteAregionfock1]), grouping=[10,6])[2]
    # emptyseedsA2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteAregionfock2,shiftedsiteAregionfock2]), grouping=[10,6])[2]
    # emptyseedsA3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteAregionfock3,shiftedsiteAregionfock3]), grouping=[10,6])[2]
    # emptyseedsB1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteBregionfock1,shiftedsiteBregionfock1]), grouping=[10,6])[2]
    # emptyseedsB2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteBregionfock2,shiftedsiteBregionfock2]), grouping=[10,6])[2]
    # emptyseedsB3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteBregionfock3,shiftedsiteBregionfock3]), grouping=[10,6])[2]

    # emptyseeds = (emptyseedsA1+emptyseedsA2+emptyseedsA3+emptyseedsB1+emptyseedsB2+emptyseedsB3)|>FockMap
    # localwannierempty = localwannierization(localemptyisometry, emptyseeds)

    courierseedsA1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteAregionfock1,shiftedsiteAregionfock1]), grouping=[90,10])[2]
    courierseedsA2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteAregionfock2,shiftedsiteAregionfock2]), grouping=[90,10])[2]
    courierseedsA3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteAregionfock3,shiftedsiteAregionfock3]), grouping=[90,10])[2]
    courierseedsB1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteBregionfock1,shiftedsiteBregionfock1]), grouping=[90,10])[2]
    courierseedsB2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteBregionfock2,shiftedsiteBregionfock2]), grouping=[90,10])[2]
    courierseedsB3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteBregionfock3,shiftedsiteBregionfock3]), grouping=[90,10])[2]

    courierseeds = (courierseedsA1+courierseedsA2+courierseedsA3+courierseedsB1+courierseedsB2+courierseedsB3)|>FockMap
    localwanniercourier = localwannierization(localcourierisometry, courierseeds)


    shiftedlocalwannierresults =  Dict(:localwanniercourier => localwanniercourier, :localwannierfilled => localfilledisometry,
                                        :localwannierempty => localemptyisometry, :courierseeds => courierseeds)
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

wanniercourierisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
wannierfilledisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
wannieremptyisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedempty[:,refunictcellfockempty])

@info "Computing local courier states..."
leftrestrict = fourier(wanniercourierisometry|>getoutspace, firsthexagonalregionfock) / (blockedcrystal|>vol|>sqrt)
rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
wanniercourierstate = leftrestrict' * wanniercourierisometry * rightrestrict

couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
display(couriercorrelations|>eigspech|>visualize)

filledcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
display(filledcorrelations|>eigspech|>visualize)

emptycorrelations = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry

couriercorrelationspectrum = couriercorrelations |> crystalspectrum

purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

firstgmeracorrelations = purifiedcouriercorrelations

firstgmeracrystalfock = firstgmeracorrelations|>getoutspace
firstgmeracrystal::Crystal = firstgmeracrystalfock|>getcrystal
firstgmeraspace::RealSpace = firstgmeracrystal|>getspace

noofflavourpermode = 10

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
localsiteB1correlations|>eigspec|>visualize

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[18,24,18])
localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[3]|>FockMap

(localfilledisometry*localfilledisometry')[siteAregionfock1,siteAregionfock1]|>eigspech|>visualize
(localcourierisometry*localcourierisometry')[siteAregionfock1,siteAregionfock1]|>eigspech|>visualize

# filledseedsA1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock1,siteAregionfock1]), grouping=[2,2])[2]
# filledseedsA2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock2,siteAregionfock2]), grouping=[2,2])[2]
# filledseedsA3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock3,siteAregionfock3]), grouping=[2,2])[2]
# filledseedsB1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteBregionfock1,siteBregionfock1]), grouping=[2,2])[2]
# filledseedsB2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteBregionfock2,siteBregionfock2]), grouping=[2,2])[2]
# filledseedsB3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteBregionfock3,siteBregionfock3]), grouping=[2,2])[2]

# filledseeds = (filledseedsA1+filledseedsA2+filledseedsA3+filledseedsB1+filledseedsB2+filledseedsB3)|>FockMap
# localwannierfilled = localwannierization(localfilledisometry, filledseeds)

# emptyseedsA1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteAregionfock1,siteAregionfock1]), grouping=[2,2])[2]
# emptyseedsA2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteAregionfock2,siteAregionfock2]), grouping=[2,2])[2]
# emptyseedsA3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteAregionfock3,siteAregionfock3]), grouping=[2,2])[2]
# emptyseedsB1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteBregionfock1,siteBregionfock1]), grouping=[2,2])[2]
# emptyseedsB2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteBregionfock2,siteBregionfock2]), grouping=[2,2])[2]
# emptyseedsB3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[siteBregionfock3,siteBregionfock3]), grouping=[2,2])[2]

# emptyseeds = (emptyseedsA1+emptyseedsA2+emptyseedsA3+emptyseedsB1+emptyseedsB2+emptyseedsB3)|>FockMap
# localwannierempty = localwannierization(localemptyisometry, emptyseeds)

courierseedsA1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteAregionfock1,siteAregionfock1]), grouping=[6,4])[2]
courierseedsA2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteAregionfock2,siteAregionfock2]), grouping=[6,4])[2]
courierseedsA3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteAregionfock3,siteAregionfock3]), grouping=[6,4])[2]
courierseedsB1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteBregionfock1,siteBregionfock1]), grouping=[6,4])[2]
courierseedsB2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteBregionfock2,siteBregionfock2]), grouping=[6,4])[2]
courierseedsB3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteBregionfock3,siteBregionfock3]), grouping=[6,4])[2]

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

    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[18,24,18])
    localcourierisometry = localstates[2]|>FockMap
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[3]|>FockMap

    # filledseedsA1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteAregionfock1,shiftedsiteAregionfock1]), grouping=[2,2])[2]
    # filledseedsA2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteAregionfock2,shiftedsiteAregionfock2]), grouping=[2,2])[2]
    # filledseedsA3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteAregionfock3,shiftedsiteAregionfock3]), grouping=[2,2])[2]
    # filledseedsB1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteBregionfock1,shiftedsiteBregionfock1]), grouping=[2,2])[2]
    # filledseedsB2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteBregionfock2,shiftedsiteBregionfock2]), grouping=[2,2])[2]
    # filledseedsB3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[shiftedsiteBregionfock3,shiftedsiteBregionfock3]), grouping=[2,2])[2]

    # filledseeds = (filledseedsA1+filledseedsA2+filledseedsA3+filledseedsB1+filledseedsB2+filledseedsB3)|>FockMap
    # localwannierfilled = localwannierization(localfilledisometry, filledseeds)

    # emptyseedsA1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteAregionfock1,shiftedsiteAregionfock1]), grouping=[2,2])[2]
    # emptyseedsA2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteAregionfock2,shiftedsiteAregionfock2]), grouping=[2,2])[2]
    # emptyseedsA3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteAregionfock3,shiftedsiteAregionfock3]), grouping=[2,2])[2]
    # emptyseedsB1 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteBregionfock1,shiftedsiteBregionfock1]), grouping=[2,2])[2]
    # emptyseedsB2 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteBregionfock2,shiftedsiteBregionfock2]), grouping=[2,2])[2]
    # emptyseedsB3 = getregionstates(localcorrelations=((localemptyisometry*localemptyisometry')[shiftedsiteBregionfock3,shiftedsiteBregionfock3]), grouping=[2,2])[2]

    # emptyseeds = (emptyseedsA1+emptyseedsA2+emptyseedsA3+emptyseedsB1+emptyseedsB2+emptyseedsB3)|>FockMap
    # localwannierempty = localwannierization(localemptyisometry, emptyseeds)

    courierseedsA1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteAregionfock1,shiftedsiteAregionfock1]), grouping=[6,4])[2]
    courierseedsA2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteAregionfock2,shiftedsiteAregionfock2]), grouping=[6,4])[2]
    courierseedsA3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteAregionfock3,shiftedsiteAregionfock3]), grouping=[6,4])[2]
    courierseedsB1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteBregionfock1,shiftedsiteBregionfock1]), grouping=[6,4])[2]
    courierseedsB2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteBregionfock2,shiftedsiteBregionfock2]), grouping=[6,4])[2]
    courierseedsB3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[shiftedsiteBregionfock3,shiftedsiteBregionfock3]), grouping=[6,4])[2]

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

wanniercourierisometry = globalwannierfunction(firstgmeracorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
wannierfilledisometry = globalwannierfunction(firstgmeracorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
wannieremptyisometry = globalwannierfunction(firstgmeracorrelations,extendedwannierizedempty[:,refunictcellfockempty])

@info "Computing local courier states..."
leftrestrict = fourier(wanniercourierisometry|>getoutspace,secondhexagonalregionfock) / (firstgmeracrystal|>vol|>sqrt)
rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
wanniercourierstate = leftrestrict' * wanniercourierisometry * rightrestrict

couriercorrelations = wanniercourierisometry' * firstgmeracorrelations * wanniercourierisometry
display(couriercorrelations|>eigspech|>visualize)

filledcorrelations = wannierfilledisometry' * firstgmeracorrelations * wannierfilledisometry
display(filledcorrelations|>eigspech|>visualize)

emptycorrelations = wannieremptyisometry' * firstgmeracorrelations * wannieremptyisometry

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
noofflavourpermode=4

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

localsiteA1correlations = regioncorrelations(secondgmeracorrelations,siteAregionfock1)
localsiteA2correlations = regioncorrelations(secondgmeracorrelations,siteAregionfock2)
localsiteA3correlations = regioncorrelations(secondgmeracorrelations,siteAregionfock3)
localsiteB1correlations = regioncorrelations(secondgmeracorrelations,siteBregionfock1)
localsiteB2correlations = regioncorrelations(secondgmeracorrelations,siteBregionfock2)
localsiteB3correlations = regioncorrelations(secondgmeracorrelations,siteBregionfock3)
localsiteA1correlations|>eigspec|>visualize

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[9,6,9])
localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[3]|>FockMap

# (localfilledisometry*localfilledisometry')[siteAregionfock1,siteAregionfock1]|>eigspech|>visualize
(localcourierisometry*localcourierisometry')[siteAregionfock1,siteAregionfock1]|>eigspech|>visualize


courierseedsA1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteAregionfock1,siteAregionfock1]), grouping=[3,1])[2]
courierseedsA2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteAregionfock2,siteAregionfock2]), grouping=[3,1])[2]
courierseedsA3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteAregionfock3,siteAregionfock3]), grouping=[3,1])[2]
courierseedsB1 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteBregionfock1,siteBregionfock1]), grouping=[3,1])[2]
courierseedsB2 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteBregionfock2,siteBregionfock2]), grouping=[3,1])[2]
courierseedsB3 = getregionstates(localcorrelations=((localcourierisometry*localcourierisometry')[siteBregionfock3,siteBregionfock3]), grouping=[3,1])[2]

courierseeds = (courierseedsA1+courierseedsA2+courierseedsA3+courierseedsB1+courierseedsB2+courierseedsB3)|>FockMap
localwanniercourier = localwannierization(localcourierisometry, courierseeds)

localwannierresults =  Dict(:localwanniercourier => localwanniercourier, :localwannierfilled => localfilledisometry,
                                :localwannierempty => localemptyisometry, :courierseeds => courierseeds)

wannierinfos =  Dict(thirdhexagonalregionfock=>localwannierresults)

wanniercourierisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localwanniercourier])
wannierfilledisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localwannierfilled])
wannieremptyisometry = globalwannierfunction(secondgmeracorrelations,wannierinfos[thirdhexagonalregionfock][:localwannierempty])
filledcorrelations = wannierfilledisometry' * secondgmeracorrelations * wannierfilledisometry
visualize(emptycorrelations|>eigspech)
emptycorrelations = wannieremptyisometry' * secondgmeracorrelations * wannieremptyisometry
couriercorrelations = wanniercourierisometry' * secondgmeracorrelations * wanniercourierisometry
couriercorrelations|>eigspech|>visualize

couriercorrelationspectrum = couriercorrelations |> crystalspectrum
purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap


coarsegrainedcrystal = purifiedcouriercorrelations|>getoutspace|>getcrystal

firstcenter = [0,0] ∈ blockedspace
firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=coarsegrainedcrystal, center=firstcenter, metricspace=blockedspace)
firsthexagonalregionfock = quantize(firsthexagonalregion,1)
localcorrelations = regioncorrelations(purifiedcouriercorrelations,firsthexagonalregionfock)
localspectrum = localcorrelations|>eigspec
localcorrelations[[m for m in firsthexagonalregionfock][5],[m for m in firsthexagonalregionfock][5]]|>eigspech|>visualize
localcorrelations[[m for m in firsthexagonalregionfock][3],[m for m in firsthexagonalregionfock][3]]|>eigspech|>geteigenvalues