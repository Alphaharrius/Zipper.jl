using Zipper
using LinearAlgebra
using Plots

setmaxthreads(8)

systemsize=32
correlations,H = generatesystem(  -0.3+0im, 0.3+0im,-1 + 0im,0im,systemsize)
crystalfock = correlations|>getoutspace

scale = Scale([8 0; 0 8], crystalfock|>getcrystal|>getspace)
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
sort([(key,norm(value)) for (key,value) in localspectrum|>geteigenvalues],by=x->x[2])[1:147]

localsiteA1correlations = regioncorrelations(blockedcorrelations,siteAregionfock1)
localsiteA2correlations = regioncorrelations(blockedcorrelations,siteAregionfock2)
localsiteA3correlations = regioncorrelations(blockedcorrelations,siteAregionfock3)
localsiteB1correlations = regioncorrelations(blockedcorrelations,siteBregionfock1)
localsiteB2correlations = regioncorrelations(blockedcorrelations,siteBregionfock2)
localsiteB3correlations = regioncorrelations(blockedcorrelations,siteBregionfock3)

localsiteA1correlations|>eigspech|>visualize


# localA1states = getregionstates(localcorrelations=localsiteA1correlations, grouping=[64])
# localA2states = getregionstates(localcorrelations=localsiteA2correlations, grouping=[64])
# localA3states = getregionstates(localcorrelations=localsiteA3correlations, grouping=[64])
# localB1states = getregionstates(localcorrelations=localsiteB1correlations, grouping=[64])
# localB2states = getregionstates(localcorrelations=localsiteB2correlations, grouping=[64])
# localB3states = getregionstates(localcorrelations=localsiteB3correlations, grouping=[64])

# localcourierA1isometry = localA1states[1]|>FockMap
# localcourierA2isometry = localA2states[1]|>FockMap
# localcourierA3isometry = localA3states[1]|>FockMap
# localcourierB1isometry = localB1states[1]|>FockMap
# localcourierB2isometry = localB2states[1]|>FockMap
# localcourierB3isometry = localB3states[1]|>FockMap

# localfilledA1isometry = localA1states[1]|>FockMap
# localfilledA2isometry = localA2states[1]|>FockMap
# localfilledA3isometry = localA3states[1]|>FockMap
# localfilledB1isometry = localB1states[1]|>FockMap
# localfilledB2isometry = localB2states[1]|>FockMap
# localfilledB3isometry = localB3states[1]|>FockMap

# localemptyA1isometry = localA1states[3]|>FockMap
# localemptyA2isometry = localA2states[3]|>FockMap
# localemptyA3isometry = localA3states[3]|>FockMap
# localemptyB1isometry = localB1states[2]|>FockMap
# localemptyB2isometry = localB2states[2]|>FockMap
# localemptyB3isometry = localB3states[2]|>FockMap

# courierseeds = localcourierA1isometry+localcourierA2isometry+localcourierA3isometry+localcourierB1isometry+localcourierB2isometry+localcourierB3isometry
# filledseeds = localfilledA1isometry+localfilledA2isometry+localfilledA3isometry
# emptyseeds = localemptyB1isometry+localemptyB2isometry+localemptyB3isometry

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[150,84,150])
localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[3]|>FockMap

# (localcourierisometry*localcourierisometry')[siteAregionfock1,siteAregionfock1]|>eigspech|>visualize

(localfilledisometry*localfilledisometry')[siteBregionfock1,siteBregionfock1]|>eigspech|>visualize

# (localfilledisometry*localfilledisometry')[firsthexagonalregionfock,firsthexagonalregionfock]|>eigspech|>visualize


subhexagonalregion = rankedandgroupoffsets(firsthexagonalregion,16)
subhexagonalregionfock = quantize(subhexagonalregion,1)
# visualize(subhexagonalregion)

(localfilledisometry*localfilledisometry')[subhexagonalregionfock,subhexagonalregionfock]|>eigspech|>visualize

# (localcourierisometry*localcourierisometry')[subhexagonalregionfock,subhexagonalregionfock]|>eigspech|>visualize

filledseedsA1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock1,siteAregionfock1]), grouping=[49,15])[2]
filledseedsA2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock2,siteAregionfock2]), grouping=[49,15])[2]
filledseedsA3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteAregionfock3,siteAregionfock3]), grouping=[49,15])[2]

filledseedsB1 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteBregionfock1,siteBregionfock1]), grouping=[42,22])[2]
filledseedsB2 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteBregionfock2,siteBregionfock2]), grouping=[42,22])[2]
filledseedsB3 = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[siteBregionfock3,siteBregionfock3]), grouping=[42,22])[2]

filledseedscenter = getregionstates(localcorrelations=((localfilledisometry*localfilledisometry')[subhexagonalregionfock,subhexagonalregionfock]), grouping=[57,39])[2]

filledseeds = (filledseedscenter)|>FockMap
filledseeds|>getoutspace

localwannierization(localfilledisometry[filledseeds|>getoutspace,:], filledseeds)

localwanniercourier = localwannierization(localcourierisometry, courierseeds)
scatter(sort([norm(ele) for ele in (localwanniercourier'*localcorrelations*localwanniercourier).rep|>diag]))
localwannierfilled = localwannierization(localfilledisometry[filledseeds|>getoutspace|>RegionFock,:], filledseeds)
localwannierempty = localwannierization(localemptyisometry[emptyseeds|>getoutspace|>RegionFock,:], emptyseeds)