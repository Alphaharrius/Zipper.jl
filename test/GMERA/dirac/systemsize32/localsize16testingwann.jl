using Pkg
Pkg.resolve()
Pkg.Registry.update()
Pkg.instantiate()
Pkg.activate("../../../../")
Pkg.resolve()

using Zipper
using LinearAlgebra,Plots

setmaxthreads(setmaxthreads(Threads.nthreads()))
usecrystaldensemap()
# plotlyjs()

power = 5
onsitepotential = 0
nnhopping = 0
systemsize=2^power
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystalfock = correlations|>getoutspace

scaling = 16
fiodir("../../../../../../../storage/data/slwongag/GMERA/dirac/systemsize$systemsize/localsize$scaling/testingwann")

blockedcorrelations,blockedH,blocker = blocking(correlations,H,scaling)
noofflavourpermode=1

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

wannierresult = Dict()
for nooffilledmodes in range(192,528,57)
    nooffilledmodes = nooffilledmodes|>Int
    noofemptymodes = nooffilledmodes
    nooffrozenmodes = nooffilledmodes+noofemptymodes
    noofcouriermodes = length(firsthexagonalregionfock|>orderedmodes)-nooffrozenmodes
    @info("no of filled modes = ",nooffilledmodes)
    @info("no of courier modes = ",noofcouriermodes)
    @info("no of empty modes = ",noofemptymodes)

    localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
    localspectrum = localcorrelations|>eigspec

    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofcouriermodes, noofemptymodes])
    localcourierisometry = localstates[2]|>FockMap
    courierproj = localcourierisometry*localcourierisometry'

    nooflocalcourierseeds = div(noofcouriermodes,6)
    reduandancy = length((siteAregionfock1|>orderedmodes))-nooflocalcourierseeds

    localAcorrelations1 = columns(rows(courierproj,siteAregionfock1),siteAregionfock1)
    localAcorrelations1|>eigspec|>visualize
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
    wannierresult[noofcouriermodes] = minsvd
    @info("no of courier seeds = ",noofcouriermodes)
    @info("minsvd = ",minsvd)
end
fiosave(wannierresultsecond,name="wannierresultsecond")
