using Zipper
using LinearAlgebra,Plots
BLAS.set_num_threads(8) 

setmaxthreads(setmaxthreads(Threads.nthreads()))
usecrystaldensemap()
plotlyjs()

power = 5
onsitepotential = 0
nnhopping = 0.1im
systemsize=2^power
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystalfock = correlations|>getoutspace

scaling = 2

blockedcorrelations,blockedH,blocker = blocking(correlations,H,scaling)
noofmodesinlocalreg = (6*(scaling^2))|>Int
noofdistillablemodes = ((1/4)*noofmodesinlocalreg)|>Int
noofcouriermodesinfirststep = noofmodesinlocalreg - noofdistillablemodes
noofcouriermodes = noofcouriermodesinfirststep
noofcouriermodesinsecondstep = noofcouriermodesinfirststep - noofdistillablemodes
noofcouriermodesinthirdstep = noofcouriermodesinsecondstep - noofdistillablemodes
noofflavourpermode = 1

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

firstinnerhexagonalregion = rankedandgroupoffsets(firsthexagonalregion,10)
firstrgshiftedinnerhexagonalregion1 = firstinnerhexagonalregion.+rgshiftedcenter1
firstrgshiftedinnerhexagonalregion2 = firstinnerhexagonalregion.+rgshiftedcenter2

localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
localspectrum = localcorrelations|>eigspec
display(localspectrum|>visualize)

nooffrozenmodes = length(firsthexagonalregionfock|>orderedmodes)-noofcouriermodes
nooffilledmodes = div(nooffrozenmodes,2)
noofemptymodes = nooffilledmodes
@info("no of filled modes = ",nooffilledmodes)
@info("no of courier modes = ",noofcouriermodes)
@info("no of empty modes = ",noofemptymodes)

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3,18, 3])
localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[3]|>FockMap
courierproj = localcourierisometry*localcourierisometry'
filledproj = localfilledisometry*localfilledisometry'
emptyproj = localemptyisometry*localemptyisometry'

siteAregion1 = intersect(intersect(firsthexagonalregion,firstrgshiftedhexagonalregion1),firstrgshiftedhexagonalregion2)
siteAregionfock1 = quantize(siteAregion1,noofflavourpermode)
siteAregion2 = c3*siteAregion1
siteAregionfock2 = quantize(siteAregion2,noofflavourpermode)
siteAregion3 = c3*siteAregion2
siteAregionfock3 = quantize(siteAregion3,noofflavourpermode)

siteBregion1 = (c6)*siteAregion1
siteBregionfock1 = quantize(siteBregion1,noofflavourpermode)
siteBregion2 = c3*siteBregion1
siteBregionfock2 = quantize(siteBregion2,noofflavourpermode)
siteBregion3 = c3*siteBregion2
siteBregionfock3 = quantize(siteBregion3,noofflavourpermode)


# nooflocalcourierseeds = 16
# localAcorrelations11 = columns(rows(courierproj,siteAregionfock11),siteAregionfock11)
# localAstates11 = getregionstates(localcorrelations=localAcorrelations11, grouping=[nooflocalcourierseeds])
# localAseeds11 = localAstates11[1]
# localAcorrelations12 = columns(rows(courierproj,siteAregionfock12),siteAregionfock12)
# localAstates12 = getregionstates(localcorrelations=localAcorrelations12, grouping=[nooflocalcourierseeds])
# localAseeds12 = localAstates12[1]
# localAcorrelations1c = columns(rows(courierproj,siteAregionfock1c),siteAregionfock1c)
# localAstates1c = getregionstates(localcorrelations=localAcorrelations1c, grouping=[5,11])
# localAseeds1c = localAstates1c[2]
# localAcorrelations1u = columns(rows(courierproj,siteAregionfock1u),siteAregionfock1u)
# localAstates1u = getregionstates(localcorrelations=localAcorrelations1u, grouping=[11,5])
# localAseeds1u = localAstates1u[2]

# localAcorrelations21 = columns(rows(courierproj,siteAregionfock21),siteAregionfock21)
# localAstates21 = getregionstates(localcorrelations=localAcorrelations21, grouping=[nooflocalcourierseeds])
# localAseeds21 = localAstates21[1]
# localAcorrelations22 = columns(rows(courierproj,siteAregionfock22),siteAregionfock22)
# localAstates22 = getregionstates(localcorrelations=localAcorrelations22, grouping=[nooflocalcourierseeds])
# localAseeds22 = localAstates22[1]
# localAcorrelations2c = columns(rows(courierproj,siteAregionfock2c),siteAregionfock2c)
# localAstates2c = getregionstates(localcorrelations=localAcorrelations2c, grouping=[5,11])
# localAseeds2c = localAstates2c[2]
# localAcorrelations2u = columns(rows(courierproj,siteAregionfock2u),siteAregionfock2u)
# localAstates2u = getregionstates(localcorrelations=localAcorrelations2u, grouping=[11,5])
# localAseeds2u = localAstates2u[2]

# localAcorrelations31 = columns(rows(courierproj,siteAregionfock31),siteAregionfock31)
# localAstates31 = getregionstates(localcorrelations=localAcorrelations31, grouping=[nooflocalcourierseeds])
# localAseeds31 = localAstates31[1]
# localAcorrelations32 = columns(rows(courierproj,siteAregionfock32),siteAregionfock32)
# localAstates32 = getregionstates(localcorrelations=localAcorrelations32, grouping=[nooflocalcourierseeds])
# localAseeds32 = localAstates32[1]
# localAcorrelations3c = columns(rows(courierproj,siteAregionfock3c),siteAregionfock3c)
# localAstates3c = getregionstates(localcorrelations=localAcorrelations3c, grouping=[5,11])
# localAseeds3c = localAstates3c[2]
# localAcorrelations3u = columns(rows(courierproj,siteAregionfock3u),siteAregionfock3u)
# localAstates3u = getregionstates(localcorrelations=localAcorrelations3u, grouping=[11,5])
# localAseeds3u = localAstates3u[2]

# allAseedsstate = localAseeds11+localAseeds12+localAseeds1c+localAseeds1u+localAseeds21+localAseeds22+localAseeds2c+localAseeds2u+localAseeds31+localAseeds32+localAseeds3c+localAseeds3u

# localBcorrelations11 = columns(rows(courierproj,siteBregionfock11),siteBregionfock11)
# localBstates11 = getregionstates(localcorrelations=localBcorrelations11, grouping=[nooflocalcourierseeds])
# localBseeds11 = localBstates11[1]
# localBcorrelations12 = columns(rows(courierproj,siteBregionfock12),siteBregionfock12)
# localBstates12 = getregionstates(localcorrelations=localBcorrelations12, grouping=[nooflocalcourierseeds])
# localBseeds12 = localBstates12[1]
# localBcorrelations1c = columns(rows(courierproj,siteBregionfock1c),siteBregionfock1c)
# localBstates1c = getregionstates(localcorrelations=localBcorrelations1c, grouping=[5,11])
# localBseeds1c = localBstates1c[2]
# localBcorrelations1u = columns(rows(courierproj,siteBregionfock1u),siteBregionfock1u)
# localBstates1u = getregionstates(localcorrelations=localBcorrelations1u, grouping=[11,5])
# localBseeds1u = localBstates1u[2]

# localBcorrelations21 = columns(rows(courierproj,siteBregionfock21),siteBregionfock21)
# localBstates21 = getregionstates(localcorrelations=localBcorrelations21, grouping=[nooflocalcourierseeds])
# localBseeds21 = localBstates21[1]
# localBcorrelations22 = columns(rows(courierproj,siteBregionfock22),siteBregionfock22)
# localBstates22 = getregionstates(localcorrelations=localBcorrelations22, grouping=[nooflocalcourierseeds])
# localBseeds22 = localBstates22[1]
# localBcorrelations2c = columns(rows(courierproj,siteBregionfock2c),siteBregionfock2c)
# localBstates2c = getregionstates(localcorrelations=localBcorrelations2c, grouping=[5,11])
# localBseeds2c = localBstates2c[2]
# localBcorrelations2u = columns(rows(courierproj,siteBregionfock2u),siteBregionfock2u)
# localBstates2u = getregionstates(localcorrelations=localBcorrelations2u, grouping=[11,5])
# localBseeds2u = localBstates2u[2]

# localBcorrelations31 = columns(rows(courierproj,siteBregionfock31),siteBregionfock31)
# localBstates31 = getregionstates(localcorrelations=localBcorrelations31, grouping=[nooflocalcourierseeds])
# localBseeds31 = localBstates31[1]
# localBcorrelations32 = columns(rows(courierproj,siteBregionfock32),siteBregionfock32)
# localBstates32 = getregionstates(localcorrelations=localBcorrelations32, grouping=[nooflocalcourierseeds])
# localBseeds32 = localBstates32[1]
# localBcorrelations3c = columns(rows(courierproj,siteBregionfock3c),siteBregionfock3c)
# localBstates3c = getregionstates(localcorrelations=localBcorrelations3c, grouping=[5,11])
# localBseeds3c = localBstates3c[2]
# localBcorrelations3u = columns(rows(courierproj,siteBregionfock3u),siteBregionfock3u)
# localBstates3u = getregionstates(localcorrelations=localBcorrelations3u, grouping=[11,5])
# localBseeds3u = localBstates3u[2]

# allBseedsstate = localBseeds11+localBseeds12+localBseeds1c+localBseeds1u+localBseeds21+localBseeds22+localBseeds2c+localBseeds2u+localBseeds31+localBseeds32+localBseeds3c+localBseeds3u


# allseedsstates = (allAseedsstate+allBseedsstate)
# courierseeds = allseedsstates|>FockMap

siteAregionfock = siteAregionfock1withoutu+siteAregionfock2withoutu+siteAregionfock3withoutu
siteBregionfock = siteBregionfock1withoutu+siteBregionfock2withoutu+siteBregionfock3withoutu
totalregionfock = siteAregionfock+siteBregionfock 

# test = Subset([m for m in totalregionfock|>orderedmodes])
svd(localcourierisometry[totalregionfock,:])

localwanniercourier,minsvd = localwannierization(localcourierisometry, courierseeds)

(localwanniercourier'*localcorrelations4*localwanniercourier)|>eigspec|>visualize

evalspair4 = [[key,norm(val)] for (key,val) in localspectrum4|>geteigenvalues]


evecs4 =  localspectrum4|>geteigenvectors

localstates4 = getregionstates(localcorrelations=localcorrelations4, grouping=[144, 96, 144])
localcourierisometry4 = localstates4[2]|>FockMap

truncatedevalspair4 = [m for m in localcourierisometry4|>getinspace]

eecontour = []
for rmode in firsthexagonalregionfock4
    save = []
    for (imode,eval) in evalspair4 
        if 0<round(eval,digits=8)<1
            append!(save,-norm(evecs4[rmode,imode])^2*(eval*log(eval)+(1-eval)*log(1-eval)) )
        end
    end
    append!(eecontour,tuple([rmode,sum(save)]))
end
eecontour
sortedeecontour = sort(eecontour,by=x->x[2])
Subset([(m|>getattr(:r))+(m|>getattr(:b)) for (m,_) in sortedeecontour[end-287:end]])|>visualize

gmera1firststepdata = gmerafirststepbycount(blockedcorrelations,blockedH,18,1)
gmera1firststepdata[:emptyH]|>getinspace|>getcrystal|>getunitcell|>visualize
gmera1secondstepdata4 = gmerasecondstepbycount(gmera1firststepdata4[:couriercorrelations],gmera1firststepdata4[:courierH],48,(96/6)|>Int)
gmera1thirdstepdata4 = gmerathirdstepbycount(gmera1secondstepdata4[:couriercorrelations],gmera1secondstepdata4[:courierH],24,(48/6)|>Int)
