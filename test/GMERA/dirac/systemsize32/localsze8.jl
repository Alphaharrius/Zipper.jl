using Zipper
using LinearAlgebra,Plots

setmaxthreads(setmaxthreads(Threads.nthreads()))
usecrystaldensemap()
plotlyjs()

power = 5
onsitepotential = 0
nnhopping = 0
systemsize=2^power
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystalfock = correlations|>getoutspace

scaling = 8

blockedcorrelations,blockedH,blocker = blocking(correlations,H,scaling)
noofflavourpermode = 1

blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal
blockedspace::RealSpace = blockedcrystal|>getspace
    
refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

firstcenter = [0,0] ∈ blockedspace
firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=firstcenter, metricspace=blockedspace)
firsthexagonalregionfock = quantize(firsthexagonalregion,noofflavourpermode)

firstinnerhexagonalregion = rankedandgroupoffsets(firsthexagonalregion,29)
(firsthexagonalregion-firstinnerhexagonalregion)|>visualize
# offstlist = [offset for offset in firstinnerhexagonalregion]
# firstinnerhexagonalregionhalf1 = Subset([offstlist[1],c3*offstlist[1],(c3)^2*offstlist[1]])
# firstinnerhexagonalregionhalf2 = Subset([c6*offstlist[1],c6*c3*offstlist[1],c6*(c3)^2*offstlist[1]])
firstinnerhexagonalregionfock = quantize(firstinnerhexagonalregion,1)
# firstinnerhexagonalregionhalf1fock = quantize(firstinnerhexagonalregionhalf1,1)
# firstinnerhexagonalregionhalf2fock = quantize(firstinnerhexagonalregionhalf2,1)

firstboundaryhexagonalregion = firsthexagonalregion-firstinnerhexagonalregion
firstboundaryhexagonalregionfock = quantize(firstboundaryhexagonalregion,1)

localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
localspectrum = localcorrelations|>eigspec
display(localspectrum|>visualize)
mevalpairs = [(m,norm(eval)) for (m,eval) in localspectrum|>geteigenvalues]
sortedmevalpairs = sort(mevalpairs,by=x->x[2])
sortedmevalpairs[1:144]
sortedmevalpairs[end-143:end]
chosenfrozenmevalpairs = sortedmevalpairs[144:end-145]
localevectors = localspectrum|>geteigenvectors
eecontour = []
for rsm in localcorrelations|>getinspace|>orderedmodes
    eesave = []
    for (m,eval) in chosenfrozenmevalpairs
        if 0<round(eval,digits=10)<1
            append!(eesave,norm(localevectors[rsm,m])^2*(-eval*log(eval)-(1-eval)*log(1-eval)))
        end
    end
    append!(eecontour,tuple([rsm,sum(eesave)]))
end
sortedeecontour = (sort(eecontour,by=x->x[2]))
firstboundaryhexagonalregionfock = quantize(Subset([(m|>getattr(:b))+(m|>getattr(:r)) for (m,_) in sortedeecontour[1:216]]),1)
# sortedeecontour[1:193]
Subset([(m|>getattr(:b))+(m|>getattr(:r)) for (m,_) in sortedeecontour[1:180]])|>visualize

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[108,168,108])
localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[3]|>FockMap
localfrozenisometry = (localstates[1]+localstates[3])|>FockMap

identity = idmap(localcorrelations|>getoutspace, localcorrelations|>getinspace)
localfrozenseeds = identity[:,firstboundaryhexagonalregionfock]
U,svdvals,Vd = svd(localfrozenseeds'*localfilledisometry)

localrot = U*Vd
localwannierfilled = localfilledisometry*localrot'
rotlocalcorrelations = localwannierfilled'*localcorrelations*localwannierfilled
inspace = localfrozenseeds|>getinspace
diff = localfrozenseeds-localwannierfilled
tr(diff'*diff)

sortedpair = sort([tuple(norm((rotlocalcorrelations)[m,m]|>rep),m) for m in rotlocalcorrelations|>getinspace|>orderedmodes],by=x->x[1])
Subset([((m|>getattr(:b))+(m|>getattr(:r))) for (_,m) in sortedpair[1:48]])|>visualize


U,svdvals,Vd = svd(localfrozenseeds'*localcourierisometry)
[val for val in svdvals]
localrot = U*Vd
localemptyisometry
localwannierempty = localemptyisometry*localrot'
rotlocalcorrelations = localemptyisometry'*localcorrelations*localemptyisometry
rotlocalcorrelations|>eigspec|>visualize

sortedpair = sort([tuple(norm((rotlocalcorrelations)[m,m]|>rep),m) for m in rotlocalcorrelations|>getinspace|>orderedmodes],by=x->x[1])
Subset([((m|>getattr(:b))+(m|>getattr(:r))) for (_,m) in sortedpair[1:144]])|>visualize


identity = idmap(localcorrelations|>getoutspace, localcorrelations|>getinspace)

localfilledseeds = identity[:,firstinnerhexagonalregionhalf1fock]
localemptyseeds = identity[:,firstinnerhexagonalregionhalf2fock]
localcourierseeds = identity[:,firstboundaryhexagonalregionfock]
Ufilled,svdvals,Vdfilled = svd(localfilledseeds'*localfilledisometry)
Uempty,svdvals,Vdempty = svd(localemptyseeds'*localemptyisometry)

localwanniercourier,minsvd = localwannierization(localcourierisometry, localcourierseeds)